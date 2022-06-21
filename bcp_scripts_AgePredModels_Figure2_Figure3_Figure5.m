% Kardan, O., Kaplan, S., ..., & Rosenberg, M.D. (2022) 
% "Resting-state functional connectivity identifies individuals and
% predicts age in 8-to-26-month-olds" Developmental Cognitive Neuroscience

% This script runs the SVR for predicting age using full connectomes or
% only within-net or only between-net edges with 10-fold cross validation
% and 500 resamples of the 170 sessions (with n = 110 unique subjects in each resample)
% and later make figures 2 and 3 from the paper

clear all
load('goodGMDec2.mat')
   Vol = bcp_volumes('BCP_Brain_Segmentation_Volume_in_mm3.csv');
   vol_name = strcat(string(Vol.ID),'_', string(Vol.visit));
    load('mask.mat'); numnets = 11;
    load('network.mat')
    load('networknames.mat')
    goodareas1 = [1:333]; subcort=0;% only corticals
    All_parcels = goodareas1(mask);   %ordered 
    All_edges = find(tril(ones(333),-1)==1);
    innet_parcels = cell(1,numnets);
    innet_edges = cell(1,numnets);
    for ntw=1:numnets
    innet_parcels{ntw} =  mask(network ==ntw);
    featsz = length(find(network ==ntw));
    innet_edges{ntw} = find(tril(ones(featsz),-1));
    end
    outnet_parcels = All_parcels;
    outnet_edges = zeros(333);
    for ntw1=1:numnets-1
        for ntw2=ntw1+1:numnets
            
            for ii=1:333
                for jj=1:333
                    if network(ii)==ntw1 & network(jj)==ntw2
                        outnet_edges(jj,ii) = 1;
                    end
                end
            end
        end
    end
  
Session =[]; Dat =[]; Age =[]; Pred_Age_totFC =[]; Pred_Age_betFC=[]; Pred_Age_witFC=[];
FDs =[]; aQC =[]; fQC =[]; Volumes=[];TRs=[];
 for zx = 1:500   
    all_braincondat_fullnet=[]; all_ages=[]; all_subid=[]; all_Vol=[]; all_fqc =[];
    all_aqc=[]; all_fd =[]; all_TR=[];
    all_braincondat_betwnet =[]; all_braincondat_withnet =[];

    % ##########
%     rng 'default'
    rng 'shuffle'
    % ##########
    subs =[];
    for i =1:length(goodGMDec)
        subs = [subs; string(goodGMDec(i).name(5:10))];
    end
    
    uniqs = unique(subs);
    for j=1:length(uniqs)
        idx = find(subs==uniqs(j));
        ii=randi(length(idx));
       sub_data = goodGMDec(idx(ii));
      
       try
       tseries = [sub_data.GordonMatAP(1:333,1:300) sub_data.GordonMatPA(1:333,1:300)];
       catch

       tseries = [sub_data.GordonMat(1:333,1:600)];  
       end
       totalcondat1 =  corr(tseries(All_parcels,1:600)');
       totalconndat = totalcondat1(All_edges);
       outconndat = totalcondat1(find(outnet_edges==1)); 
       inconndat =[];
       for k=1:11
           temp_condat = corr(tseries(innet_parcels{k},1:600)');
           inconndat = [inconndat temp_condat(innet_edges{k})'];
       end
        all_braincondat_fullnet = [all_braincondat_fullnet;totalconndat'];
        all_braincondat_betwnet =[all_braincondat_betwnet; outconndat'];
        all_braincondat_withnet =[all_braincondat_withnet; inconndat];
        all_ages = [all_ages; sub_data.Scanage];
        if ~isempty(find(vol_name == sub_data.name))
            vv = Vol.Brain_Segmentation_Volume_in_mm3(find(vol_name == sub_data.name));
        else
            vv = NaN;
        end
        all_Vol = [all_Vol; vv];    
        all_fqc =[all_fqc; sub_data.qcrate_func];
        all_aqc=[all_aqc; sub_data.qcrate_anat]; 
        all_fd =[all_fd; sub_data.FDmean];
        all_subid = [all_subid; string(sub_data.name)];
        all_TR = [all_TR; sub_data.TR];
       
    end
 
   xylabel = 'Age (months)';

        Y = all_ages;
        
        X = [atanh(all_braincondat_fullnet)];
        MdlLin1 = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        
        mseLin10 = kfoldLoss(MdlLin1,'Mode','individual');
        mseLin = mean(mseLin10);
        YHat1 = kfoldPredict(MdlLin1); n=size(Y,1);

        X = [atanh(all_braincondat_withnet)]; 
        MdlLin2 = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        YHat2 = kfoldPredict(MdlLin2); n=size(Y,1);


        X = [atanh(all_braincondat_betwnet)]; 
        MdlLin3 = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        YHat3 = kfoldPredict(MdlLin3); n=size(Y,1);

        Session =[Session; all_subid]; 
        Dat =[Dat; zx.*ones(n,1)];
        Age =[Age; all_ages]; 
        Volumes =[Volumes; all_Vol];
        Pred_Age_totFC = [Pred_Age_totFC; YHat1]; 
        Pred_Age_betFC = [Pred_Age_betFC; YHat3]; 
        Pred_Age_witFC = [Pred_Age_witFC; YHat2];
        FDs =[FDs; all_fd]; 
        aQC =[aQC; all_aqc]; 
        fQC =[fQC; all_fqc];
        TRs =[TRs; all_TR];
    
    zx
 end
T = table(Session,Dat,Age,Volumes,...
    Pred_Age_totFC,Pred_Age_betFC,Pred_Age_witFC,...
    FDs, aQC, fQC, TRs);
writetable(T,'AgePredStacked.csv');
%
Session =[]; Dat =[]; Age =[]; Pred_Age_totFC =[]; Pred_Age_betFC=[]; Pred_Age_witFC=[];
FDs =[]; aQC =[]; fQC =[]; Volumes=[];TRs=[]; shuff_inds=[];
 for zx = 1:500   
    all_braincondat_fullnet=[]; all_ages=[]; all_subid=[]; all_Vol=[]; all_fqc =[];
    all_aqc=[]; all_fd =[]; all_TR=[];
    all_braincondat_betwnet =[]; all_braincondat_withnet =[];

    % ##########
%     rng 'default'
    rng 'shuffle'
    % ##########
    subs =[];
    for i =1:length(goodGMDec)
        subs = [subs; string(goodGMDec(i).name(5:10))];
    end
    
    uniqs = unique(subs);
    for j=1:length(uniqs)
        idx = find(subs==uniqs(j));
        ii=randi(length(idx));
       sub_data = goodGMDec(idx(ii));
      
       try
       tseries = [sub_data.GordonMatAP(1:333,1:300) sub_data.GordonMatPA(1:333,1:300)];
       catch

       tseries = [sub_data.GordonMat(1:333,1:600)];  
       end
       totalcondat1 =  corr(tseries(All_parcels,1:600)');
       totalconndat = totalcondat1(All_edges);
       outconndat = totalcondat1(find(outnet_edges==1)); 
       inconndat =[];
       for k=1:11
           temp_condat = corr(tseries(innet_parcels{k},1:600)');
           inconndat = [inconndat temp_condat(innet_edges{k})'];
       end
        all_braincondat_fullnet = [all_braincondat_fullnet;totalconndat'];
        all_braincondat_betwnet =[all_braincondat_betwnet; outconndat'];
        all_braincondat_withnet =[all_braincondat_withnet; inconndat];
        all_ages = [all_ages; sub_data.Scanage];
        if ~isempty(find(vol_name == sub_data.name))
            vv = Vol.Brain_Segmentation_Volume_in_mm3(find(vol_name == sub_data.name));
        else
            vv = NaN;
        end
        all_Vol = [all_Vol; vv];    
        all_fqc =[all_fqc; sub_data.qcrate_func];
        all_aqc=[all_aqc; sub_data.qcrate_anat]; 
        all_fd =[all_fd; sub_data.FDmean];
        all_subid = [all_subid; string(sub_data.name)];
        all_TR = [all_TR; sub_data.TR];
       
    end
 
   xylabel = 'Age (months)';

        shuff_ids = randperm(length(uniqs));
        Y = all_ages(shuff_ids);
        
        X = [atanh(all_braincondat_fullnet)];
        MdlLin1 = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        
        mseLin10 = kfoldLoss(MdlLin1,'Mode','individual');
        mseLin = mean(mseLin10);
        YHat1 = kfoldPredict(MdlLin1); n=size(Y,1);

        X = [atanh(all_braincondat_withnet)]; 
        MdlLin2 = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        YHat2 = kfoldPredict(MdlLin2); n=size(Y,1);


        X = [atanh(all_braincondat_betwnet)]; 
        MdlLin3 = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        YHat3 = kfoldPredict(MdlLin3); n=size(Y,1);


        Session =[Session; all_subid]; 
        Dat =[Dat; zx.*ones(n,1)];
        Age =[Age; all_ages]; 
        Volumes =[Volumes; all_Vol];
        Pred_Age_totFC = [Pred_Age_totFC; YHat1]; 
        Pred_Age_betFC = [Pred_Age_betFC; YHat3]; 
        Pred_Age_witFC = [Pred_Age_witFC; YHat2];
        FDs =[FDs; all_fd]; 
        aQC =[aQC; all_aqc]; 
        fQC =[fQC; all_fqc];
        TRs =[TRs; all_TR];
        shuff_inds=[shuff_inds; shuff_ids'];
    
    zx
 end
T = table(Session,Dat,Age,Volumes,...
    Pred_Age_totFC,Pred_Age_betFC,Pred_Age_witFC,...
    FDs, aQC, fQC, TRs, shuff_inds);
writetable(T,'nullAgePredStacked.csv');
%% stats for visualization of the distributions
clear all
T = readtable('AgePredStacked.csv');
Tnull = readtable('nullAgePredStacked.csv');
dist_stats =[]; dist_nullstats=[];
for k=1:500

datsamp = T(T.Dat==k,:);
datsamp.volsq = (datsamp.Volumes-nanmean(datsamp.Volumes)).^2;
Y = datsamp.Age;
YHat = datsamp.Pred_Age_totFC;
n = length(~isnan(YHat));
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq1 = 1 -(mseLin./errvar);
[r1 p ] = corr(Y,YHat,'Rows','complete');
[r2 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');

% between net
YHat = datsamp.Pred_Age_betFC;
mseLinb = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq2 = 1 -(mseLinb./errvar);
[r3 p ] = corr(Y,YHat,'Rows','complete');
[r4 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');

% within net

YHat = datsamp.Pred_Age_witFC;
mseLinw = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq3 = 1 -(mseLinw./errvar);
[r5 p ] = corr(Y,YHat,'Rows','complete');
[r6 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');

% regressions
glm1 = fitglm(datsamp,'Age ~ FDs +aQC +fQC +TRs');
glm2 = fitglm(datsamp,'Age ~ FDs +aQC +fQC +TRs+ Volumes + volsq');
glm3 = fitglm(datsamp,'Age ~ Pred_Age_totFC + FDs +aQC +fQC +TRs +Volumes +volsq');
glm4 = fitglm(datsamp,'Age ~ Pred_Age_betFC + FDs +aQC +fQC +TRs +Volumes +volsq');
glm5 = fitglm(datsamp,'Age ~ Pred_Age_witFC + FDs +aQC +fQC +TRs +Volumes +volsq');
deltaRsq1 = glm3.Rsquared.Adjusted - glm2.Rsquared.Adjusted;
deltaRsq2 = glm5.Rsquared.Adjusted - glm4.Rsquared.Adjusted;

dist_stats = [dist_stats; [r1 r3 r5 rsq1 rsq2 rsq3 r2 r4 r6 ...
    glm1.Rsquared.Adjusted glm2.Rsquared.Adjusted ...  %10,11
    glm3.Rsquared.Adjusted glm4.Rsquared.Adjusted ...
    glm5.Rsquared.Adjusted, ...
    deltaRsq1  deltaRsq2,...
    mseLin mseLinb mseLinw]];
k
end

for k=1:500

datsamp = Tnull(Tnull.Dat==k,:);
datsamp.volsq = (datsamp.Volumes-nanmean(datsamp.Volumes)).^2;
Y = datsamp.Age(datsamp.shuff_inds);
YHat = datsamp.Pred_Age_totFC;
n = length(~isnan(YHat));
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq1 = 1 -(mseLin./errvar);
[r1 p ] = corr(Y,YHat,'Rows','complete');
[r2 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');

% between net
YHat = datsamp.Pred_Age_betFC;
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq2 = 1 -(mseLin./errvar);
[r3 p ] = corr(Y,YHat,'Rows','complete');
[r4 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');

% within net

YHat = datsamp.Pred_Age_witFC;
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq3 = 1 -(mseLin./errvar);
[r5 p ] = corr(Y,YHat,'Rows','complete');
[r6 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');

% regressions
glm1 = fitglm(datsamp,'Age ~ FDs +aQC +fQC +TRs');
glm2 = fitglm(datsamp,'Age ~ FDs +aQC +fQC +TRs+ Volumes + volsq');
glm3 = fitglm(datsamp,'Age ~ Pred_Age_totFC + FDs +aQC +fQC +TRs +Volumes +volsq');
glm4 = fitglm(datsamp,'Age ~ Pred_Age_betFC + FDs +aQC +fQC +TRs +Volumes +volsq');
glm5 = fitglm(datsamp,'Age ~ Pred_Age_witFC + FDs +aQC +fQC +TRs +Volumes +volsq');
deltaRsq1 = glm3.Rsquared.Adjusted - glm2.Rsquared.Adjusted;
deltaRsq2 = glm5.Rsquared.Adjusted - glm4.Rsquared.Adjusted;

dist_nullstats = [dist_nullstats; [r1 r3 r5 rsq1 rsq2 rsq3 r2 r4 r6 ... % r then rsq then part.r
    glm1.Rsquared.Adjusted glm2.Rsquared.Adjusted ...  %10 nui ,11 vol, 12 fc
    glm3.Rsquared.Adjusted glm4.Rsquared.Adjusted ...
    glm5.Rsquared.Adjusted,...
    deltaRsq1  deltaRsq2]]; % 15 16
k
end
%% comparing within and between network SVR models 
pwithbetw = length(find([dist_stats(:,5) - dist_stats(:,6)]>=0))/500;
%% Figure 2 Distributions of age-prediction model performance metrics
figure
subplot(3,1,1)
distributionPlot({dist_nullstats(:,1)},'color',{[0.7 0.7 0.7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,1)},'color',{[0.70,0.75,0.93]},'xyOri','flipped','histOri','right','showMM',0)
legend({'Null Model','Functional Connectivty Model'},'AutoUpdate','off','FontSize',16);
xlim([-1,1]);xlabel('Correlation between predicted and true age');ylabel('Count (normalized)');set(gca,'FontSize',16)

subplot(3,1,2)
distributionPlot({dist_nullstats(:,4)},'color',{[0.7 0.7 0.7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,4)},'color',{[0.70,0.75,0.93]},'xyOri','flipped','histOri','right','showMM',0)
legend({'Null Model','Functional Connectivty Model'},'AutoUpdate','off','FontSize',16);
xlim([-1,1]);xlabel('Prediction R^2');ylabel('Count (normalized)');set(gca,'FontSize',16)

subplot(3,1,3)
distributionPlot({dist_nullstats(:,7)},'color',{[0.7 0.7 0.7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,7)},'color',{[0.70,0.75,0.93]},'xyOri','flipped','histOri','right','showMM',0)
legend({'Null Model','Functional Connectivty Model'},'AutoUpdate','off','FontSize',16);
xlim([-1,1]);xlabel('Partial correlation');ylabel('Count (normalized)');set(gca,'FontSize',16)

%% Figure 3 Distributions of R2 for age-prediction models with nuisance variables 
figure
subplot(3,3,1:2)
distributionPlot({dist_stats(:,10)},'color',{[0.7 0.7 0.7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,11)},'color',{[0.20,0.75,0.23]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,12)},'color',{[0.70,0.75,0.93]},'xyOri','flipped','histOri','right','showMM',0)
legend({'Nuisance only','Nuisance & Volume','Nuisance & Volume & FC'},'AutoUpdate','off','FontSize',16);
xlim([0,1]);xlabel('Regression Models R^2');ylabel('Count (normalized)');set(gca,'FontSize',16)

subplot(3,3,3)
distributionPlot({dist_stats(:,15)},'color',{[0.70,0.75,0.93]},'xyOri','flipped','histOri','right','showMM',0)
legend({'(Nuisance & Volume) vs (Nuisance & Volume & FC)'},'AutoUpdate','off','FontSize',16);
xlim([-.5,.5]);xlabel('\Delta R^2');ylabel('Count (normalized)');set(gca,'FontSize',16,'XTick',[-.4:.2:.4])
p1 = 1 - (length(find(dist_stats(:,15)>0))/500); text(0,.5,['P(\leq 0) = ',num2str(p1)],'FontSize',16);
%% Figure 5 Performance of the rsFC age-prediction models using only within-net or between-net connections 

figure
subplot(2,3,1:3)
distributionPlot({dist_nullstats(:,5)},'color',{[.3 .9 .7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_nullstats(:,6)},'color',{[.9 .7 .3]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,5)},'color',{[.3 .9 .7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,6)},'color',{[.9 .7 .3]},'xyOri','flipped','histOri','right','showMM',0)
legend({'Null Model:Between Net FC',...
    'Null Model:Within Net FC',...
    'Between Net FC Model','Withinn Net FC Model'},'AutoUpdate','off','FontSize',16);
xlim([-1,1]);xlabel('Prediction R^2');ylabel('Count (normalized)');set(gca,'FontSize',16)

subplot(2,3,4:5)
distributionPlot({dist_stats(:,10)},'color',{[0.7 0.7 0.7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,11)},'color',{[0.20,0.75,0.23]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,13)},'color',{[.3 .9 .7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,14)},'color',{[.9 .7 .3]},'xyOri','flipped','histOri','right','showMM',0)
legend({'Nuisance only','Nuisance & Volume','Nuisance & Volume & between Net FC',...
    'Nuisance & Volume & within Net FC'},'AutoUpdate','off','FontSize',16);
xlim([0,1]);xlabel('Regression Models R^2');ylabel('Count (normalized)');set(gca,'FontSize',16)

subplot(2,3,6)
distributionPlot({dist_stats(:,16)},'color',{[.9 .7 .3]},'xyOri','flipped','histOri','right','showMM',0)
legend({['Difference between models']'},'AutoUpdate','off','FontSize',16);
xlim([-.5,.5]);xlabel('\Delta R^2');ylabel('Count (normalized)');set(gca,'FontSize',16,'XTick',[-.4:.2:.4])
p2 = 1 - (length(find(dist_stats(:,16)>0))/500); text(0,.5,['P(\leq 0) = ',num2str(p2)],'FontSize',16);

