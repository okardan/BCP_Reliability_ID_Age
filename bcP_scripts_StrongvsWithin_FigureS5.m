% Kardan, O., Kaplan, S., ..., & Rosenberg, M.D. (2022) 
% "Resting-state functional connectivity identifies individuals and
% predicts age in 8-to-26-month-olds" Developmental Cognitive Neuroscience

% Script to compare age prediction power of 6017 strong between-net
% edges with the 6017 within-net edges and make Figure S5

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

   
%%   determine strong between-net edges
load('BCP_averageFC_babynetassignment.mat');
outnet_conns = GG3.*outnet_edges;
% figure;imagesc(outnet_conns);
outconn = GG3(outnet_edges==1);
[sortedoutconn,I] = sort(outconn,'descend');
top12kidx = I(1:6017);
%%
Session =[]; Dat =[]; Age =[]; Pred_Age_totFC =[]; Pred_Age_strbetFC=[]; Pred_Age_witFC=[];
FDs =[]; aQC =[]; fQC =[]; Volumes=[];TRs=[];FCval=[];FCvalsd=[];FCvaltstat=[];
 for zx = 1:500   
    all_braincondat_fullnet=[]; all_ages=[]; all_subid=[]; all_Vol=[]; all_fqc =[];
    all_aqc=[]; all_fd =[]; all_TR=[]; all_ttest =[];
    all_braincondat_Strongbetwnet =[]; all_braincondat_withnet =[];

    % ##########
%     rng 'default'
    rng 'shuffle'
    % ##########
    subs =[];
    for i =1:length(goodGMDec)
        subs = [subs; string(goodGMDec(i).name(5:10))];
    end
%     qs = randperm(12000,6017);
     qs = 1:6017;
    uniqs = unique(subs);
    for j=1:length(uniqs)
        idx = find(subs==uniqs(j));
        ii=randi(length(idx));
       sub_data = goodGMDec(idx(ii));
      
       try
       tseries = [sub_data.GordonMatAP(1:333,1:300) sub_data.GordonMatPA(1:333,1:300)];
       catch
%            disp(sub_data.name); 
%            a= size(sub_data.GordonMatAP,2)
%            b = size(sub_data.GordonMatPA,2)
       tseries = [sub_data.GordonMat(1:333,1:600)];  
       end
       totalcondat1 =  corr(tseries(All_parcels,1:600)');
%        totalconndat = totalcondat1(All_edges);
       outconndat = totalcondat1(find(outnet_edges==1)); 
             
       inconndat =[];
       for k=1:11
           temp_condat = corr(tseries(innet_parcels{k},1:600)');
           inconndat = [inconndat temp_condat(innet_edges{k})'];
       end
%         all_braincondat_fullnet = [all_braincondat_fullnet;totalconndat'];

 fcs = outconndat(top12kidx(qs))';
        all_braincondat_Strongbetwnet =[all_braincondat_Strongbetwnet; fcs];
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
        [h,p,ci,stats] = ttest2(fcs',inconndat');
       all_ttest = [all_ttest; stats.tstat];
    end
 
   xylabel = 'Age (months)';

        Y = all_ages;
        
       n=size(Y,1);
       FCval = [FCval; [mean(mean(all_braincondat_Strongbetwnet,2)) mean(mean(all_braincondat_withnet,2))]];
       FCvalsd = [FCvalsd; [std(mean(all_braincondat_Strongbetwnet,2)) std(mean(all_braincondat_withnet,2))]];
       FCvaltstat = [FCvaltstat;all_ttest];
 %%%
        X = [atanh(all_braincondat_withnet)]; 
        MdlLin2 = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        YHat2 = kfoldPredict(MdlLin2); n=size(Y,1);


        X = [atanh(all_braincondat_Strongbetwnet)]; 
        MdlLin3 = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        YHat3 = kfoldPredict(MdlLin3); n=size(Y,1);


        Session =[Session; all_subid]; 
        Dat =[Dat; zx.*ones(n,1)];
        Age =[Age; all_ages]; 
        Volumes =[Volumes; all_Vol];
%         Pred_Age_totFC = [Pred_Age_totFC; YHat1]; 
        Pred_Age_strbetFC = [Pred_Age_strbetFC; YHat3]; 
        Pred_Age_witFC = [Pred_Age_witFC; YHat2];
        FDs =[FDs; all_fd]; 
        aQC =[aQC; all_aqc]; 
        fQC =[fQC; all_fqc];
        TRs =[TRs; all_TR];
%  %%%   
    zx
 end
T = table(Session,Dat,Age,Volumes,...
    Pred_Age_strbetFC,Pred_Age_witFC,...
    FDs, aQC, fQC, TRs,FCvaltstat);
writetable(T,'AgePredStacked_withORstrongest_1_500.csv');
T2 = table(FCval,FCvalsd);
writetable(T2,'AgeStats_withORstrongest_1_500.csv');
%%
% clear all
T = readtable('AgePredStacked_withORstrongest_1_500.csv');
dist_stats =[]; 
for k=1:500

datsamp = T(T.Dat==k,:);
datsamp.volsq = (datsamp.Volumes-nanmean(datsamp.Volumes)).^2;
Y = datsamp.Age;
n=size(Y,1);
% strong between net
YHat = datsamp.Pred_Age_strbetFC;
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq1 = 1 -(mseLin./errvar);
[r1 p ] = corr(Y,YHat,'Rows','complete');
[r2 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');

% within net

YHat = datsamp.Pred_Age_witFC;
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq2 = 1 -(mseLin./errvar);
[r3 p ] = corr(Y,YHat,'Rows','complete');
[r4 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');

% regressions
glm1 = fitglm(datsamp,'Age ~ FDs +aQC +fQC +TRs');
glm2 = fitglm(datsamp,'Age ~ FDs +aQC +fQC +TRs+ Volumes + volsq');

glm4 = fitglm(datsamp,'Age ~ Pred_Age_strbetFC + FDs +aQC +fQC +TRs +Volumes +volsq');
glm5 = fitglm(datsamp,'Age ~ Pred_Age_witFC + FDs +aQC +fQC +TRs +Volumes +volsq');

deltaRsq2 = glm5.Rsquared.Adjusted - glm4.Rsquared.Adjusted;

dist_stats = [dist_stats; [r1 r3 rsq1 rsq2 r2 r4 ...
    glm1.Rsquared.Adjusted glm2.Rsquared.Adjusted, ...  %10,11
    glm4.Rsquared.Adjusted, ...
    glm5.Rsquared.Adjusted, ...
    deltaRsq2]];
k
end
%% Making Supplementary Figure S5
figure
subplot(3,3,4:6)
distributionPlot({dist_stats(:,5)},'color',{[.3 .9 .7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,6)},'color',{[.9 .7 .3]},'xyOri','flipped','histOri','right','showMM',0)
legend({'Strong positive 6017 between net FC Model','Withinn net FC Model'},'AutoUpdate','off','FontSize',16);
xlim([-1,1]);xlabel('Prediction R^2');ylabel('Count (normalized)');set(gca,'FontSize',16)

subplot(3,3,7:8)
distributionPlot({dist_stats(:,7)},'color',{[0.7 0.7 0.7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,8)},'color',{[0.20,0.75,0.23]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,9)},'color',{[.3 .9 .7]},'xyOri','flipped','histOri','right','showMM',0)
distributionPlot({dist_stats(:,10)},'color',{[.9 .7 .3]},'xyOri','flipped','histOri','right','showMM',0)
legend({'Nuisance only','Nuisance & Brain Volume','Nuisance & Brain Volume & 6017 between Net FC',...
    'Nuisance & Volume & within Net FC'},'AutoUpdate','off','FontSize',16);
xlim([0,1]);xlabel('Regression Models R^2');ylabel('Count (normalized)');set(gca,'FontSize',16)
p1 = length(find(dist_stats(:,5)>=dist_stats(:,6)))/500;text(0,.5,['P(\leq 0) = ',num2str(p1)],'FontSize',16);

subplot(3,3,9)
distributionPlot({dist_stats(:,11)},'color',{[.9 .7 .3]},'xyOri','flipped','histOri','right','showMM',0)
legend({['Difference between models']'},'AutoUpdate','off','FontSize',16);
xlim([-.3,.5]);xlabel('\Delta R^2');ylabel('Count (normalized)');set(gca,'FontSize',16,'XTick',[-.4:.2:.4])
p2 = length(find(dist_stats(:,11)<=0))/500; text(0,.5,['P(\leq 0) = ',num2str(p2)],'FontSize',16);
p1 = length(find(dist_stats(:,5)>=dist_stats(:,6)))/500;