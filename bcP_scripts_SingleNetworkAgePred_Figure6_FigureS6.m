% Kardan, O., Kaplan, S., ..., & Rosenberg, M.D. (2022) 
% "Resting-state functional connectivity identifies individuals and
% predicts age in 8-to-26-month-olds" Developmental Cognitive Neuroscience

% Script for individual network models predciting age (Figure 6 and Figure
% S6)

clear all
load('goodGMDec2.mat')
   Vol = bcp_volumes('CrossVal_Dec_BCP\BCP_Brain_Segmentation_Volume_in_mm3.csv');
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

figure;imagesc(outnet_edges)   
  
Session =[]; Dat =[]; Age =[]; permAge=[];

Pred_Age_nullFC1=[]; Pred_Age_nullFC2=[]; Pred_Age_nullFC3=[]; Pred_Age_nullFC4=[];
Pred_Age_nullFC5=[]; Pred_Age_nullFC6=[]; Pred_Age_nullFC7=[]; Pred_Age_nullFC8=[];
Pred_Age_nullFC9=[]; Pred_Age_nullFC10=[]; Pred_Age_nullFC11=[];
Pred_Age_betFC1=[]; Pred_Age_witFC1=[]; Pred_Age_betFC2=[]; Pred_Age_witFC2=[];
Pred_Age_betFC3=[]; Pred_Age_witFC3=[]; Pred_Age_betFC4=[]; Pred_Age_witFC4=[];
Pred_Age_betFC5=[]; Pred_Age_witFC5=[]; Pred_Age_betFC6=[]; Pred_Age_witFC6=[];
Pred_Age_betFC7=[]; Pred_Age_witFC7=[]; Pred_Age_betFC8=[]; Pred_Age_witFC8=[];
Pred_Age_betFC9=[]; Pred_Age_witFC9=[]; Pred_Age_betFC10=[]; Pred_Age_witFC10=[];
Pred_Age_betFC11=[]; Pred_Age_witFC11=[]; 
FDs =[]; aQC =[]; fQC =[]; Volumes=[];TRs=[];
 for zx = 1:500   
    all_ages=[]; all_subid=[]; all_Vol=[]; all_fqc =[];
    all_aqc=[]; all_fd =[]; all_TR=[];
    
    all_braincondat_withnet1=[]; all_braincondat_withnet2=[]; all_braincondat_withnet3=[];
all_braincondat_withnet4=[]; all_braincondat_withnet5=[]; all_braincondat_withnet6=[];
all_braincondat_withnet7=[]; all_braincondat_withnet8=[]; all_braincondat_withnet9=[];
all_braincondat_withnet10=[]; all_braincondat_withnet11=[]; 

    all_braincondat_betwnet1=[]; all_braincondat_betwnet2=[]; all_braincondat_betwnet3=[];
all_braincondat_betwnet4=[]; all_braincondat_betwnet5=[]; all_braincondat_betwnet6=[];
all_braincondat_betwnet7=[]; all_braincondat_betwnet8=[]; all_braincondat_betwnet9=[];
all_braincondat_betwnet10=[]; all_braincondat_betwnet11=[]; 
    % ##########
%     rng 'default'
    rng 'shuffle'
    % ##########
    randouts = cell(11,1);
    for k=1:11
        
        randouts{k} = randperm(49261,length(innet_edges{k}));
        
    end

    
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
%            disp(sub_data.name); 
%            a= size(sub_data.GordonMatAP,2)
%            b = size(sub_data.GordonMatPA,2)
       tseries = [sub_data.GordonMat(1:333,1:600)];  
       end
       totalcondat1 =  corr(tseries(All_parcels,1:600)');
       totalconndat = totalcondat1(All_edges);
       outconndat = totalcondat1(find(outnet_edges==1)); 
       for k=1:11
           temp_condat = corr(tseries(innet_parcels{k},1:600)');
           tmp = temp_condat(innet_edges{k})';
           eval(['inconndat',num2str(k),' =  tmp;']);      
           
%            tmp = outconndat(randperm(49261,length(innet_edges{k})))';
           tmp = outconndat(randouts{k})';
           eval(['outconndat',num2str(k),' =  tmp;']);
       end
       for k=1:11
       eval(['all_braincondat_withnet',num2str(k),' =[all_braincondat_withnet',num2str(k),'; inconndat',num2str(k),'];']);
       eval(['all_braincondat_betwnet',num2str(k),' =[all_braincondat_betwnet',num2str(k),'; outconndat',num2str(k),'];']);
       end

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

        Y = all_ages; n=size(Y,1); all_permages = Y(randperm(length(Y)));
        
        for k=1:11
        X = eval(['[atanh(all_braincondat_withnet',num2str(k),')];' ]);
        MdlLin = fitrsvm(X,all_permages,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        eval(['YHat',num2str(k),' = kfoldPredict(MdlLin);']); 
        end
        
        for k=1:11
        eval(['Pred_Age_nullFC',num2str(k),' = [Pred_Age_nullFC',num2str(k),'; YHat',num2str(k),'];']); 
        end
        
        for k=1:11
        X = eval(['[atanh(all_braincondat_withnet',num2str(k),')];' ]);
        MdlLin = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        eval(['YHat',num2str(k),' = kfoldPredict(MdlLin);']); 
        end
        
        for k=1:11
        eval(['Pred_Age_witFC',num2str(k),' = [Pred_Age_witFC',num2str(k),'; YHat',num2str(k),'];']); 
        end
        
        for k=1:11
        X = eval(['[atanh(all_braincondat_betwnet',num2str(k),')];' ]);
        MdlLin = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
        eval(['YHat',num2str(k),' = kfoldPredict(MdlLin);']); 
        end

        for k=1:11
        eval(['Pred_Age_betFC',num2str(k),' = [Pred_Age_betFC',num2str(k),'; YHat',num2str(k),'];']); 
        end
        
        Session =[Session; all_subid]; 
        Dat =[Dat; zx.*ones(n,1)];
        Age =[Age; all_ages]; 
        permAge=[permAge; all_permages];
        Volumes =[Volumes; all_Vol];
        FDs =[FDs; all_fd]; 
        aQC =[aQC; all_aqc]; 
        fQC =[fQC; all_fqc];
        TRs =[TRs; all_TR];
    
    zx
 end
T = table(Session,Dat,Age,permAge,Volumes,...
    Pred_Age_nullFC1,Pred_Age_nullFC2,Pred_Age_nullFC3,Pred_Age_nullFC4,Pred_Age_nullFC5,...
    Pred_Age_nullFC6,Pred_Age_nullFC7,Pred_Age_nullFC8,Pred_Age_nullFC9,Pred_Age_nullFC10,Pred_Age_nullFC11,...
    Pred_Age_witFC1,Pred_Age_witFC2,Pred_Age_witFC3,Pred_Age_witFC4,Pred_Age_witFC5,...
    Pred_Age_witFC6,Pred_Age_witFC7,Pred_Age_witFC8,Pred_Age_witFC9,Pred_Age_witFC10,Pred_Age_witFC11,...
    Pred_Age_betFC1,Pred_Age_betFC2,Pred_Age_betFC3,Pred_Age_betFC4,Pred_Age_betFC5,...
    Pred_Age_betFC6,Pred_Age_betFC7,Pred_Age_betFC8,Pred_Age_betFC9,Pred_Age_betFC10,Pred_Age_betFC11,...
    FDs, aQC, fQC, TRs);
writetable(T,'singleNet_AgePredStacked_1_500.csv');

%% visualization of the distributions
% clear all
T = readtable('singleNet_AgePredStacked_1_500.csv');
dist_stats_adj =[]; dist_stats_nadj =[]; 
 load('networknames.mat');
 netsizes = [50 48 26 26 18 17 28 54 28 33 5];
 [kl ord] = sort(netsizes,'Descend');
xnames =[networknames(ord(1:10));"None"]; 
for k=1:length(unique(T.Dat))

datsamp = T(T.Dat==k,:);
datsamp.volsq = (datsamp.Volumes-nanmean(datsamp.Volumes)).^2;
Y = datsamp.Age;

% random equal size from between net
r0s = [];
r00s = [];
for jj=1:11
    j = ord(jj);
eval(['YHat = datsamp.Pred_Age_betFC',num2str(j),';']);
n = length(~isnan(YHat));
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq1 = 1 -(mseLin./errvar);
[r0 p ] = corr(Y,YHat,'Rows','complete');
[r00 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');
r0s = [r0s r0];
r00s = [r00s r00];
end

% permuted null
r1s = [];
r2s = [];
for jj=1:11
    j = ord(jj);
eval(['YHat = datsamp.Pred_Age_nullFC',num2str(j),';']);
n = length(~isnan(YHat));
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq1 = 1 -(mseLin./errvar);
[r1 p ] = corr(Y,YHat,'Rows','complete');
[r2 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');
r1s = [r1s r1];
r2s = [r2s r2];
end

% single within net

r3s = [];
r4s = [];
for jj=1:11
    j = ord(jj);
eval(['YHat = datsamp.Pred_Age_witFC',num2str(j),';']);
n = length(~isnan(YHat));
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq1 = 1 -(mseLin./errvar);
[r3 p ] = corr(Y,YHat,'Rows','complete');
[r4 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');
r3s = [r3s r3];
r4s = [r4s r4];
end

dist_stats_nadj = [dist_stats_nadj; [r1s r3s r0s]];
dist_stats_adj = [dist_stats_adj; [r2s r4s r00s]];
k
end
%% Figure 6 single networks and null models predicting age
for k=1:11
pvals(k) = 1 - length(find(dist_stats_adj(:,k) < dist_stats_adj(:,k+11)))./length(unique(T.Dat));
xnames(k)
end
figure;
distributionPlot({dist_stats_adj(:,12)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','right','showMM',0)
distributionPlot({dist_stats_adj(:,1)},...
    'color',{[0.7 0.7 0.7]},'xyOri','normal','histOri','right','showMM',0,'xNames',xnames)

legend({'Single-network FC Model','Single-network Null Model'},'AutoUpdate','off','FontSize',16);
distributionPlot({dist_stats_adj(:,12);dist_stats_adj(:,13);dist_stats_adj(:,14);...
    dist_stats_adj(:,15);dist_stats_adj(:,16);dist_stats_adj(:,17);dist_stats_adj(:,18);...
    dist_stats_adj(:,19);dist_stats_adj(:,20);dist_stats_adj(:,21);dist_stats_adj(:,22)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','right','showMM',0,'xNames',xnames)
hold on

distributionPlot({dist_stats_adj(:,1);dist_stats_adj(:,2);dist_stats_adj(:,3);...
    dist_stats_adj(:,4);dist_stats_adj(:,5);dist_stats_adj(:,6);dist_stats_adj(:,7);...
    dist_stats_adj(:,8);dist_stats_adj(:,9);dist_stats_adj(:,10);dist_stats_adj(:,11)},...
    'color',{[0.7 0.7 0.7]},'xyOri','normal','histOri','right','showMM',0)

line([0 ,12],[.76 ,.76],'Color','blue'); % all within-network edges median partial r
ylim([-.4,1])
xlabel('Network in the FC model');ylabel('Correlation between predicted and true age (partial r)');
set(gca,'FontSize',16,'XTickLabelRotation',45)
plot(0:12,zeros(1,13),'color','k');

meds = median(dist_stats_adj(:,12:22));
for k =1:11
    if pvals(k)<=.01
text(k-.1,meds(k),'*','FontSize',25);
    end
end

%% Figure S6 single networks and control models predicting age
for k=1:11
pvals(k) = 1 - length(find(dist_stats_adj(:,k) < dist_stats_adj(:,k+11)))./length(unique(T.Dat));
xnames(k)
end
figure;
distributionPlot({dist_stats_adj(:,23)},...
    'color',{[0.60,0.75,0.73]},'xyOri','normal','histOri','right','showMM',0)
distributionPlot({dist_stats_adj(:,12)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','right','showMM',0,'xNames',xnames)

legend({'Size-matched Control Model','Single-network FC Model'},'AutoUpdate','off','FontSize',16);
distributionPlot({dist_stats_adj(:,23);dist_stats_adj(:,24);dist_stats_adj(:,25);...
    dist_stats_adj(:,26);dist_stats_adj(:,27);dist_stats_adj(:,28);dist_stats_adj(:,29);...
    dist_stats_adj(:,30);dist_stats_adj(:,31);dist_stats_adj(:,32);dist_stats_adj(:,33)},...
    'color',{[0.60,0.75,0.73]},'xyOri','normal','histOri','right','showMM',0,'xNames',xnames)
hold on

distributionPlot({dist_stats_adj(:,12);dist_stats_adj(:,13);dist_stats_adj(:,14);...
    dist_stats_adj(:,15);dist_stats_adj(:,16);dist_stats_adj(:,17);dist_stats_adj(:,18);...
    dist_stats_adj(:,19);dist_stats_adj(:,20);dist_stats_adj(:,21);dist_stats_adj(:,22)},...
    'color',{[0.70,0.75,0.93]},'xyOri','normal','histOri','right','showMM',0,'xNames',xnames)

xlabel('Network in the FC model');ylabel('Correlation between predicted and true age (partial r)');
set(gca,'FontSize',16,'XTickLabelRotation',45)
plot(0:12,zeros(1,13),'color','k');

line([0 ,12],[.76 ,.76],'Color','blue'); % all within-network edges median partial r
ylim([-.4,1])

meds = median(dist_stats_adj(:,12:22));
for k =1:11
    if pvals(k)<=.01
text(k-.1,meds(k),'*','FontSize',25);
    end
end