% Kardan, O., Kaplan, S., ..., & Rosenberg, M.D. (2022) 
% "Resting-state functional connectivity identifies individuals and
% predicts age in 8-to-26-month-olds" Developmental Cognitive Neuroscience

% Intraclass correlation of edges (ICC), Differential power, and scripts
% for generating supplemtary Figure S1 and supplemtary Figure S4

%  Thomas Zoeller (2022). Intraclass correlation coefficient with confidence intervals 
% (https://www.mathworks.com/matlabcentral/fileexchange/26885-intraclass-correlation-coefficient-with-confidence-intervals),
% MATLAB Central File Exchange. Retrieved June 15, 2022.

clear all
addpath(genpath('icc21')); % the icc21 matlab package by Thomas Zoeller referenced above
load('goodGMDec2.mat')
   load('mask.mat'); numnets = 11;
    goodareas1 = [1:333]; subcort=0;% only corticals
    All_parcels = goodareas1(mask);   %ordered 
    All_edges = find(tril(ones(333),-1)==1);
  

    all_braincondat_fullnetAP=[]; all_braincondat_fullnetPA=[]; all_ages=[]; all_subid=[]; 

    subs =[];
    for i =1:length(goodGMDec)
        subs = [subs; string(goodGMDec(i).name(5:10))];
    end
    
    uniqs = unique(subs);
    for j=1:length(uniqs)
        idx = find(subs==uniqs(j));
        ii=randi(length(idx));
       sub_data = goodGMDec(idx(ii));
      

       APconns = corr([sub_data.GordonMatAP(All_parcels,1:min([300,size(sub_data.GordonMatAP,2)]) )]');
       PAconns = corr([sub_data.GordonMatPA(All_parcels,1:min([300,size(sub_data.GordonMatPA,2)]) )]');
 
 
       totalconndatAP = APconns(All_edges);
       totalconndatPA = PAconns(All_edges);
       
        all_braincondat_fullnetAP = [all_braincondat_fullnetAP; atanh(totalconndatAP)'];
        all_braincondat_fullnetPA = [all_braincondat_fullnetPA; atanh(totalconndatPA)'];
        
        all_ages = [all_ages; sub_data.Scanage];
        
      j
    end
    %% edge-wise ICC
ICCs = [];
    for k=1:55278
        edat = [all_braincondat_fullnetAP(:,k)  all_braincondat_fullnetPA(:,k)];
        R = icc21(edat);
        ICCs = [ICCs; R];
        
        if k/1000 == floor(k/1000)
            disp(num2str(k));
        end
    end
    ICCs(ICCs<0)=0;
 mean(ICCs)
 std(ICCs)
    
    %% correlation of ICC of edges with their age-prediction SVR beta
    age_ICC_r = [];
    for mm=1:50
        load(['Example_trained_Models\mdl_',num2str(mm),'.mat']);
        C1 = MdlLin1;
        for j=1:10
            age_ICC_r = [age_ICC_r ; corr(ICCs,C1.Trained{j, 1}.Beta)];
            
        end
    end
    mean(age_ICC_r)
    std(age_ICC_r)
    
    %% Differential Power (Section 5 of the supplemetary results)
    DP = NaN(55278,112);
    
    all_braincondat_fullnetAP = zscore(all_braincondat_fullnetAP,0,2);
    all_braincondat_fullnetPA = zscore(all_braincondat_fullnetPA,0,2);
    
    betweenR1 = zeros(112,112,55278);
    for j=1:112
        for k=1:112
    betweenR1(j,k,:) = all_braincondat_fullnetAP(j,:).*all_braincondat_fullnetAP(k,:);
        end
    end
    
    betweenR2 = zeros(112,112,55278);
    for j=1:112
        for k=1:112
    betweenR2(j,k,:) = all_braincondat_fullnetPA(j,:).*all_braincondat_fullnetPA(k,:);
        end
    end
    
    withinR = all_braincondat_fullnetAP.*all_braincondat_fullnetPA;
    
    for e =1:55278
        for i =1:112
            count = length(find(squeeze(betweenR1(i,setdiff(1:112,i),e)) > withinR(i,e) |...
                squeeze(betweenR2(i,setdiff(1:112,i),e)) > withinR(i,e) )); 
  DP(e,i) = -log(count/(2*111));
        end
    end
 %% Differential Power (Section 5 of the supplemetary results)
 DP(isinf(DP)) = NaN;
 DPe = nanmean(DP')';
  
     age_DP_r = [];
    for mm=1:50
        load(['Example_trained_Models\mdl_',num2str(mm),'.mat']);
        C1 = MdlLin1;
        for j=1:10
            age_DP_r = [age_DP_r ; corr(DPe,C1.Trained{j, 1}.Beta)];
            
        end
    end
    mean(age_DP_r)
    std(age_DP_r)   
    %% Age-predictive edges and DP
    age_svr_beta = zeros(55278,1);
    for mm=1:50
        load(['Example_trained_Models\mdl_',num2str(mm),'.mat']);
        C1 = MdlLin1;
        for j=1:10
            age_svr_beta = age_svr_beta + C1.Trained{j, 1}.Beta;
            
        end
    end
    age_svr_beta = age_svr_beta/500;
    
    %% make a figure with DPs and age Betas (Figure S4)
        load('network.mat')
    load('networknames.mat')
    GG1 = zeros(333);
GG1(All_edges) = zscore(age_svr_beta);

    GG2 = zeros(333);
GG2(All_edges) = zscore(DPe);
%   GG3(All_edges) = mean(adultnet_FCs);

GG3 = GG2 + GG1';
%GG3(333,333)=-1;GG3(1,1)=1; % for normalized colorbar
figure;
imagesc(GG3);colormap(bluewhitered), colorbar;
xnames =[networknames(1:10);string('None')]; 
netsizes1 = [0 50 48 26 26 18 17 28 54 28 33];
netsizes = [50 48 26 26 18 17 28 54 28 33 5];
%%%
% xnames = GAnetnames;
% netsizes1 = [0 24    40     5    41    32    24    47     8    38     8     4    23];
% netsizes = [24    40     5    41    32    24    47     8    38     8     4    23    39];
%%%

set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',13);
set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',13);

makans = cumsum(netsizes);
for j=1:length(netsizes)
line([makans(j),makans(j)],[0 ,333],'Color','black');
end
for j=1:length(netsizes)
line([0 ,333],[makans(j),makans(j)],'Color','black');
end
title('DP (below diagonal) and age Beta (above diagonal)');
% title('Average of rsFC (Gordon 2017 network assignments)');
line([0 ,333],[0 ,333],'Color','black');
ylabel('z-scored DP\Beta value','Position',[381,167,1]); axis square

    %% make a figure with AP and PA FCs of random babies (Figure S1)
    figure;
    randsubs = [1,2,12,19, 30,61]; % picked psudo-random subs but from different ages
    for k =1:6
        subplot(2,3,k)
        GG1 = zeros(333);
        GG1(All_edges) = all_braincondat_fullnetAP(randsubs(k),:) ;
        
        GG2 = zeros(333);
        GG2(All_edges) = all_braincondat_fullnetPA(randsubs(k),:);
        %   GG3(All_edges) = mean(adultnet_FCs);
        rsph = corr(all_braincondat_fullnetAP(randsubs(k),:)',all_braincondat_fullnetPA(randsubs(k),:)');
        GG3 = GG2 + GG1';
        GG3(333,333)=-4;GG3(1,1)=7; % for normalized colorbar
        
        imagesc(GG3);colormap(bluewhitered), colorbar;
        xnames =[networknames(1:10);string('None')];
        netsizes1 = [0 50 48 26 26 18 17 28 54 28 33];
        netsizes = [50 48 26 26 18 17 28 54 28 33 5];
        %%%
        % xnames = GAnetnames;
        % netsizes1 = [0 24    40     5    41    32    24    47     8    38     8     4    23];
        % netsizes = [24    40     5    41    32    24    47     8    38     8     4    23    39];
        %%%
        
        set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
            'XtickLabelRotation',45,'FontSize',13);
        set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
            'XtickLabelRotation',45,'FontSize',13);
        
        makans = cumsum(netsizes);
        for j=1:length(netsizes)
            line([makans(j),makans(j)],[0 ,333],'Color','black');
        end
        for j=1:length(netsizes)
            line([0 ,333],[makans(j),makans(j)],'Color','black');
        end
        title(['Subject ',num2str(k),' Age = ',num2str(all_ages(randsubs(k))),...
            ' mo', '  r_{split-half} = ',num2str(round(rsph*100)/100)]);
        % title('Average of rsFC (Gordon 2017 network assignments)');
        line([0 ,333],[0 ,333],'Color','black');
        ylabel('atanh(r)','Position',[381,167,1]); axis square
    end
    %%
  %export_fig -m 4 -transparent  