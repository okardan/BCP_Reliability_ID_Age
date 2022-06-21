% Kardan, O., Kaplan, S., ..., & Rosenberg, M.D. (2022) 
% "Resting-state functional connectivity identifies individuals and
% predicts age in 8-to-26-month-olds" Developmental Cognitive Neuroscience

% This script runs the SVR on permuted age using full connectomes or
% only within-net or only between-net edges with 10-fold cross validation
% to make null distributions nullAgePredStacked.csv which is used in
% bcp_scripts_AgePredModels_Figure2_Figure3_Figure5.m

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

   
%%   
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
%            disp(sub_data.name); 
%            a= size(sub_data.GordonMatAP,2)
%            b = size(sub_data.GordonMatPA,2)
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
%         mseLin = kfoldLoss(MdlLin);
        YHat2 = kfoldPredict(MdlLin2); n=size(Y,1);
%         errvar = sum(((Y-mean(Y)).^2))/(n-1);
%         rsq = 1 -(mseLin./errvar)

        X = [atanh(all_braincondat_betwnet)]; 
        MdlLin3 = fitrsvm(X,Y,'Standardize',true,'BoxConstraint',Inf,'Epsilon',0.00001,'KFold',10);
%         mseLin = kfoldLoss(MdlLin);
        YHat3 = kfoldPredict(MdlLin3); n=size(Y,1);
%         errvar = sum(((Y-mean(Y)).^2))/(n-1);
%         rsq = 1 -(mseLin./errvar)

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