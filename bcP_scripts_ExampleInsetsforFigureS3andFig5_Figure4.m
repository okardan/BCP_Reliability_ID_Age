% Kardan, O., Kaplan, S., ..., & Rosenberg, M.D. (2022) 
% "Resting-state functional connectivity identifies individuals and
% predicts age in 8-to-26-month-olds" Developmental Cognitive Neuroscience

% 1 example for visualization used in Figure 5 and Figure S3 Scatterplot insets 
T = readtable('AgePredStacked.csv');

figure
   xylabel = 'Age (months)';
   
datsamp = T(T.Dat==1,:); Y = datsamp.Age;

% subplot(1,3,1)
c1 = [0.5,0.5,.91];c2 =  [0.70,0.75,0.93];
YHat = datsamp.Pred_Age_totFC;
n = length(~isnan(YHat));
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq = 1 -(mseLin./errvar);
scatter(Y,YHat,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
xlabel(['True ',xylabel],'FontSize',15);ylabel(['Predicted ',xylabel],'FontSize',15);
xlim([.9*min(Y),1.1*max(Y)]);ylim([.9*min(Y),1.1*max(Y)]); set(gca,'FontSize',15);axis square;
title('total FC');

[r p ] = corr(Y,YHat,'Rows','complete');
[r2 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');

text(1.2*min(Y),1.7*mean(YHat),...
    {['prediction R^2 = ',num2str(.001*round(1000*rsq)),'; r = ',num2str(.001*round(1000*r))]...
    ,['RMSE = ',num2str(.01*round(100*sqrt(mseLin))),'; n = ',num2str(n)],...
    ['partial r = ',num2str(.001*round(1000*(r2))),'; p = ',num2str(p2)]},'FontSize',14);

% between net
figure 
YHat = datsamp.Pred_Age_betFC;
c1 = [0,1,0];c2 = [.3 .9 .7];
mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq = 1 -(mseLin./errvar);
scatter(Y,YHat,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
xlabel(['True ',xylabel],'FontSize',15);ylabel(['Predicted ',xylabel],'FontSize',15);
xlim([.9*min(Y),1.1*max(Y)]);ylim([.9*min(Y),1.1*max(Y)]); set(gca,'FontSize',15);axis square;
title('between-net FC');

[r p ] = corr(Y,YHat,'Rows','complete'); 
[r2 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC datsamp.TRs],'Rows','complete');
text(1.2*min(Y),1.7*mean(YHat),...
    {['prediction R^2 = ',num2str(.001*round(1000*rsq)),'; r = ',num2str(.001*round(1000*r))]...
    ,['RMSE = ',num2str(.01*round(100*sqrt(mseLin))),'; n = ',num2str(n)],...
    ['partial r = ',num2str(.001*round(1000*(r2))),'; p = ',num2str(p2)]},'FontSize',14);

% within net
figure  
 YHat = datsamp.Pred_Age_witFC;
    c1 = [1,0,0];c2 = [.9 .7 .3];
    mseLin = sum((Y-YHat).^2)/n;
errvar = sum(((Y-mean(Y)).^2))/n;
rsq = 1 -(mseLin./errvar);
     scatter(Y,YHat,'MarkerEdgeColor',c2,'MarkerFaceColor',c2);
    xlabel(['True ',xylabel],'FontSize',15);ylabel(['Predicted ',xylabel],'FontSize',15);
    xlim([.9*min(Y),1.1*max(Y)]);ylim([.9*min(Y),1.1*max(Y)]); set(gca,'FontSize',15);axis square;
 title('within-net FC');  
 
    [r p ] = corr(Y,YHat,'Rows','complete');
  [r2 p2 ]= partialcorr(Y,YHat,[datsamp.FDs datsamp.aQC datsamp.fQC  datsamp.TRs],'Rows','complete');  
text(1.2*min(Y),1.7*mean(YHat),...
    {['prediction R^2 = ',num2str(.001*round(1000*rsq)),'; r = ',num2str(.001*round(1000*r))]...
    ,['RMSE = ',num2str(.01*round(100*sqrt(mseLin))),'; n = ',num2str(n)],...
    ['partial r = ',num2str(.001*round(1000*(r2))),'; p = ',num2str(p2)]},'FontSize',14);
datsamp.volsq = (datsamp.Volumes-nanmean(datsamp.Volumes)).^2;

%% network plots avergae FC over all sample in Figure 4

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
    
    
        load('GAmask.mat'); 
    load('GAnetnames.mat');
    load('GAnetwork.mat')
        goodareas1 = [1:333]; subcort=0;% only corticals
   All_parcels_adult = goodareas1(GAmask);   %ordered 
    All_edges_adult = find(tril(ones(333),-1)==1);   
    numnets_adults = 13;


    innet_parcels_adults = cell(1,numnets_adults);
    innet_edges_adults = cell(1,numnets_adults);
    for ntw=1:numnets_adults
    innet_parcels_adults{ntw} =  GAmask(GAnetwork ==ntw);
    featsz = length(find(GAnetwork ==ntw));
    innet_edges_adults{ntw} = find(tril(ones(featsz),-1));
    end
    outnet_parcels_adults = All_parcels_adult;
    outnet_edges_adults = zeros(333);
    for ntw1=1:numnets_adults-1
        for ntw2=ntw1+1:numnets_adults
            
            for ii=1:333
                for jj=1:333
                    if GAnetwork(ii)==ntw1 & GAnetwork(jj)==ntw2
                        outnet_edges_adults(jj,ii) = 1;
                    end
                end
            end
        end
    end
%% average rsFC plots for each networks assignments in Figure 4
load('goodGMDec2.mat')
% The colormap used is from:
% Nathan Childress (2022). bluewhitered (https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered)
% download colormap above and add to path before running this section
    addpath(genpath('bluewhitered')) 
    babynet_FCs=[]; adultnet_FCs =[];
    for k=1:170
        ses_data = goodGMDec(k);
        try
            tseries = [ses_data.GordonMatAP(1:333,1:300) ses_data.GordonMatPA(1:333,1:300)];
        catch
            tseries = [ses_data.GordonMat(1:333,1:600)];
        end
        totalcondat1 =  corr(tseries(All_parcels,1:600)');
        totalconndat = totalcondat1(All_edges);
        babynet_FCs = [babynet_FCs; totalconndat'];
        totalcondat1 =  corr(tseries(All_parcels_adult,1:600)');
        totalconndat = totalcondat1(All_edges_adult);
        adultnet_FCs = [adultnet_FCs; totalconndat'];
    end
  GG3 = zeros(333);
GG3(All_edges) = mean(babynet_FCs);
% GG3(All_edges) = mean(adultnet_FCs); % uncomment to use Gordon 2016 adult networks instead

  
GG3(333,333)=-1;GG3(1,1)=1; % for normalized colorbar
figure;
imagesc(GG3);colormap(bluewhitered), colorbar;
xnames =[networknames(1:10);string('None')]; 
netsizes1 = [0 50 48 26 26 18 17 28 54 28 33];
netsizes = [50 48 26 26 18 17 28 54 28 33 5];
%%%
% uncomment the next three lines to use Gordon 2016 adult networks instead
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
line([makans(j),makans(j)],[makans(j) ,333],'Color','black');
end
for j=1:length(netsizes)
line([0 ,makans(j)],[makans(j),makans(j)],'Color','black');
end
title('Average of rsFC (Baby network assignments)');
% uncomment the next line to use Gordon 2016 adult networks instead
% title('Average of rsFC (Gordon 2016 network assignments)');
line([0 ,333],[0 ,333],'Color','black');
ylabel('average correlation value','Position',[381,167,1]); axis square

wids = find(tril(ones(333),-1)==1 & outnet_edges ==0);
bids = find(outnet_edges ==1);
% uncomment the next two lines to use Gordon 2016 adult networks instead
% wids = find(tril(ones(333),-1)==1 & outnet_edges_adults ==0);
% bids = find(outnet_edges_adults ==1);
mean(GG3(wids))
std(GG3(wids))
mean(GG3(bids))
std(GG3(bids))
figure;subplot(1,2,1);imagesc(outnet_edges);subplot(1,2,2);imagesc(outnet_edges_adults);
%% within and between net masked feature images for the example sample in Figure S3 inset
load('example_trainedModel_withinFC.mat')
load('example_trainedModel_betwFC.mat')
Xwith = MdlLin2.X; Y = MdlLin2.Y;
Xbetw = MdlLin3.X; Y = MdlLin3.Y;

wids = find(tril(ones(333),-1)==1 & outnet_edges ==0);
bids = find(outnet_edges ==1);

Bmapsw = zeros(size(MdlLin2.Trained{1, 1}.Beta));
for f=1:10
betamap = MdlLin2.Trained{f, 1}.Beta;
Bmapsw = Bmapsw + betamap/10;
end
Bmapsb = zeros(size(MdlLin3.Trained{1, 1}.Beta));
for f=1:10
betamap = MdlLin3.Trained{f, 1}.Beta;
Bmapsb = Bmapsb + betamap/10;
end

GGw = zeros(333);GGw(wids) = zscore(Bmapsw); GGw(1,1)=-4;GGw(333,333)=4.5;

GGb = zeros(333);GGb(bids) = zscore(Bmapsb); GGb(1,1)=-4;GGb(333,333)=4.5;
xnames =[networknames(1:10);string('None')]; 
  netsizes1 = [0 50 48 26 26 18 17 28 54 28 33];
 netsizes = [50 48 26 26 18 17 28 54 28 33 5];
 
figure;
imagesc(GGw);colormap(bluewhitered), colorbar;
% imagesc(GGb);colormap(bluewhitered), colorbar;


set(gca,'XTick',round(cumsum(netsizes1)+ netsizes/2),'Xticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',13);
set(gca,'YTick',round(cumsum(netsizes1)+ netsizes/2),'Yticklabel',xnames,...
         'XtickLabelRotation',45,'FontSize',13);

makans = cumsum(netsizes);
for j=1:numnets
line([makans(j),makans(j)],[makans(j),333],'Color','black');
end
for j=1:numnets
line([0 ,makans(j)],[makans(j),makans(j)],'Color','black');
end
title('Beta Map from 10 Folds of SVR with only Between-Net Edges');
line([0 ,333],[0 ,333],'Color','black');
ylabel('z-scored SVR normalized beta','Position',[381,167,1]); axis square

