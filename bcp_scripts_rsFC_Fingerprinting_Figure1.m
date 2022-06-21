% Kardan, O., Kaplan, S., ..., & Rosenberg, M.D. (2022) 
% "Resting-state functional connectivity identifies individuals and
% predicts age in 8-to-26-month-olds" Developmental Cognitive Neuroscience

% Fingerprinting analyses and  stats for making Figure 1
clear all
load('goodGMDec2.mat')  %goodGMDec2.mat contains the preprocessed fMRI parcellated timeseries for the 170 fMRI sessions
%%
sescount = 170;

inds = find(tril(ones(333),-1));
incparcs = 1:333;
subids=[];
omrs =[]; ters =[];
for j=1:sescount
    subid = str2num(goodGMDec(j).name(5:10));
    subids = [subids;subid];
    omrs = [omrs;goodGMDec(j).Scanage];
    ters = [ters;goodGMDec(j).TR];
end
accs =[];ns=[];
for frm =[5:5:25,50:50:400];
    
    Corrs2 = NaN(sescount);Ages =[];
    for j=1:sescount
        l1 = size(goodGMDec(j).GordonMatAP,2);
        l2 = size(goodGMDec(j).GordonMatPA,2);
        if l1<frm | l2<frm
            continue
        end
        firsthalf = [goodGMDec(j).GordonMatAP(incparcs,1:frm)];
        secondhalf = [goodGMDec(j).GordonMatPA(incparcs,1:frm)];
        apcor = corr(firsthalf');
        apcorr = apcor(inds);
        pacor = corr(secondhalf');
        pacorr = pacor(inds);
        dir_corr = corr(atanh(apcorr),atanh(pacorr));
        
        if j>1
            for k=1:j-1
                
                l1 = size(goodGMDec(k).GordonMatAP,2);
                l2 = size(goodGMDec(k).GordonMatPA,2);
                if l1<frm | l2<frm
                    continue
                end
                
                firsthalf = [goodGMDec(k).GordonMatAP(incparcs,1:1:frm)];
                secondhalf = [goodGMDec(k).GordonMatPA(incparcs,1:frm)];
                apcor_null = corr(firsthalf');
                apcorr_null = apcor_null(inds);
                pacor_null = corr(secondhalf');
                pacorr_null = pacor_null(inds);
                dir_corr_null1 = corr(atanh(apcorr),atanh(pacorr_null));
                dir_corr_null2 = corr(atanh(pacorr),atanh(apcorr_null));
                
                Corrs2(k,j) = dir_corr_null1; Corrs2(j,k) = dir_corr_null2;
            end
        end
        
        Corrs2(j,j) = dir_corr;
        Ages = [Ages; goodGMDec(j).Scanage];
        j
    end
    
    sames = Corrs2(find(eye(sescount)));
    
    max_null =NaN(sescount,1);
    for i=1:sescount
        elses = find( subids ~= subids(i));
        max_null(i) = max([max(Corrs2(elses,i)); max(Corrs2(i,elses))]);
    end
    
    accs = [accs; length(find([sames-max_null]>0))./length(find(~isnan(sames)))];
    ns = [ns; length(find(~isnan(sames)))];
    if frm == 300
    save('Corrs2.mat','Corrs2');
    save('Ages.mat','Ages');
    end
end
fprnt_accs_ns = [accs ns];
save('fprnt_accs_ns2.mat','fprnt_accs_ns');

%% Figure 1 Panel (B) correlation between sessions’ split-half functional connectivity patterns at 300 frames
load('fprnt_accs_ns2.mat') % accs across different inc frames
load('Corrs2.mat') % full mat for 300-frame splits
load('fingerprinting\Ages.mat');
figure
Corrs2_age = Corrs2(~isnan(Corrs2(:,1)),~isnan(Corrs2(:,1)));
[a,b] = sort(Ages);
imagesc(Corrs2_age(b,b));ag = Ages(b);
set(gca,'XTick',[],'XTickLabel',string([]),...
    'YTick',[],'YTickLabel',string([]),...
    'XTickLabelRotation',45,'FontSize',13)
%     xlabel('Sessions'); ylabel('Sessions');
agsizes =[];
for k=[8:26]
    agsizes = [agsizes;length(find(ag==k))];
end
for t=[1:10,12,15,17,19]
    tt = cumsum(agsizes);
    text(tt(t),159,[num2str(ag(tt(t))),' mo'],'rotation',45,'FontSize',11,'HorizontalAlignment','right');
    text(-1,tt(t),[num2str(ag(tt(t))),' mo'],'rotation',45,'FontSize',11,'HorizontalAlignment','right');
end
colormap(jet), colorbar; axis square

%% Figure 1 panel (C)  AP-PA halves correlations for the same vs other participants
Corrs3 = Corrs2_age(b,b);
ids = ~isnan(Corrs2(:,1));
subidss = subids(ids);
figure
scatter(1:158,Corrs3(find(eye(158))),'r','LineWidth',1.8)
ylim([.1,.8]);
hold on
mms =[];
for i=1:158
    
    elses = find( subidss(b) ~= subidss(b(i)));
    sameelses = find( subidss(b) == subidss(b(i)) & ag(b)~=ag(b(i)));
    mins = min([min(Corrs3(i,elses)); min(Corrs3(elses,i))]);
    maxs = max([max(Corrs3(elses,i)); max(Corrs3(i,elses))]);
    
    line([i,i],[mins,maxs],'LineWidth',1.2 )
    
    if length(sameelses)>0
        scatter(i*ones(1,length(sameelses)),Corrs3(i,sameelses),160,'.k');
        scatter(i*ones(1,length(sameelses)),Corrs3(sameelses,i),160,'.k');
    end
    if i==1
        legend({'same participant & age','same participant other age','other participants'...
            },'AutoUpdate','off','FontSize',16);
    end
    mms=[mms;maxs];
end
xlim([0,159]);
% xlabel('Session');
ylabel('Split-half FC correlation (r)');ylim([0,1])
set(gca,'FontSize',16,'XTick',[],'XTickLabel',string([]),...
    'XTickLabelRotation',45);
for t=[1:19]
    tt = cumsum(agsizes);
    text(tt(t),-.01,[num2str(ag(tt(t))),' mo'],'rotation',45,'FontSize',13,'HorizontalAlignment','right');
    
end
hold off

%% Figure 1 panel (A) fingerprinting ID accuracy as a function of number of frames
figure
plot(1:13,fprnt_accs_ns(:,1)'.*100); hold on
l=[];
for i=1:sescount
    l= [l; size(goodGMDec(i).GordonMat,2)];
end
nsess = [sescount*ones(1,11) length(find(l>700))    length(find(l>800)) ];
ns = [fprnt_accs_ns(:,2)' ];
nss = 100*ns./(sescount*2);

plot( 1:13, nss)
for j=[6,10,11,12,13]
%     text(j,nss(j),['n = ',num2str(ns(j)),'/',num2str(nsess(j))])
     text(j,nss(j),['n = ',num2str(ns(j))])
end
line(1:13, 100./ns)
hold off
xlabel('Number of frames in each split-half FC');ylabel('Identification accuracy (%)');
set(gca,'FontSize',16,'XTick',1:13,'XTickLabel',string([5:5:25,50:50:400]),...
    'XTickLabelRotation',45);
%% idenitifiability and age relations  
% (stats reported in paragraphs 2 and 3 of the result section:
% Functional connectome fingerprinting and reliability

Age_grad = NaN(158);
for j=1:158
    Age_grad(:,j) = [ag - ag(j)];
end
figure
imagesc(abs(Age_grad))
c1 = Corrs3(:); c1(eye(158)==1) = NaN; c2 = Age_grad(:);
[R,P] = corrcoef(abs(Age_grad(:)),c1,'Rows','pairwise');

gg=[];
for i=1:1000
    e = corrcoef(c1,abs(c2(randperm(24964))),'Rows','pairwise');
    gg = [gg; e(2,1)];
end
length(find(gg>R(2,1)))

same_corrs = Corrs3(find(eye(164)));
[r,p]=corr(same_corrs,ag);
idable = ones(164,1); idable(mms>=same_corrs)=0;
tts = table(idable,ag);
fitglm(tts,'idable ~ ag','Distribution','binomial')