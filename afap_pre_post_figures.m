function [] = afap_pre_post_figures(datadir)
% Generates several behaviour figures related to pre and post conditions
% from the NI ablation experiment.
%
% Input is directory with pre a and post folders
%
% Pedro Henriques

if nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\NI project\Ablations', ...
        'Select directory to analyse');
end

%% Load structures

disp('Loading structures')
pre = load(fullfile(datadir,'pre','aq'),'aq');
post = load(fullfile(datadir,'post','aq'),'aq');
disp('Processing data')

if pre.aq.settings.analysisversion >= 180119
    stimtypcol = 2; stimlatcol = 3; stimstcol = 4; stimedcol = 5;
else
    stimtypcol = 1; stimlatcol = 2; stimstcol = 3; stimedcol = 4;
end

if ~exist(fullfile(datadir,'figures'),'dir')
    mkdir(fullfile(datadir,'figures'))
end

% Spatial distribution

sigma = 15;
ang=0:0.01:2*pi;
xp=(pre.aq.settings.CAMERApxpermm*17.5)*cos(ang);
yp=(pre.aq.settings.CAMERApxpermm*17.5)*sin(ang);

fig1 = figure('Name','Spatial distribution', ...
    'position',[400 400 900 400]);
colormap hot
subplot(1,2,1)
[n,~] = hist3(pre.aq.datainterp(:,2:3),{0:1:950 0:1:950});
B = imgaussfilt(n,sigma);
imagesc(B,[0 prctile(reshape(B,numel(B),1),99)]);
axis off square; hold on;
plot((pre.aq.settings.FOV(1)/2)+xp,(pre.aq.settings.FOV(2)/2)+yp, ...
    'Color', [1 1 1], 'LineWidth', 3);
title('Pre')

subplot(1,2,2)
[n,~] = hist3(post.aq.datainterp(:,2:3),{0:1:950 0:1:950});
B = imgaussfilt(n,sigma);
imagesc(B,[0 prctile(reshape(B,numel(B),1),99)]);
axis off square; hold on;
plot((post.aq.settings.FOV(1)/2)+xp,(post.aq.settings.FOV(2)/2)+yp, ...
    'Color', [1 1 1], 'LineWidth', 3);
title('Post 24h')

saveas(fig1,fullfile(datadir,'figures',[fig1.Name '.fig']));
close(fig1)

% Hunting

fig2 = figure('name','Vergence angle');
histogram(pre.aq.Verg,-20:1:120,'Normalization','probability'); hold on
histogram(post.aq.Verg,-20:1:120,'Normalization','probability'); hold off
legend pre post24
ylabel('Probability'); xlabel('Eye vergence angle');

saveas(fig2,fullfile(datadir,'figures',[fig2.Name '.fig']));
close(fig2)

%% 1st hunt bout dori

hbt1pre = find([0 diff([pre.aq.btf.Convbout])] > 0);
doripre = cumsum(pre.aq.dori_fix,'omitnan');
hbt1post = find([0 diff([post.aq.btf.Convbout])] > 0);
doripost = cumsum(post.aq.dori_fix,'omitnan');

fig3 = figure('name','1st hunt bout dori');
histogram(doripre([pre.aq.btf(hbt1pre).ited])-doripre([pre.aq.btf(hbt1pre).itst]), ...
    -100:3:100,'Normalization','probability'); hold on
histogram(doripost([post.aq.btf(hbt1post).ited])-doripost([post.aq.btf(hbt1post).itst]), ...
    -100:3:100,'Normalization','probability'); hold off
xlabel('Delta ori (deg)'); ylabel('probability');
legend({'pre','post 24h'});

saveas(fig3,fullfile(datadir,'figures',[fig3.Name '.fig']));
close(fig3)

% 2nd bout

hbt1pre = find([pre.aq.btf.Convbout] ~= 0 & [pre.aq.btf.stimbout] == 2);
doripre = cumsum(pre.aq.dori_fix,'omitnan');
hbt1post = find([post.aq.btf.Convbout] ~= 0 & [post.aq.btf.stimbout] == 2);
doripost = cumsum(post.aq.dori_fix,'omitnan');

fig3 = figure('name','2nd hunt bout dori');
histogram(doripre([pre.aq.btf(hbt1pre).ited])-doripre([pre.aq.btf(hbt1pre).itst]), ...
    -100:3:100,'Normalization','probability'); hold on
histogram(doripost([post.aq.btf(hbt1post).ited])-doripost([post.aq.btf(hbt1post).itst]), ...
    -100:3:100,'Normalization','probability'); hold off
xlabel('Delta ori (deg)'); ylabel('probability');
legend({'pre','post 24h'});

saveas(fig3,fullfile(datadir,'figures',[fig3.Name '.fig']));
close(fig3)


%% Mean hunting duration

fig4 = figure('name','Mean hunting duration');
histogram(pre.aq.huntep(:,3)/pre.aq.settings.camerarate, ...
    'BinWidth',0.2,'Normalization','probability'); hold on
histogram(post.aq.huntep(:,3)/post.aq.settings.camerarate, ...
    'BinWidth',0.2,'Normalization','probability');
xlabel('Hunting duration (s)'); ylabel('probability');
legend({'pre','post 24h'});

saveas(fig4,fullfile(datadir,'figures',[fig4.Name '.fig']));
close(fig4)

%% Hunting duration

hbtnpre = zeros(size(pre.aq.huntep,1),1);
for i = 1:size(hbtnpre,1)
    cst = pre.aq.huntep(i,1);
    ced = pre.aq.huntep(i,2);
    hbtnpre(i) = sum([pre.aq.btf.itst] >= cst-70 & [pre.aq.btf.ited] <= ced);
end
hbtnpost = zeros(size(post.aq.huntep,1),1);
for i = 1:size(hbtnpost,1)
    cst = post.aq.huntep(i,1);
    ced = post.aq.huntep(i,2);
    hbtnpost(i) = sum([post.aq.btf.itst] >= cst-70 & [post.aq.btf.ited] <= ced);
end

hbtfitpre = fit(pre.aq.huntep(:,3)./pre.aq.settings.camerarate,hbtnpre,'poly1');
hbtfitpost = fit(post.aq.huntep(:,3)./post.aq.settings.camerarate,hbtnpost,'poly1');

X = [pre.aq.huntep(:,3); post.aq.huntep(:,3)]./pre.aq.settings.camerarate;
Y = [hbtnpre; hbtnpost];
group = [zeros(length(hbtnpre),1); ones(length(hbtnpost),1)];

fig5 = figure('name','Hunting duration');
scatterhist(X,Y,'Group',group,'Kernel','on'); hold on;
plot(hbtfitpre,'b')
plot(hbtfitpost,'r')
legend pre post24
xlabel('Hunting bout duration (s)')
ylabel('Number of bouts')

saveas(fig5,fullfile(datadir,'figures',[fig5.Name '.fig']));
close(fig5)

%% Hunting bouts progression

nbts = 3;
X = NaN(nbts,2,6,2); % bt/mean+std/feature/pre+post
ix1 = find([pre.aq.btf.Convbout] & [pre.aq.btf.stimbout] == nbts);
ix1 = ix1-nbts+1;
ix2 = find([post.aq.btf.Convbout] & [post.aq.btf.stimbout] == nbts);
ix2 = ix2-nbts+1;

for b = 1:nbts
    bix1 = ix1+b-1;
    bix2 = ix2+b-1;
    
    X(b,1,1,1) = nanmean([pre.aq.btf(bix1).dur]*1000);
    X(b,2,1,1) = nanstd([pre.aq.btf(bix1).dur]*1000);
    X(b,1,1,2) = nanmean([post.aq.btf(bix2).dur]*1000);
    X(b,2,1,2) = nanstd([post.aq.btf(bix2).dur]*1000);
    
    X(b,1,2,1) = nanmean(abs([pre.aq.btf(bix1).max_angl]));
    X(b,2,2,1) = nanstd(abs([pre.aq.btf(bix1).max_angl]));
    X(b,1,2,2) = nanmean(abs([post.aq.btf(bix2).max_angl]));
    X(b,2,2,2) = nanstd(abs([post.aq.btf(bix2).max_angl]));
    
    X(b,1,3,1) = nanmean(abs([pre.aq.btf(bix1).max_vel]));
    X(b,2,3,1) = nanstd(abs([pre.aq.btf(bix1).max_vel]));
    X(b,1,3,2) = nanmean(abs([post.aq.btf(bix2).max_vel]));
    X(b,2,3,2) = nanstd(abs([post.aq.btf(bix2).max_vel]));
    
    X(b,1,4,1) = nanmean([pre.aq.btf(bix1).mean_TBF]);
    X(b,2,4,1) = nanstd([pre.aq.btf(bix1).mean_TBF]);
    X(b,1,4,2) = nanmean([post.aq.btf(bix2).mean_TBF]);
    X(b,2,4,2) = nanstd([post.aq.btf(bix2).mean_TBF]);
    
    X(b,1,5,1) = nanmean([pre.aq.btf(bix1).morphAI2]);
    X(b,2,5,1) = nanstd([pre.aq.btf(bix1).morphAI2]);
    X(b,1,5,2) = nanmean([post.aq.btf(bix2).morphAI2]);
    X(b,2,5,2) = nanstd([post.aq.btf(bix2).morphAI2]);
    
    % ignore last bout for IBI
    if ~isempty(bix2) && bix2(end) == length(post.aq.btf)
        bix2(end) = [];
    end
    if ~isempty(bix1) && bix1(end) == length(pre.aq.btf)
        bix1(end) = [];
    end
    X(b,1,6,1) = nanmean([pre.aq.btf(bix1+1).st]*1000-[pre.aq.btf(bix1).ed]*1000);
    X(b,2,6,1) = nanstd([pre.aq.btf(bix1+1).st]*1000-[pre.aq.btf(bix1).ed]*1000);
    X(b,1,6,2) = nanmean([post.aq.btf(bix2+1).st]*1000-[post.aq.btf(bix2).ed]*1000);
    X(b,2,6,2) = nanstd([post.aq.btf(bix2+1).st]*1000-[post.aq.btf(bix2).ed]*1000);
end

btftlab = {'Duration (ms)','Max angle (deg)','Max vel (deg/s)', ...
    'Mean TBF','MorphAI2','IBI (ms)'};
fig6 = figure('name','Hunting bout features', ...
    'position',[400 400 900 400]);

for f = 1:6
    subplot(2,3,f)
    errorbar(1:nbts,X(:,1,f,1),X(:,2,f,1),'o-'); hold on;
    errorbar(1:nbts,X(:,1,f,2),X(:,2,f,2),'o-');
    set(gca,'XTick',1:nbts); xlim([0 nbts+1])
    xlabel('Bout #'); ylabel(btftlab{f});
end

saveas(fig6,fullfile(datadir,'figures',[fig6.Name '.fig']));
close(fig6)

%% Paramecia curve

if exist(fullfile(datadir,'pre','Paramecia_count','counts.xlsx'),'file')
    parapre = xlsread(fullfile(datadir,'pre','Paramecia_count','counts'));
    parapost = xlsread(fullfile(datadir,'post','Paramecia_count','counts'));
    
    fig7 = figure('name','Paramecia curve');
    plot(parapre(:,1)/17/60,parapre(:,2),'o-'); hold on
    plot(parapost(:,1)/17/60,parapost(:,2),'o-');
    ylabel('# Paramecia'); xlabel('Time (s)');
    
    saveas(fig7,fullfile(datadir,'figures',[fig7.Name '.fig']));
    close(fig7)
end

%% Loom contrast-response rate

% loomcont = flip(unique(pre.aq.loomz(:,14)));
% L = zeros(5,5,2);
% for i = 1:length(loomcont)
%     rix = ~isnan(pre.aq.loomz(:,6));
%     six = pre.aq.loomz(:,14) == loomcont(i);
%     L(i,1,1) = sum(rix & six)/sum(six);
%     L(i,2,1) = sum(rix & six);
%     L(i,3,1) = sum(six);
%     L(i,4,1) = nanmean(pre.aq.loomz(rix & six,6)/pre.aq.settings.datarate);
%     L(i,5,1) = nanstd(pre.aq.loomz(rix & six,6)/pre.aq.settings.datarate);
%     
%     rix = ~isnan(post.aq.loomz(:,6));
%     six = post.aq.loomz(:,14) == loomcont(i);
%     L(i,1,2) = sum(rix & six)/sum(six);
%     L(i,2,2) = sum(rix & six);
%     L(i,3,2) = sum(six);
%     L(i,4,2) = nanmean(post.aq.loomz(rix & six,6)/post.aq.settings.datarate);
%     L(i,5,2) = nanstd(post.aq.loomz(rix & six,6)/post.aq.settings.datarate);
% end
% 
% x = 100-(mod(loomcont, 1000));
% 
% fig8 = figure('name','Loom contrast-response rate (~hunting)');
% plot(x,L(:,1,1),'o-'); hold on;
% plot(x,L(:,1,2),'o-')
% ylim([0 1]);
% ylabel('Response rate'); xlabel('Contrast (%)')
% legend({'pre','post 24h'});
% 
% saveas(fig8,fullfile(datadir,'figures',[fig8.Name '.fig']));
% close(fig8)

%% Loom response size

% fig9 = figure('name','Loom response size');
% subplot(1,2,1)
% boxplot(pre.aq.loomz(:,7),100-(mod(pre.aq.loomz(:,14),1000)), ...
%     'Notch','on'); hold on
% plot((100-(mod(pre.aq.loomz(:,14),1000)))./20,pre.aq.loomz(:,7),'*');
% ylim([0 110])
% ylabel('Stimulus size'); xlabel('Contrast (%)'); title('Pre')
% 
% subplot(1,2,2)
% boxplot(post.aq.loomz(:,7),100-(mod(post.aq.loomz(:,14),1000)), ...
%     'Colors','r','Notch','on'); hold on
% plot((100-(mod(post.aq.loomz(:,14),1000)))./20,post.aq.loomz(:,7),'r*');
% ylim([0 110])
% xlabel('Contrast (%)'); title('Post 24')
% 
% saveas(fig9,fullfile(datadir,'figures',[fig9.Name '.fig']));
% close(fig9)

%% Loom response laterality

LLpre = pre.aq.loomz(:,12) == 0;
LRpre = pre.aq.loomz(:,12) == 1;
LLpost = post.aq.loomz(:,12) == 0;
LRpost = post.aq.loomz(:,12) == 1;
Epre = ~isnan(pre.aq.loomz(:,7));
Epost = ~isnan(post.aq.loomz(:,7));

LLRpre = pre.aq.loomz(:,13) < 0;
LRRpre = pre.aq.loomz(:,13) > 0;
LLRpost = post.aq.loomz(:,13) < 0;
LRRpost = post.aq.loomz(:,13) > 0;

fig10 = figure('name','Loom response laterality');
bfig = bar([1,2,3,4,5], ...
    [sum(LLRpre & LLpre & Epre)/sum(Epre & LLpre), ...
    sum(LRRpre & LLpre & Epre)/sum(Epre & LLpre); ...
    sum(LLRpre & LRpre & Epre)/sum(Epre & LRpre), ...
    sum(LRRpre & LRpre & Epre)/sum(Epre & LRpre); ...
    [0 0]; ...
    sum(LLRpost & LLpost & Epost)/sum(Epost & LLpost), ...
    sum(LRRpost & LLpost & Epost)/sum(Epost & LLpost); ...
    sum(LLRpost & LRpost & Epost)/sum(Epost & LRpost), ...
    sum(LRRpost & LRpost & Epost)/sum(Epost & LRpost)], ...
    0.8);
xlim([0 6])
bfig(1).FaceColor = 'b'; bfig(2).FaceColor = 'r';
set(gca,'XTick',[1,2,4,5],'XTickLabel',{'Left (L/R)','Right (L/R)', ...
    'Left (L/R)','Right (L/R)'}, ...
    'YLim',[0,1]);
fix_xticklabels();
ylabel('Reponse rate (of max)');

saveas(fig10,fullfile(datadir,'figures',[fig10.Name '.fig']));
close(fig10)

%% Hunting loom size reponse

washuntpre = pre.aq.loomz(:,8) == 1 | pre.aq.loomz(:,9) == 1;
washuntpost = post.aq.loomz(:,8) == 1 | post.aq.loomz(:,9) == 1;

fig11 = figure('name','Hunting loom size reponse');
errorbar([nanmean(pre.aq.loomz(~washuntpre,7)) nanmean(pre.aq.loomz(washuntpre,7))], ...
    [nanstd(pre.aq.loomz(~washuntpre,7)) nanstd(pre.aq.loomz(washuntpre,7))],'-o');
hold on;
errorbar([nanmean(post.aq.loomz(~washuntpost,7)) nanmean(post.aq.loomz(washuntpost,7))], ...
    [nanstd(post.aq.loomz(~washuntpost,7)) nanstd(post.aq.loomz(washuntpost,7))],'-o');
set(gca,'XTick',1:2,'XTickLabel',{'Not hunting','Hunting'}); axis([0 3 0 110]);
ylabel('Stimulus size');
legend({'pre','post 24h'});

saveas(fig11,fullfile(datadir,'figures',[fig11.Name '.fig']));
close(fig11)

%% Loom bout features

loomrspixpre = ~isnan(pre.aq.loomz(:,6));
loomrspstpre = pre.aq.loomz(loomrspixpre,1)+pre.aq.loomz(loomrspixpre,3);
loomrspixpost = ~isnan(post.aq.loomz(:,6));
loomrspstpost = post.aq.loomz(loomrspixpost,1)+post.aq.loomz(loomrspixpost,3);

loomrspbtspre = zeros(length(loomrspstpre),1);
for i = 1:length(loomrspstpre)
    loomrspbtspre(i) = findnearest([pre.aq.btf.itst],loomrspstpre(i),1);
end
loomrspbtspost = zeros(length(loomrspstpost),1);
for i = 1:length(loomrspstpost)
    loomrspbtspost(i) = findnearest([post.aq.btf.itst],loomrspstpost(i),1);
end

fig12 = figure('name','Loom bout features', ...
    'position',[400 300 900 600]);

ypre = abs(doripre([pre.aq.btf(loomrspbtspre).ited]) - ...
    doripre([pre.aq.btf(loomrspbtspre).itst]));
ypost = abs(doripost([post.aq.btf(loomrspbtspost).ited]) - ...
    doripost([post.aq.btf(loomrspbtspost).itst]));
subplot(2,3,1)
violinplot(padcat(ypre,ypost));
set(gca,'XTickLabel',{'Pre','Post 24h'});
ylabel('Delta ori (deg)');

ypre = [pre.aq.btf(loomrspbtspre).dur]*1000;
ypost = [post.aq.btf(loomrspbtspost).dur]*1000;
subplot(2,3,2)
violinplot(padcat(ypre,ypost)');
set(gca,'XTickLabel',{'Pre','Post 24h'});
ylabel('Bout duration');

ypre = abs([pre.aq.btf(loomrspbtspre).max_angl]);
ypost = abs([post.aq.btf(loomrspbtspost).max_angl]);
subplot(2,3,3)
violinplot(padcat(ypre,ypost)');
set(gca,'XTickLabel',{'Pre','Post 24h'});
ylabel('Bout max angle');

ypre = abs([pre.aq.btf(loomrspbtspre).max_vel]);
ypost = abs([post.aq.btf(loomrspbtspost).max_vel]);
subplot(2,3,4)
violinplot(padcat(ypre,ypost)');
set(gca,'XTickLabel',{'Pre','Post 24h'});
ylabel('Bout max velocity');

ypre = abs([pre.aq.btf(loomrspbtspre).mean_TBF]);
ypost = abs([post.aq.btf(loomrspbtspost).mean_TBF]);
subplot(2,3,5)
violinplot(padcat(ypre,ypost)');
set(gca,'XTickLabel',{'Pre','Post 24h'});
ylabel('Bout mean TBF');

ypre = abs([pre.aq.btf(loomrspbtspre).morphAI2]);
ypost = abs([post.aq.btf(loomrspbtspost).morphAI2]);
subplot(2,3,6)
violinplot(padcat(ypre,ypost)');
set(gca,'XTickLabel',{'Pre','Post 24h'});
ylabel('MorphAI2');

saveas(fig12,fullfile(datadir,'figures',[fig12.Name '.fig']));
close(fig12)

%% Eye convergence

for i = 1:size(pre.aq.huntep,1)
    ixb = find([pre.aq.btf.Convbout] == i & [0 diff([pre.aq.btf.Convbout])] >= 1);
    if ~isempty(ixb)
        ixb = ixb(1);
        pre.aq.huntep(i,9) = pre.aq.cumori(pre.aq.btf(ixb).ited)-pre.aq.cumori(pre.aq.btf(ixb).itst);
    else
        pre.aq.huntep(i,9) = 0;
    end
end

for i = 1:size(post.aq.huntep,1)
    ixb = find([post.aq.btf.Convbout] == i & [0 diff([post.aq.btf.Convbout])] >= 1);
    if ~isempty(ixb)
        ixb = ixb(1);
        post.aq.huntep(i,9) = post.aq.cumori(post.aq.btf(ixb).ited)-post.aq.cumori(post.aq.btf(ixb).itst);
    else
        post.aq.huntep(i,9) = 0;
    end
end

lixpre = pre.aq.huntep(:,9) < 0;
lixpost = post.aq.huntep(:,9) < 0;

warning('off','all')
Cpre = NaN(size(pre.aq.huntep,1),601,2);
Cpost = NaN(size(post.aq.huntep,1),601,2);
for i = 1:size(Cpre,1)
    ix = pre.aq.huntep(i,1);
    if ix-300 > 0 && ix+300 < length(pre.aq.Rangl)
        L = pre.aq.Langl(ix-300:ix+300);
        R = pre.aq.Rangl(ix-300:ix+300);
        Lv = abs(diff(smooth(L,140,'loess')));
        Rv = abs(diff(smooth(R,140,'loess')));
        [Lpks,Llocs] = findpeaks(Lv,'MinPeakHeight',0.25);
        [Rpks,Rlocs] = findpeaks(Rv,'MinPeakHeight',0.25);
        Lvix = find(abs(Llocs-300) <= 100);
        Rvix = find(abs(Rlocs-300) <= 100);
        if (~isempty(Lvix) && ~isempty(Rvix) && max(Lpks(Lvix)) >= max(Rpks(Rvix))) || ...
                (~isempty(Lvix) && isempty(Rvix))
            [~, mix] = max(Lpks(Lvix));
            ix = ix+(Llocs(Lvix(mix))-300);
        elseif (~isempty(Lvix) && ~isempty(Rvix) && max(Lpks(Lvix)) <= max(Rpks(Rvix))) || ...
                (isempty(Lvix) && ~isempty(Rvix))
            [~, mix] = max(Rpks(Rvix));
            ix = ix+(Rlocs(Rvix(mix))-300);
        end
        if ix-300 > 0 && ix+300 < length(pre.aq.Rangl)
            Cpre(i,:,1) = pre.aq.Langl(ix-300:ix+300);
            Cpre(i,:,2) = pre.aq.Rangl(ix-300:ix+300);
        end
    end
end

for i = 1:size(Cpost,1)
    ix = post.aq.huntep(i,1);
    if ix-300 > 0 && ix+300 < length(post.aq.Rangl)
        L = post.aq.Langl(ix-300:ix+300);
        R = post.aq.Rangl(ix-300:ix+300);
        Lv = abs(diff(smooth(L,140,'loess')));
        Rv = abs(diff(smooth(R,140,'loess')));
        [Lpks,Llocs] = findpeaks(Lv,'MinPeakHeight',0.25);
        [Rpks,Rlocs] = findpeaks(Rv,'MinPeakHeight',0.25);
        Lvix = find(abs(Llocs-300) <= 100);
        Rvix = find(abs(Rlocs-300) <= 100);
        if (~isempty(Lvix) && ~isempty(Rvix) && max(Lpks(Lvix)) >= max(Rpks(Rvix))) || ...
                (~isempty(Lvix) && isempty(Rvix))
            [~, mix] = max(Lpks(Lvix));
            ix = ix+(Llocs(Lvix(mix))-300);
        elseif (~isempty(Lvix) && ~isempty(Rvix) && max(Lpks(Lvix)) <= max(Rpks(Rvix))) || ...
                (isempty(Lvix) && ~isempty(Rvix))
            [~, mix] = max(Rpks(Rvix));
            ix = ix+(Rlocs(Rvix(mix))-300);
        end
        if ix-300 > 0 && ix+300 < length(post.aq.Rangl)
            Cpost(i,:,1) = post.aq.Langl(ix-300:ix+300);
            Cpost(i,:,2) = post.aq.Rangl(ix-300:ix+300);
        end
    end
end
warning('on','all')

fig13 = figure('name','Eye convergence left turns');
shadedErrorBar(1:601, ...
    nanmean(Cpre(lixpre,:,1)-nanmean(Cpre(lixpre,100:200,1),2)) + ...
    nanmean(nanmean(Cpre(lixpre,100:200,1),2)), ...
    nansem(Cpre(lixpre,:,1)-nanmean(Cpre(lixpre,100:200,1),2)),'-b')
hold on
shadedErrorBar(1:601, ...
    nanmean(Cpre(lixpre,:,2)-nanmean(Cpre(lixpre,100:200,2),2)) + ...
    nanmean(nanmean(Cpre(lixpre,100:200,2),2)), ...
    nansem(Cpre(lixpre,:,2)-nanmean(Cpre(lixpre,100:200,2),2)),'-r')
shadedErrorBar(1:601, ...
    nanmean(Cpost(lixpost,:,1)-nanmean(Cpost(lixpost,100:200,1),2)) + ...
    nanmean(nanmean(Cpost(lixpost,100:200,1),2)), ...
    nansem(Cpost(lixpost,:,1)-nanmean(Cpost(lixpost,100:200,1),2)),'-c')
shadedErrorBar(1:601, ...
    nanmean(Cpost(lixpost,:,2)-nanmean(Cpost(lixpost,100:200,2),2)) + ...
    nanmean(nanmean(Cpost(lixpost,100:200,2),2)), ...
    nansem(Cpost(lixpost,:,2)-nanmean(Cpost(lixpost,100:200,2),2)),'-m')
ylabel('Eye angle (deg)'); xlabel('Time (frames)'); title('Left turns');
line([301 301],[-50 50],'LineStyle','--','Color','k')

fig14 = figure('name','Eye convergence right turns');
shadedErrorBar(1:601, ...
    nanmean(Cpre(~lixpre,:,1)-nanmean(Cpre(~lixpre,100:200,1),2)) + ...
    nanmean(nanmean(Cpre(~lixpre,100:200,1),2)), ...
    nansem(Cpre(~lixpre,:,1)-nanmean(Cpre(~lixpre,100:200,1),2)),'-b')
hold on
shadedErrorBar(1:601, ...
    nanmean(Cpre(~lixpre,:,2)-nanmean(Cpre(~lixpre,100:200,2),2)) + ...
    nanmean(nanmean(Cpre(~lixpre,100:200,2),2)), ...
    nansem(Cpre(~lixpre,:,2)-nanmean(Cpre(~lixpre,100:200,2),2)),'-r')
shadedErrorBar(1:601, ...
    nanmean(Cpost(~lixpost,:,1)-nanmean(Cpost(~lixpost,100:200,1),2)) + ...
    nanmean(nanmean(Cpost(lixpost,100:200,1),2)), ...
    nansem(Cpost(~lixpost,:,1)-nanmean(Cpost(~lixpost,100:200,1),2)),'-c')
shadedErrorBar(1:601, ...
    nanmean(Cpost(~lixpost,:,2)-nanmean(Cpost(~lixpost,100:200,2),2)) + ...
    nanmean(nanmean(Cpost(~lixpost,100:200,2),2)), ...
    nansem(Cpost(~lixpost,:,2)-nanmean(Cpost(~lixpost,100:200,2),2)),'-m')
ylabel('Eye angle (deg)'); xlabel('Time (frames)'); title('Left turns');
line([301 301],[-50 50],'LineStyle','--','Color','k')
ylabel('Eye angle (deg)'); xlabel('Time (frames)'); title('Right turns');
line([301 301],[-50 50],'LineStyle','--','Color','k')

saveas(fig13,fullfile(datadir,'figures',[fig13.Name '.fig']));
saveas(fig14,fullfile(datadir,'figures',[fig14.Name '.fig']));
close(fig13); close(fig14);

%% Hunting bout features

nbts = 3;
for b = 1:nbts
    
    hixpre = find([pre.aq.btf.Convbout] ~= 0 & [pre.aq.btf.stimbout] == b);
    hixpost = find([post.aq.btf.Convbout] ~= 0 & [post.aq.btf.stimbout] == b);        
    
    ypre = pre.aq.cumori(([pre.aq.btf(hixpre).ited])) - ...
        pre.aq.cumori([pre.aq.btf(hixpre).itst]);
    ypost = post.aq.cumori([post.aq.btf(hixpost).ited]) - ...
        post.aq.cumori([post.aq.btf(hixpost).itst]);
    
    
    fig13 = figure('name',sprintf('Hunting bout %d features',b), ...
        'position',[400 300 900 600]);
    subplot(2,3,1)
    violinplot(padcat(ypre,ypost));
    set(gca,'XTickLabel',{'Pre','Post 24h'});
    ylabel('Delta ori (deg)');
    
    ypre = [pre.aq.btf(hixpre).dur]*1000;
    ypost = [post.aq.btf(hixpost).dur]*1000;
    subplot(2,3,2)
    violinplot(padcat(ypre,ypost)');
    set(gca,'XTickLabel',{'Pre','Post 24h'});
    ylabel('Bout duration');
    
    ypre = [pre.aq.btf(hixpre).max_angl];
    ypost =[post.aq.btf(hixpost).max_angl];
    subplot(2,3,3)
    violinplot(padcat(ypre,ypost)');
    set(gca,'XTickLabel',{'Pre','Post 24h'});
    ylabel('Bout max angle');
    
    ypre = [pre.aq.btf(hixpre).max_vel];
    ypost = [post.aq.btf(hixpost).max_vel];
    subplot(2,3,4)
    violinplot(padcat(ypre,ypost)');
    set(gca,'XTickLabel',{'Pre','Post 24h'});
    ylabel('Bout max velocity');
    
    ypre = [pre.aq.btf(hixpre).mean_TBF];
    ypost = [post.aq.btf(hixpost).mean_TBF];
    subplot(2,3,5)
    violinplot(padcat(ypre,ypost)');
    set(gca,'XTickLabel',{'Pre','Post 24h'});
    ylabel('Bout mean TBF');
    
    ypre = [pre.aq.btf(hixpre).morphAI2];
    ypost = [post.aq.btf(hixpost).morphAI2];
    subplot(2,3,6)
    violinplot(padcat(ypre,ypost)');
    set(gca,'XTickLabel',{'Pre','Post 24h'});
    ylabel('MorphAI2');
    
    saveas(fig13,fullfile(datadir,'figures',[fig13.Name '.fig']));
    close(fig13)
end

%% OMR

omrpre = pre.aq.visstim(:,stimtypcol) == 3;
omrLpre = pre.aq.visstim(:,stimlatcol) == 0;
maxitLpre = max(diff(pre.aq.visstim(omrpre & omrLpre,stimstcol:stimedcol),1,2));
maxitRpre = max(diff(pre.aq.visstim(omrpre & ~omrLpre,stimstcol:stimedcol),1,2));
OLpre = NaN(sum(omrpre & omrLpre),maxitLpre);
ORpre = NaN(sum(omrpre & ~omrLpre),maxitRpre);

ix = find(omrpre & omrLpre);
for i = 1:length(ix)
    x = pre.aq.cumori(pre.aq.visstim(ix(i),stimstcol):pre.aq.visstim(ix(i),stimedcol));
    x = x-x(1);
    OLpre(i,1:length(x)) = x;
end
ix = find(omrpre & ~omrLpre);
for i = 1:length(ix)
    x = pre.aq.cumori(pre.aq.visstim(ix(i),stimstcol):pre.aq.visstim(ix(i),stimedcol));
    x = x-x(1);
    ORpre(i,1:length(x)) = x;
end

omrpost = post.aq.visstim(:,stimtypcol) == 3;
omrLpost = post.aq.visstim(:,stimlatcol) == 0;
maxitLpost = max(diff(post.aq.visstim(omrpost & omrLpost,stimstcol:stimedcol),1,2));
maxitRpost = max(diff(post.aq.visstim(omrpost & ~omrLpost,stimstcol:stimedcol),1,2));
OLpost = NaN(sum(omrpost & omrLpost),maxitLpost);
ORpost = NaN(sum(omrpost & ~omrLpost),maxitRpost);

ix = find(omrpost & omrLpost);
for i = 1:length(ix)
    x = post.aq.cumori(post.aq.visstim(ix(i),stimstcol):post.aq.visstim(ix(i),stimedcol));
    x = x-x(1);
    OLpost(i,1:length(x)) = x;
end
ix = find(omrpost & ~omrLpost);
for i = 1:length(ix)
    x = post.aq.cumori(post.aq.visstim(ix(i),stimstcol):post.aq.visstim(ix(i),stimedcol));
    x = x-x(1);
    ORpost(i,1:length(x)) = x;
end

fig14 = figure('Name', 'OMR responses', ...
    'position',[400 300 900 500]);

subplot(1,2,1)
plot(OLpre','b'); hold on
ix = find(any(isnan(OLpre)));
if isempty(ix); ix = size(OLpre,2); end
plot(nanmean(OLpre(:,1:ix(1))),'b','LineWidth',5);
plot(ORpre','r'); hold on
ix = find(any(isnan(ORpre)));
if isempty(ix); ix = size(ORpre,2); end
plot(nanmean(ORpre(:,1:ix(1))),'R','LineWidth',5);
axis([0 5000 -500 500])
line([0 5000],[0 0],'LineStyle','--','Color','k')
ylabel('Cumulative angle (deg)')
xlabel('Frames')
title('Pre')

subplot(1,2,2)
plot(OLpost','b'); hold on
ix = find(any(isnan(OLpost)));
if isempty(ix); ix = size(OLpost,2); end
plot(nanmean(OLpost(:,1:ix(1))),'b','LineWidth',5);
plot(ORpost','r'); hold on
ix = find(any(isnan(ORpost)));
if isempty(ix); ix = size(ORpost,2); end
plot(nanmean(ORpost(:,1:ix(1))),'R','LineWidth',5);
axis([0 5000 -500 500])
line([0 5000],[0 0],'LineStyle','--','Color','k')
xlabel('Frames')
title('Pre')
axis([0 5000 -500 500])
xlabel('Frames')
title('Post')

saveas(fig14,fullfile(datadir,'figures',[fig14.Name '.fig']));
close(fig14)


%% Taps

tapzpre = pre.aq.tapz;
tapzpost = post.aq.tapz;

fig15 = figure('Name','Tap responses', ...
    'position',[400 300 900 500]);
subplot(1,3,1)
denompre = size(tapzpre, 1);
numerpre = sum(~isnan(tapzpre(:,6)));
denompost = size(tapzpost, 1);
numerpost = sum(~isnan(tapzpost(:,6)));
bar(1,numerpre/denompre, .6, 'FaceColor', 'b'); hold on;
bar(2,numerpost/denompost, .6, 'FaceColor', 'r')
set(gca, 'XTick', 1:2, 'XTickLabels',{'Pre','Post'}, ...
    'XLim',[0 3],'YLim',[0 1],'Box', 'off')
ylabel('Response Rate')

subplot(1,3,2)
lmeanpre = 1000.*nanmean(tapzpre(:,6))./pre.aq.settings.datarate;
lstdpre = 1000.*nanstd(tapzpre(:,6))./pre.aq.settings.datarate;
lmeanpost = 1000.*nanmean(tapzpost(:,6))./post.aq.settings.datarate;
lstdpost = 1000.*nanstd(tapzpost(:,6))./post.aq.settings.datarate;
bar(1,lmeanpre,0.6,'FaceColor','b'); hold on;
errorbar(1, lmeanpre, lstdpre, 'LineWidth',2, 'Color','k')
bar(2,lmeanpost,0.6,'FaceColor','r');
errorbar(2, lmeanpost, lstdpost, 'LineWidth',2, 'Color','k')
set(gca, 'XTick',1:2, 'XTickLabels',{'Pre','Post'},'Box', 'off','XLim',[0 3])
ylabel('Latency [ms]')

subplot(1,3,3)
rtrpre = ~isnan(tapzpre(:,6));
nLpre = sum(tapzpre(rtrpre,13)<0)/sum(rtrpre);
nRpre = sum(tapzpre(rtrpre,13)>0)/sum(rtrpre);
rtrpost = ~isnan(tapzpost(:,6));
nLpost = sum(tapzpost(rtrpost,13)<0)/sum(rtrpost);
nRpost = sum(tapzpost(rtrpost,13)>0)/sum(rtrpost);

bar(1:2,[nLpre nRpre], .6,'FaceColor','b'); hold on;
bar(4:5,[nLpost nRpost], .6,'FaceColor','r'); hold off
set(gca, 'XTick', [1.5,4.5], 'XTickLabels',{'Pre (L/R)','Post (L/R)'}, ...
    'Box','off','YLim',[0 1],'XLim',[0 6])
ylabel('Escape rate')
xlabel('Laterality')


saveas(fig15,fullfile(datadir,'figures',[fig15.Name '.fig']));
close(fig15)
end