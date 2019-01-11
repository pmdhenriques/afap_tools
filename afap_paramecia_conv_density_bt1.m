function [] = afap_paramecia_conv_density_bt1(datadir)
% Builds paramecia density plots for the first convergence bout during
% hunting routines
%
% Pedro Henriques

imstr = 'im';
[frid, x_px, y_px] = afaptdms_GetInfo(datadir, imstr);
npx = x_px*y_px;

load(fullfile(datadir,'aq'));

%%

hix = find([aq.btf.Convbout] ~= 0 & [aq.btf.stimbout] == 1);

% Use only bouts that are present in the movie
ix = [aq.btf(hix).ited];
ixx = aq.datainterp(ix,1);
ixxx = ismember(ixx,frid);
hix = hix(ixxx);

imst = zeros(x_px,y_px,length(hix),'uint8');
imend = zeros(x_px,y_px,length(hix),'uint8');

h = waitbar(0,'Initializing');
k = 1;
for i = hix
    waitbar(k/length(hix),h,'Processing')
    itst = aq.btf(i).itst;
    ited = aq.btf(i).ited;
    
    fst = aq.datainterp(itst,1);
    fend = aq.datainterp(ited,1);
    bst = afaptdms_GetFrLocation_ph(fst, frid);
    bend = afaptdms_GetFrLocation_ph(fend, frid);
    len = bend-bst;
    
    im = TDMS_readTDMSFile(fullfile(datadir, ...
        [imstr '.tdms']), ...
        'SUBSET_GET', [bst*npx bend*npx], ...
        'SUBSET_IS_LENGTH', false);
    
    im = reshape(im.data{3}(1:npx*len)', [x_px y_px len]);
    
    imst(:,:,k) = afap_im_rotate(im(:,:,1));
    imend(:,:,k) = afap_im_rotate(im(:,:,end));
    
    k = k+1;
end

%%

[optimizer,metric] = imregconfig('monomodal');

immean_st = mean(imst,3);
imreg_st = zeros(size(imst),'uint8');
immean_end = mean(imend,3);
imreg_end = zeros(size(imend),'uint8');

for i = 1:size(imst,3)
    waitbar(i/size(imst,3),h,'Processing')
    imreg_st(:,:,i) = imregister(imst(:,:,i),immean_st, ...
        'rigid',optimizer,metric,'PyramidLevels',4);
    imreg_end(:,:,i) = imregister(imend(:,:,i),immean_end, ...
        'rigid',optimizer,metric,'PyramidLevels',4);
end
close(h)

%%

Cst = struct;
Cend = struct;
for i = 1:size(imreg_st,3)
    Cst(i).cnt = afap_paramecia_segment(imreg_st(:,:,i));
    Cend(i).cnt = afap_paramecia_segment(imreg_end(:,:,i));
end

%%


C = cat(1,Cst.cnt);
[N,~] = hist3(C,{1:1:180,1:1:180});
N = N/size(C,1)/size(Cst,2);
N = imgaussfilt(N',4);
M = scaleif(mean(imreg_st,3),min(N(:)),max(N(:)));

C2 = cat(1,Cend.cnt);
[N2,~] = hist3(cat(1,C2),{1:1:180,1:1:180});
N2 = N2/size(C2,1)/size(Cend,2);
N2 = imgaussfilt(N2',4);
M2 = scaleif(mean(imreg_end,3),min(N2(:)),max(N2(:)));

fig = figure;
subplot(1,2,1)
imagesc(N+M); axis square off; title('Conv bout #1 start')
subplot(1,2,2)
imagesc(N2+M2); axis square off; title('Conv bout #1 end')
colormap(lbmap(256))

mkdir(fullfile(datadir,'figures'));
saveas(fig,fullfile(datadir,'figures','conv_bout_1_density.fig'));
close(fig);

%%

%     cp = [90.5666 66.1129];
%     hvec = [-0.0141 -18.5538 0];
%
%     for c = 1:size(C,1)
%         fish2target_vec = [C(c,1)-cp(1), C(c,2)-cp, 0];
%         C(c,3) = atan2d(hvec(1)*fish2target_vec(2)-hvec(2)*fish2target_vec(1), ...
%             hvec(1)*fish2target_vec(1)+hvec(2)*fish2target_vec(2));
%         C(c,4) = sqrt((C(c,1)-cp(1))^2 + ...
%             (C(c,2)-cp(2))^2);
%     end
%
%     for c = 1:size(C2,1)
%         fish2target_vec = [C2(c,1)-cp(1), C2(c,2)-cp, 0];
%         C2(c,3) = atan2d(hvec(1)*fish2target_vec(2)-hvec(2)*fish2target_vec(1), ...
%             hvec(1)*fish2target_vec(1)+hvec(2)*fish2target_vec(2));
%         C2(c,4) = sqrt((C2(c,1)-cp(1))^2 + ...
%             (C2(c,2)-cp(2))^2);
%     end
%
%     fig = figure;
%     polarhistogram(deg2rad(C(:,3)),'BinWidth',deg2rad(10),'Normalization','probability')
%     hold on
%     polarhistogram(deg2rad(C2(:,3)),'BinWidth',deg2rad(10),'Normalization','probability')
%     legend({'Start','End'},'Location','best')
%     saveas(fig,fullfile(datadir,'figures','conv_bout_1_polarhist.fig'));
%     close(fig);
save(fullfile(datadir,'figures','conv_start_analy.mat'), ...
    'Cst','Cend','imst','imend','imreg_st','imreg_end','C','C2');
end

%% Group stats
%
% folders = ctl_pre;
% h = waitbar(0,'Initializing');
% for f = 1:length(folders)
%     waitbar(f/length(folders),h,'Processing')
%     datadir = folders{f};
%     load(fullfile(datadir,'figures','conv_start_analy.mat'));
%
%     [N,~] = hist3(C(:,1:2),{1:2:180,1:2:180});
%     [N2,~] = hist3(C2(:,1:2),{1:2:180,1:2:180});
%     F.ctl.pre(f).N = N/size(Cst,2)/size(C,1);
%     F.ctl.pre(f).N2 = N2/size(Cend,2)/size(C2,1);
%
% end
% close(h)
%
% figure
% subplot(1,2,1)
% imagesc(imgaussfilt(prctile(cat(3,F.abl.post.N),75,3),4)'); axis square off
% title('Bout Start')
% subplot(1,2,2)
% imagesc(imgaussfilt(prctile(cat(3,F.abl.post.N2),75,3),4)'); axis square off
% title('Bout End')