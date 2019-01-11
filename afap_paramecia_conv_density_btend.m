function [] = afap_paramecia_conv_density_btend(datadir)
% Builds paramecia density plots for the last convergence bout during
% hunting routines
%
% Pedro Henriques

imstr = 'im';
[frid, x_px, y_px] = afaptdms_GetInfo(datadir, imstr);
npx = x_px*y_px;

load(fullfile(datadir,'aq'));

%%  Run the loop

% Find last convercence bout
hix = find([aq.btf.Convbout] ~= 0 & [diff([aq.btf.stimbout]) 0] < 0);

% Use only bouts that are present in the movie
ix = [aq.btf(hix).ited];
ixx = aq.datainterp(ix,1);
ixxx = ismember(ixx,frid);
hix = hix(ixxx);

imst = zeros(x_px,y_px,length(hix),'uint8');
imend = zeros(x_px,y_px,length(hix),'uint8');
imst2 = zeros(x_px,y_px,length(hix),'uint8');
imend2 = zeros(x_px,y_px,length(hix),'uint8');
imcend = zeros(x_px,y_px,length(hix),'uint8');

h = waitbar(0,'Initializing');
k = 1;
for i = hix
    waitbar(k/length(hix),h,'Processing')
    itst = aq.btf(i).itst;
    ited = aq.btf(i).ited;
    itst2 = aq.btf(i+1).itst;
    itend2 = aq.btf(i+1).ited;
    itcix = aq.btf(i).Convbout;
    itced = aq.huntep(itcix,2);
    
    fst = aq.datainterp(itst,1);
    fend = aq.datainterp(ited,1);
    fst2 = aq.datainterp(itst2,1);
    fend2 = aq.datainterp(itend2,1);
    fcend = aq.datainterp(itced,1);
    
    bst = afaptdms_GetFrLocation_ph(fst, frid);
    bend = afaptdms_GetFrLocation_ph(fend, frid);
    len = bend-bst;
    bst2 = afaptdms_GetFrLocation_ph(fst2, frid);
    bend2 = afaptdms_GetFrLocation_ph(fend2, frid);
    len2 = bend2-bst2;
    bcst = afaptdms_GetFrLocation_ph(fcend, frid);
    
    im = TDMS_readTDMSFile(fullfile(datadir, ...
        [imstr '.tdms']), ...
        'SUBSET_GET', [bst*npx bend*npx], ...
        'SUBSET_IS_LENGTH', false);
    
    im = reshape(im.data{3}(1:npx*len)', [x_px y_px len]);
    
    im2 = TDMS_readTDMSFile(fullfile(datadir, ...
        [imstr '.tdms']), ...
        'SUBSET_GET', [bst2*npx bend2*npx], ...
        'SUBSET_IS_LENGTH', false);
    
    im2 = reshape(im2.data{3}(1:npx*len2)', [x_px y_px len2]);
    
    imc = TDMS_readTDMSFile(fullfile(datadir, ...
        [imstr '.tdms']), ...
        'SUBSET_GET', [bcst*npx (bcst+1)*npx], ...
        'SUBSET_IS_LENGTH', false);
    
    imc = reshape(imc.data{3}(1:npx)', [x_px y_px]);
    
    imst(:,:,k) = afap_rotate_im(im(:,:,1));
    imend(:,:,k) = afap_rotate_im(im(:,:,end));
    imst2(:,:,k) = afap_rotate_im(im2(:,:,1));
    imend2(:,:,k) = afap_rotate_im(im2(:,:,end));
    imcend(:,:,k) = afap_rotate_im(imc);
    
    k = k+1;
end

%%  Rigid registrataion for more accurate results

[optimizer,metric] = imregconfig('monomodal');

imstmean = mean(imst,3);
imstreg = zeros(size(imst),'uint8');
imendmean = mean(imend,3);
imendreg = zeros(size(imend),'uint8');

imst2mean = mean(imst2,3);
imst2reg = zeros(size(imst2),'uint8');
imend2mean = mean(imend2,3);
imend2reg = zeros(size(imend2),'uint8');
imcendmean = mean(imcend,3);
imcendreg = zeros(size(imcend),'uint8');

for i = 1:size(imst,3)
    waitbar(i/size(imst,3),h,'Registering')
    imstreg(:,:,i) = imregister(imst(:,:,i),imstmean, ...
        'rigid',optimizer,metric,'PyramidLevels',4);
    imendreg(:,:,i) = imregister(imend(:,:,i),imendmean, ...
        'rigid',optimizer,metric,'PyramidLevels',4);
    imst2reg(:,:,i) = imregister(imst2(:,:,i),imst2mean, ...
        'rigid',optimizer,metric,'PyramidLevels',4);
    imend2reg(:,:,i) = imregister(imend2(:,:,i),imend2mean, ...
        'rigid',optimizer,metric,'PyramidLevels',4);
    imcendreg(:,:,i) = imregister(imcend(:,:,i),imcendmean, ...
        'rigid',optimizer,metric,'PyramidLevels',4);
end

%%  Segment paramecia

PD = struct;

for i = 1:size(imstreg,3)
    waitbar(i/size(imstreg,3),h,'Segmenting')
    PD.Cst(i).cnt = afap_paramecia_segment(imstreg(:,:,i));
    PD.Cend(i).cnt = afap_paramecia_segment(imendreg(:,:,i));
    PD.Cst2(i).cnt = afap_paramecia_segment(imst2reg(:,:,i));
    PD.Cend2(i).cnt = afap_paramecia_segment(imend2reg(:,:,i));
    PD.Ccend(i).cnt = afap_paramecia_segment(imcendreg(:,:,i));
end
close(h)

%%  Filter

C = cat(1,PD.Cst.cnt);
[N,~] = hist3(C,{1:1:x_px,1:1:x_px});
N = N/size(C,1)/size(PD.Cst,2);
N = imgaussfilt(N',4);
M = scaleif(mean(imstreg,3),min(N(:)),max(N(:)));

C2 = cat(1,PD.Cend.cnt);
[N2,~] = hist3(cat(1,C2),{1:1:x_px,1:1:x_px});
N2 = N2/size(C2,1)/size(PD.Cend,2);
N2 = imgaussfilt(N2',4);
M2 = scaleif(mean(imendreg,3),min(N2(:)),max(N2(:)));

C3 = cat(1,PD.Cst2.cnt);
[N3,~] = hist3(cat(1,C3),{1:1:x_px,1:1:x_px});
N3 = N3/size(C3,1)/size(PD.Cst2,2);
N3 = imgaussfilt(N3',4);
M3 = scaleif(mean(imst2reg,3),min(N3(:)),max(N3(:)));

C4 = cat(1,PD.Cend2.cnt);
[N4,~] = hist3(cat(1,C4),{1:1:x_px,1:1:x_px});
N4 = N4/size(C4,1)/size(PD.Cend2,2);
N4 = imgaussfilt(N4',4);
M4 = scaleif(mean(imend2reg,3),min(N4(:)),max(N4(:)));

C5 = cat(1,PD.Ccend.cnt);
[N5,~] = hist3(cat(1,C5),{1:1:x_px,1:1:x_px});
N5 = N5/size(C5,1)/size(PD.Ccend,2);
N5 = imgaussfilt(N5',4);
M5 = scaleif(mean(imcendreg,3),min(N5(:)),max(N5(:)));

%% Show results

fig = figure('Position',[200 100 1600 800]);
subplot(2,4,1)
imagesc(N+M); axis square off; title('Last conv bout start')
subplot(2,4,2)
imagesc(N2+M2); axis square off; title('Last conv bout end')
subplot(2,4,5)
imagesc(N3+M3); axis square off; title('Last conv bout+1 start')
subplot(2,4,6)
imagesc(N4+M4); axis square off; title('Last conv bout+1 end')
subplot(2,4,[3:4,7:8])
imagesc(N5+M5); axis square off; title('Conv end')
colormap(lbmap(256))

%%  Save figure, structure and images

if ~exist(fullfile(datadir,'figures'),'dir')
    mkdir(fullfile(datadir,'figures'));
end
saveas(fig,fullfile(datadir,'figures','conv_bout_end_density.fig'));
close(fig);
save(fullfile(datadir,'figures','conv_end_analy.mat'), ...
    'PD','imstreg','imendreg','imst2reg','imend2reg', ...
    'imcendreg');

end
