function [] = afap_paramecia_null_density(datadir)
% Used to compute the paramecia null density across the entire experiment
% 
% CODE NEEDS REVIEWING
%
% Pedro Henriques

imstr = 'im';
[frid, x_px, y_px] = afaptdms_GetInfo(datadir, imstr);
npx = x_px*y_px;

load(fullfile(datadir,'aq'),'aq');

%%
its = round(((aq.huntep(2:end,1)-aq.huntep(1:end-1,2))/2) ...
    +aq.huntep(1:end-1,2));

ix = aq.datainterp(its,1);
ixx = ismember(ix,frid);
its = its(ixx);

imn = zeros(x_px,y_px,length(its),'uint8');

h = waitbar(0,'Initializing');
for i = 1:length(its)
    waitbar(i/length(its),h,'Processing')
    it = its(i);
    fst = aq.datainterp(it,1);
    bst = afaptdms_GetFrLocation_ph(fst, frid);
    im = TDMS_readTDMSFile(fullfile(datadir, ...
        [imstr '.tdms']), ...
        'SUBSET_GET', [bst*npx (bst+1)*npx], ...
        'SUBSET_IS_LENGTH', false);
    im = reshape(im.data{3}(1:npx)', [x_px y_px]);
    imn(:,:,i) = afap_rotate_im(im);
end

%%

[optimizer,metric] = imregconfig('monomodal');

imnmean = mean(imn,3);
imnreg = zeros(size(imn),'uint8');
for i = 1:size(imn,3)
    waitbar(i/size(imn,3),h,'Registering')
    imnreg(:,:,i) = imregister(imn(:,:,i),imnmean, ...
        'rigid',optimizer,metric,'PyramidLevels',4);
end
close(h)

%%

for i = 1:size(imnreg,3)
    PDnull(i).cnt = afap_paramecia_segment(imnreg(:,:,i));
end

%%

Cn = cat(1,PDnull.cnt);
[Nn,~] = hist3(Cn,{1:1:x_px,1:1:x_px});
Nn = Nn/size(Cn,1)/size(PDnull,2);
Nn = imgaussfilt(Nn',4);
Gthr = prctile(Nn(:),99);

Dn = pdist(Cn);
Dn = squareform(Dn);
LDn = nanmean(exp(-Dn));
LDthr = prctile(LDn,99);

aq.Paramecia.GaussThr = Gthr;
aq.Paramecia.LocalDensityThr = LDthr;

save(fullfile(datadir,'aq'),'aq');

end