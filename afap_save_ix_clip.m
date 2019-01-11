function [] = afap_save_ix_clip(datadir,vis,lat,ix)

% vis = 4; % 1 = loom; 2 = tap; 3 = grat; 4 = hunt
% lat = 1;
% ix = 242;

fstoffsethi = 200;
fendoffsethi = 200;
imthr = [0 60];
highres = 1;

imstr = 'imbkd';
[frid, x_px, y_px] = afaptdms_GetInfo(datadir, imstr);
npx = x_px*y_px;

[frid2, x_px2, y_px2] = afaptdms_GetInfo(datadir, 'im');
npx2 = x_px2*y_px2;

% Load structures
disp('Loading structures');
load(fullfile(datadir,'aq'),'aq');
load(fullfile(datadir,'objects'));
t = matfile(fullfile(datadir,'tracks'));

labels = {'loom','tap','grat','hunt'};

if vis == 1
    fst = aq.datainterp(aq.loomz(ix,1)+aq.loomz(ix,6)-fstoffsethi,1);
    fend = aq.datainterp(aq.loomz(ix,2)+fendoffsethi,1);
elseif vis == 2
    fst = aq.datainterp(aq.tapz(ix,1)+aq.tapz(ix,6)-fstoffsethi,1);
    fend = aq.datainterp(aq.tapz(ix,2)+fendoffsethi,1);
elseif vis == 3
    fixs = find(aq.visstim(:,1) == vis & aq.visstim(:,2) == lat);
    fst = aq.datainterp(aq.visstim(fixs(ix),3),1);
    fend = aq.datainterp(aq.visstim(fixs(ix),4),1);
elseif vis == 4
    fst = aq.datainterp(aq.huntep(ix,1)-fstoffsethi,1);
    fend = aq.datainterp(aq.huntep(ix,2)+fendoffsethi,1);
end

% Load image frames

bst = afaptdms_GetFrLocation_ph(fst, frid);
bend = afaptdms_GetFrLocation_ph(fend, frid);
imcfrids = frid(bst:bend-1);
len = bend-bst;

im = TDMS_readTDMSFile(fullfile(datadir, ...
    [imstr '.tdms']), ...
    'SUBSET_GET', [bst*npx bend*npx], ...
    'SUBSET_IS_LENGTH', false);

im = reshape(im.data{3}(1:npx*len)', [x_px y_px len]);

% load tracks for current image frames
tracks = t.tracks(:,bst:bst+len);

fst2 = aq.data(aq.huntep(ix,1)-fstoffsethi,1);
fend2 = aq.data(aq.huntep(ix,2)+fendoffsethi,1);

bst2 = afaptdms_GetFrLocation_ph(fst2, frid2);
bend2 = afaptdms_GetFrLocation_ph(fend2, frid2);
im2cfrids = frid2(bst2:bend2-1);
len2 = bend2-bst2;

im2 = TDMS_readTDMSFile(fullfile(datadir, ...
    ['im' '.tdms']), ...
    'SUBSET_GET', [bst2*npx2 bend2*npx2], ...
    'SUBSET_IS_LENGTH', false);

im2 = reshape(im2.data{3}(1:npx2*len2)', [x_px2 y_px2 len2]);

%% Show the clip
writerObj = VideoWriter( ...
    fullfile(datadir,sprintf('%s_%d_hires%d',labels{vis},ix,highres)), ...
    'Motion JPEG AVI');
open(writerObj)

% Show high resolution video
fig = figure; colormap hot;
for i = 1:len2
    
    % Show image
    imagesc(im2(:,:,i),imthr);
    hold on
    
    cfrid2 = im2cfrids(i);
    frame = find(aq.datainterp(:,1) == cfrid2);
    fishcent = [aq.datainterp(frame,3) aq.datainterp(frame,2)];
    imcfridix = findnearest(cfrid2,imcfrids);
    imcfridix = imcfridix(1);
    bix = bst+imcfridix-1;
    
    % Find tracks in current frame
    trackids = find(tracks(:,imcfridix) ~= 0);
    tracksobjix = tracks(trackids,imcfridix);
    
    % Choose objects to draw
    if ~isempty(trackids)
        objcents = cat(1,objects(bix).centroid(tracksobjix,:));
        objcents2draw = ...
            objcents(:,1) >= fishcent(1)-x_px2/2 & ...
            objcents(:,1) <= fishcent(1)+x_px2/2 & ...
            objcents(:,2) >= fishcent(2)-y_px2/2 & ...
            objcents(:,2) <= fishcent(2)+y_px2/2;
        if any(objcents2draw)
            objcents_corr = [objcents(objcents2draw,1)+(x_px2/2)-fishcent(1), ...
                objcents(objcents2draw,2)+(y_px2/2)-fishcent(2)];
            text(objcents_corr(:,1)+5, ...
                objcents_corr(:,2), ...
                num2str(trackids(objcents2draw)), 'Color','w')
        end
    end
    
    text(10,170,sprintf('%.0f/%.0f ms', ...
        (i-fstoffsethi)/aq.settings.camerarate*1000, ...
        (len2-fstoffsethi)/aq.settings.camerarate*1000), ...
        'Color','w');
    hold off
    axis square off
    drawnow
    
    thisaviframe = getframe(gca);
    writeVideo(writerObj,thisaviframe);
end
close(writerObj)
close(fig)
end
