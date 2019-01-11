function [] = afap_im_save_clip(datadir,highres,vis,lat,ix,plotptracks,plotftracks,motioncorr,foffset,aq)
% Extracts and saves a high resolution clip of the desired index of a
% certain visual stimulis id and plots several eye, tail and paramecia
% features if prompt.
% User may also use names of people as paramecia IDs
%
% vis: 1 = loom; 2 = tap; 3 = grat; 4 = hunt
% lat: 0 = left; 1 = right;
%
% Ex: vis = 1; lat = 1; ix = 5 -- 5th right looming stimulus
%
% plotptracks - plot paramecia tracks
% plotftracks - plot fish eye and tail features
% motioncorr - motion correct the movie
% foffset - Frame offset from begining and end of routine
%
% Track names uses dataset from:
% https://www.kaggle.com/kaggle/us-baby-names/downloads/NationalNames.csv
%
% Pedro Henriques, Oct 2017
% Modified, Oct 2018
% Modified, Nov 2018

plotptext = 1; % Plot paramecia tracks as text?
usenames = 0;
bodythr = 0.01; % Fish body threshold
eyethr = 0.2;   % Fish eyes threshold
mrksz = 15; % Paramecia tracks marker size
vframerateratio = 0.1;  % Video framerate ratio
tracks2exc = [57796; 57149; 57176; 57209; 57213; 88697; 88703; 88714; 88719;
    57191; 57787]; % Tracks to exclude
tfiltersz = 7;

if usenames
    names = readtable("\\128.40.168.141\data\Pedro\Other\NationalNames.csv");
    names = unique(names.Name);
    nix = randperm(length(names));
    names = names(nix);
end

labels = {'loom','tap','grat','hunt'};

[frid1,x_px1,y_px1] = afaptdms_GetInfo(datadir,'imbkd');
npx1 = x_px1*y_px1;
[frid2, x_px2, y_px2] = afaptdms_GetInfo(datadir, 'im');
npx2 = x_px2*y_px2;

if ~highres
    plotftracks = 0;
    motioncorr = 0;
    imthr = [0 20];
    bcol = 'w';
    [XG,YG] = meshgrid(1:1:x_px1,1:y_px1);
    frmrt = 17.5;
else
    imthr = [0 40];
    bcol = 'w';
    [XG,YG] = meshgrid(1:1:x_px2,1:y_px2);
    frmrt = 700;
end

vframerate = frmrt*vframerateratio;

% Load structures
disp('Loading structures');
if nargin < 10
    load(fullfile(datadir,'aq'),'aq');
end
load(fullfile(datadir,'objects'),'objects');
if plotptracks
    t = matfile(fullfile(datadir,'tracks'));
end

if nargin < 9
    foffset = 0; % Time offset in seconds
end

% Arena boundaries
ang=0:0.01:2*pi;
xp=((aq.settings.FOV(1)/2)-5)*cos(ang);
yp=((aq.settings.FOV(2)/2)-5)*sin(ang);

if vis == 1
    fst = aq.datainterp(aq.loomz(ix,1)+aq.loomz(ix,6)-foffset,1);
    fend = aq.datainterp(aq.loomz(ix,2)+foffset,1);
elseif vis == 2
    fst = aq.datainterp(aq.tapz(ix,1)+aq.tapz(ix,6)-foffset,1);
    fend = aq.datainterp(aq.tapz(ix,2)+foffset,1);
elseif vis == 3
    fixs = find(aq.visstim(:,1) == vis & aq.visstim(:,2) == lat);
    fst = aq.datainterp(aq.visstim(fixs(ix),3),1);
    fend = aq.datainterp(aq.visstim(fixs(ix),4),1);
elseif vis == 4
    fst = aq.datainterp(aq.huntep(ix,1)-foffset,1);
    fend = aq.datainterp(aq.huntep(ix,2)+foffset,1);
end

% Load image frames

bst1 = afaptdms_GetFrLocation_ph(fst, frid1);
bend1 = afaptdms_GetFrLocation_ph(fend, frid1);
im1cfrids = frid1(bst1:bend1-1);
len1 = bend1-bst1;

if plotptracks
    % load tracks for current image frames
    tracks = t.tracks(:,bst1:bst1+len1-1);
    tlen = size(tracks,2);
end

bst2 = afaptdms_GetFrLocation_ph(fst, frid2);
bend2 = afaptdms_GetFrLocation_ph(fend, frid2);
im2cfrids = frid2(bst2:bend2-1);
len2 = bend2-bst2;

if highres
    im = TDMS_readTDMSFile(fullfile(datadir, ...
        ['im' '.tdms']), ...
        'SUBSET_GET', [bst2*npx2 bend2*npx2], ...
        'SUBSET_IS_LENGTH', false);
    im = reshape(im.data{3}(1:npx2*len2)', [x_px2 y_px2 len2]);
    if motioncorr
        [im,offset] = afap_im_motion_correct(im);
    else
        offset = zeros(len2,2);
    end
    len = len2;
else
    im = TDMS_readTDMSFile(fullfile(datadir, ...
        ['imbkd' '.tdms']), ...
        'SUBSET_GET', [bst1*npx1 bend1*npx1], ...
        'SUBSET_IS_LENGTH', false);
    im = reshape(im.data{3}(1:npx1*len1)', [x_px1 y_px1 len1]);
    offset = zeros(len1,2);
    len = len1;
end

%% Show the clip

mn = 1;
savepath = fullfile(datadir,sprintf('%s_%d_hires%d_ft%d_pt%d_%0.2fx_%d.avi', ...
    labels{vis},ix,highres,plotftracks,plotptracks,vframerateratio,mn));
while exist(savepath,'file')
    mn = mn+1;
    savepath = fullfile(datadir,sprintf('%s_%d_hires%d_ft%d_pt%d_%0.2fx_%d.avi', ...
    labels{vis},ix,highres,plotftracks,plotptracks,vframerateratio,mn));
    %     answer = input('File already exists. Overwrite? (y/n) ','s');
end

writerObj = VideoWriter( ...
    fullfile(savepath),'Motion JPEG AVI');
writerObj.Quality = 100;
writerObj.FrameRate = vframerate;
open(writerObj)

if plotptracks
    trackids = find(sum(tracks > 0,2) >= 5);
    trackids = trackids(~ismember(trackids,tracks2exc));
    ntracks = length(trackids);
    Ct = NaN(ntracks,tlen-1,2);
    for i = 1:ntracks
        for j = 1:tlen
            fridix = bst1+j-1;
            objid = tracks(trackids(i),j);
            if objid ~= 0
                Ct(i,j,:) = objects(fridix).centroid(objid,:);
            end
        end
    end
    
    if highres
        frames = find(ismember(aq.datainterp(:,1),im1cfrids));
        Ct_corr = cat(3,Ct(:,:,1)-(aq.datainterp(frames,3)-x_px2/2)', ...
            Ct(:,:,2)-(aq.datainterp(frames,2)-y_px2/2)');
        Ct_int = NaN(ntracks,len,2);
        for i = 1:ntracks
            lstix = find(im2cfrids == im1cfrids(end));
            Ct_int(i,:,:) = interp1(single(im1cfrids), ...
                squeeze(Ct_corr(i,:,:)),single(im2cfrids),'linear');
            Ct_int(i,lstix:end,:) = nan;
        end
        for i = 1:2
            Ct_int(:,:,i) = Ct_int(:,:,i)+offset(:,i)';
        end
        for i = 1:ntracks
            Ct_int(i,:,:) = smoothdata(squeeze(Ct_int(i,:,:)),1,'rloess',tfiltersz);
        end
    else
        Ct_int = Ct;
    end
end

% Show high resolution video
fig = figure;
for i = 1:len
    % Show image
    imagesc(im(:,:,i),imthr);
    hold on
    
    if highres
        colormap hot;
    else
        plot((aq.settings.FOV(1)/2)+xp,(aq.settings.FOV(2)/2)+yp, ...
            'Color',bcol,'LineWidth',2);
        %         colormap(flip(gray,1));
        colormap gray;
    end
    
    % Plot eye and tail features
    if plotftracks
        O = afap_im_segment_fish(im(:,:,i),eyethr,bodythr,1);
        
        lix = [O.stats.Side] == 0;
        rix = [O.stats.Side] == 1;
        lcp = O.stats(lix).Centroid;
        rcp = O.stats(rix).Centroid;
        lext = [O.stats(lix).ext_lower;
            O.stats(lix).ext_higher];
        rext = [O.stats(rix).ext_lower;
            O.stats(rix).ext_higher];
        
        tanlen = 10;
        ol = O.stats(lix).Orientation;
        if ol >= 0
            ol = ol+90;
        else
            ol = ol-90;
        end        
        or = O.stats(rix).Orientation;
        if or >= 0
            or = or-90;
        else
            or = or+90;
        end            
        
        lp2 = [lcp(1)+(cos(deg2rad(ol)).*tanlen),lcp(2)-(sin(deg2rad(ol)).*tanlen)];
        rp2 = [rcp(1)+(cos(deg2rad(or)).*tanlen),rcp(2)-(sin(deg2rad(or)).*tanlen)];
        
        scatter(O.tailpts(:,1),O.tailpts(:,2),'k','filled','MarkerEdgeColor',[1,1,1])
        plot(lcp(1),lcp(2),'b.')
        plot(rcp(1),rcp(2),'r.')
        %         line(lext(:,1),lext(:,2),'Color','b','LineWidth',2)
        %         line(rext(:,1),rext(:,2),'Color','r','LineWidth',2)
        line([lcp(1) lp2(1)],[lcp(2) lp2(2)],'Color','b','LineWidth',2)
        line([rcp(1) rp2(1)],[rcp(2) rp2(2)],'Color','r','LineWidth',2)
        
        eyes_ori = cat(1,O.stats(lix).Orientation,O.stats(rix).Orientation);
        eyes_malen = cat(1,O.stats(lix).MajorAxisLength,O.stats(rix).MajorAxisLength);
        eyes_cent = cat(1,O.stats(lix).Centroid,O.stats(rix).Centroid);
        
        x1 = eyes_cent(:,1)+(cos(deg2rad(eyes_ori)).*eyes_malen./2);
        x2 = eyes_cent(:,1)-(cos(deg2rad(eyes_ori)).*eyes_malen./2);
        y1 = eyes_cent(:,2)-(sin(deg2rad(eyes_ori)).*eyes_malen./2);
        y2 = eyes_cent(:,2)+(sin(deg2rad(eyes_ori)).*eyes_malen./2);
        
        [lidx,ridx] = visual_field_idx(8.5,O.cp,O.hvec, ...
            lcp,[x1(1),y1(1);x2(1),y2(1)], ...
            rcp,[x1(2),y1(2);x2(2),y2(2)], ...
            XG,YG);
        BL = bwboundaries(lidx);
        BR = bwboundaries(ridx);
        patch(BL{1}(:,2),BL{:}(:,1),[.5,.5,.5],'FaceAlpha',0.2,'EdgeColor',[.5,.5,.5]);
        patch(BR{1}(:,2),BR{:}(:,1),[.5,.5,.5],'FaceAlpha',0.2,'EdgeColor',[.5,.5,.5]);
    end
    
    % Plot paramecia tracks
    if plotptracks
        tid = aq.Hc([aq.Hc.hix] == ix).TargetID;
        tidix = trackids == tid;
        if usenames
            tstr = names(trackids);
        else
            tstr = num2str(trackids);
        end
        if plotptext
            if any(tidix)
                text(Ct_int(tidix,i,1)+5,Ct_int(tidix,i,2)-5,tstr(tidix,:),'Color','g')
            end
            text(Ct_int(~tidix,i,1)+5,Ct_int(~tidix,i,2)-5,tstr(~tidix,:),'Color',bcol)
        else
            if any(tidix)
                plot(Ct_int(tidix,i,1),Ct_int(tidix,i,2),'go','MarkerSize',mrksz)
            end
            plot(Ct_int(~tidix,i,1),Ct_int(~tidix,i,2),['o' bcol],'MarkerSize',mrksz)
        end
    end
    %             text(0.05,0.05,sprintf('%.0f/%.0f ms', ...
    %                 (i-toffset)/aq.settings.camerarate*1000, ...
    %                 (len-toffset)/aq.settings.camerarate*1000), ...
    %                 'Color',bcol,'Units','normalized');
    hold off
    axis square off
    drawnow
    
    thisaviframe = getframe(gca);
    writeVideo(writerObj,thisaviframe);
end
close(writerObj)
close(fig)
end
