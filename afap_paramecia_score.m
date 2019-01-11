function [] = afap_paramecia_score(datadir,hixst)
%
% Allows user to selects the targets of all hunting events detected through
% the afap pipeline, using the high resolution videos
%
% The user is shown a clip of the hunting events and paramecia tracks
% detected through the paramecia tracking algorithm. At the end of the
% routine the user is prompted with the target id. Returning empty will
% result in the event being shown again. 0 should signify no target
% selected/visible
%
% Optionally the user may return 'j', to join two or more paramecia tracks
% (eg., [5241, 6895, 5584]), or 's', to change the speed to the video.
%
% Requires imbkd movie to have been recorded as well as having been
% processed through the afap_paramecia_seg and _tracks functions.
%
% Pedro Henriques, Sep 2017

%%  Hardcoded parameters

fstoffset = 200;
fendoffset = 700;
fovoffset = 50; % Field of view offset
imspeed = 5;
imthr = [0 50];

%%
if nargin < 1
    datadir = uigetdir;
end

% Load imbkd info
imstr = 'imbkd';
[frid, x_px, y_px] = afaptdms_GetInfo(datadir, imstr);
npx = x_px*y_px;

[frid2, x_px2, y_px2] = afaptdms_GetInfo(datadir, 'im');
npx2 = x_px2*y_px2;

% Load structures
disp('Loading structures');
load(fullfile(datadir,'aq'),'aq');
objfile = dir(fullfile(datadir,'objects*.mat'));
if ~isempty(objfile)
    objpath = fullfile(objfile(1).folder,objfile(1).name);
    load(objpath)
    t = matfile(fullfile(datadir,'tracks'));
    
    if nargin < 2
        % Check if to append to previous targets file
        if exist(fullfile(datadir,'Hc.mat'),'file')
            loadtarget = questdlg('Append to previous file?', ...
                'Please answer', ...
                'Yes','No','Yes');
            switch loadtarget
                case 'Yes'
                    load(fullfile(datadir,'Hc'),'Hc')
                    hixst = size(Hc,2)+1;
                otherwise
                    hixst = 1;
            end
        else
            hixst = 1;
        end
    end
    
    %%  Run the loop
    
    set(0,'DefaultFigureWindowStyle','docked')
    
    figure; colormap hot;
    for hix = hixst:size(aq.huntep,1)
        fprintf('Hunting event %d/%d\n',hix,size(aq.huntep,1));
        
        % Check it offsets pass index boudaries
        if aq.huntep(hix,1)-fstoffset <= 0
            fstoffset2 = 1;
        else
            fstoffset2 = fstoffset;
        end
        
        if aq.huntep(hix,2)+fendoffset > size(aq.datainterp,1)
            fendoffset2 = size(aq.datainterp,1);
        else
            fendoffset2 = fendoffset;
        end
        
        if aq.huntep(hix,1)-fstoffset2 <= 0
            continue
        end
        
        fst = aq.data(aq.huntep(hix,1)-fstoffset2,1);
        fend = aq.data(aq.huntep(hix,2)+fendoffset2,1);
        
        bst = afaptdms_GetFrLocation_ph(fst, frid);
        bend = afaptdms_GetFrLocation_ph(fend, frid);
        imcfrids = frid(bst:bend-1);
        len = bend-bst;
        
        % Load hires frames
        
        bst2 = afaptdms_GetFrLocation_ph(fst, frid2);
        bend2 = afaptdms_GetFrLocation_ph(fend, frid2);
        im2cfrids = frid2(bst2:bend2-1);
        len2 = bend2-bst2;
        
        im2 = TDMS_readTDMSFile(fullfile(datadir, ...
            ['im' '.tdms']), ...
            'SUBSET_GET', [bst2*npx2 bend2*npx2], ...
            'SUBSET_IS_LENGTH', false);
        
        im2 = reshape(im2.data{3}(1:npx2*len2)', [x_px2 y_px2 len2]);
        
        % load tracks for current image frames
        tracks = t.tracks(:,bst:bst+len);
        
        class = 0;
        while ~class
            FS = stoploop('Stop the loop? ');
            for i = 1:imspeed:len2
                if(~FS.Stop())    % Check if the loop has to be stopped
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
                    tracksobjix = tracksobjix(tracksobjix <= size(objects(bix).centroid,1));
                    
                    % Choose objects to draw
                    if ~isempty(trackids)
                        objcents = cat(1,objects(bix).centroid(tracksobjix,:));
                        objcents2draw = ...
                            objcents(:,1) >= fishcent(1)-(x_px2/2)-fovoffset & ...
                            objcents(:,1) <= fishcent(1)+(x_px2/2)+fovoffset & ...
                            objcents(:,2) >= fishcent(2)-(y_px2/2)-fovoffset & ...
                            objcents(:,2) <= fishcent(2)+(y_px2/2)+fovoffset;
                        if any(objcents2draw)
                            objcents_corr = [objcents(objcents2draw,1)+(x_px2/2)-fishcent(1), ...
                                objcents(objcents2draw,2)+(y_px2/2)-fishcent(2)];
                            plot(objcents_corr(:,1), ...
                                objcents_corr(:,2),'og')
                            text(objcents_corr(:,1)+5, ...
                                objcents_corr(:,2), ...
                                num2str(trackids(objcents2draw)), 'Color','g')
                        end
                    end
                    
                    % Color text in yellow if conv on
                    if frame >= aq.huntep(hix,1) && ...
                            frame <= aq.huntep(hix,2)
                        color = 'y';
                    elseif aq.vergF(frame) >= aq.gintersect
                        color = 'b';
                    else
                        color = 'w';
                    end
                    
                    text(10,170,sprintf('%.0f/%.0f ms', ...
                        (i-fstoffset)/aq.settings.camerarate*1000, ...
                        (len2-fstoffset)/aq.settings.camerarate*1000), ...
                        'Color',color);
                    hold off
                    axis([-fovoffset x_px2+fovoffset -fovoffset y_px2+fovoffset])
                    drawnow
                end
            end
            FS.Clear();
            clear FS;
            
            % Capture scoring
            objid = input('ObjID (options: ''j'' - join tracks; ''s'' - change video speed)\n','s');
            if objid == 'j'
                tracks2j = input('Which targets to join? ','s');
                afap_paramecia_tracks_join(datadir,str2num(tracks2j));
                tracks = t.tracks(:,bst:bst+len);
                continue
            elseif objid == 's'
                imspeed = input(sprintf('Change video speed from %d to ? ',imspeed),'s');
                imspeed = str2num(imspeed);
                continue
            end
            objid = str2num(objid);
            if isempty(objid)
                continue
            end
            
            if objid == 0
                code = 0;
                class = 1;
            else
                code = str2num(input('Code\n','s'));
                if ~isempty(code)
                    class = 1;
                else
                    continue
                end
            end
            
            notes = input('Notes\n','s');
        end
        
        Hc(hix).hix = hix;
        Hc(hix).TargetID = objid;
        Hc(hix).Code = code;
        Hc(hix).Notes = notes;
        save(fullfile(datadir,'Hc'),'Hc');
    end
else
    disp('No object file');
end


aq.Hc = Hc;
disp('Saving...')
save(fullfile(datadir,'aq'),'aq')
end