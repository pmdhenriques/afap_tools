function [] = afap_paramecia_score_info(datadir,hispeedtrack)
% Computes several parameters of fish 2 paramecia interaction from scored
% aq structures (requires Hc field).
%
% Pedro Henriques, Dez 2017

if nargin < 2
    hispeedtrack = 1;
    if nargin < 1
        datadir = uigetdir;
    end
end

overwrite = 'Yes';
overbypass = 1; % Bypass overwrite question (for batch processing purposes)
showfig = 0;

% Hardcoded parameters
nobjfrms = 10;  % Number of frames before start of bout to extract centroids
eyethr = 0.3;   % eye thresolh for fish segmentation
bodythr = 0.01; % body threshold for fish segmentation

% Load image metadata
thisimj = 'im';
thisimj2 = 'imbkd';
[frid, x_px, y_px] = afaptdms_GetInfo(datadir, thisimj);
[frid2, x_px2, y_px2] = afaptdms_GetInfo(datadir, thisimj2);
npx = x_px*y_px;
npx2 = x_px2*y_px2;

%%  Load structures

h = waitbar(0,'Initializing');
waitbar(0,h,'Loading structures');

load(fullfile(datadir,'aq'),'aq');
if isfield(aq,'Hc')
    fldnames = fieldnames(aq.Hc);
    if length(fldnames) > 4
        if ~overbypass
            overwrite = questdlg('Overwrite?','Overwrite','Yes','No','Yes');
        end
        switch overwrite
            case 'Yes'
                aq.Hc = rmfield(aq.Hc,fldnames(5:end));
                disp('Deleting previous Hc fields')
        end
    end
else
    fprintf('No Hc structure in %s\n',datadir);
    return
end

objfile = dir(fullfile(datadir,'objects*.mat'));
if ~isempty(objfile)
    objpath = fullfile(objfile(1).folder,objfile(1).name);
    load(objpath)
else
    fprintf('No objects file found on %s\n',datadir);
    return
end

if exist(fullfile(datadir,'tracks.mat'),'file')
    tracksfile = matfile(fullfile(datadir,'tracks'));
    ntracks = size(tracksfile,'tracks',1);
else
    fprintf('No tracks file found on %s\n',datadir);
    return
end

%%  Main loop

% hunting indices
hixs = find(cellfun(@any,{aq.Hc.TargetID}));  % instances where target ID was identified

k = 1;
for hixx = hixs
    waitbar(k/length(hixs),h,'Processing...')
    tid = aq.Hc(hixx).TargetID; % target id
    if any(tid)
        
        tid = tid(1);
        if tid > 0 && tid <= ntracks && tid == round(tid)
            hix = aq.Hc(hixx).hix;
            tracks = tracksfile.tracks(tid,:);   % load track
            
            btsix = find([aq.btf.Convbout] == hix); % Bouts contained in hunting event
            
            if isempty(btsix)
                % No bout within hunting event
                continue
            end
            
            framest = aq.btf(btsix(1)).itst;
            frameend = aq.btf(btsix(end)).ited;
            
            bufferst = aq.datainterp(framest,1);
            bufferend = aq.datainterp(frameend,1);
            
            % im
            fridixst = afaptdms_GetFrLocation_ph(bufferst, frid);
            fridixend = afaptdms_GetFrLocation_ph(bufferend, frid);
            
            imbuffers = frid(fridixst:fridixend);
            len = length(imbuffers);
            imframes = find(ismember(aq.datainterp(:,1),imbuffers));
            
            % imbkd
            fridixst2 = findnearest(frid2,bufferst,1); fridixst2 = fridixst2(1);
            if frid2(fridixst2) > bufferst && fridixst2 ~= 1
                fridixst2 = fridixst2-1; end
            
            fridixend2 = findnearest(frid2,bufferend,-1); fridixend2 = fridixend2(1);
            if frid2(fridixend2) < bufferend
                fridixend2 = fridixend2+1; end
            
            im2buffers = frid2(fridixst2:fridixend2);
            len2 = length(im2buffers);
            im2frames = find(ismember(aq.datainterp(:,1),im2buffers));
            
            if fridixst2-nobjfrms > 0
                
                % Load hi res image
                
%                 try
                    pxst = npx*(fridixst-1)+1;
                    pxend = npx*fridixend;
                    im = TDMS_readTDMSFile(fullfile(datadir, ...
                        [thisimj '.tdms']), ...
                        'SUBSET_GET', [pxst pxend], ...
                        'SUBSET_IS_LENGTH', false);
                    im = reshape(im.data{3}(1:npx*len)', [x_px y_px len]);
%                 catch
%                     k = k+1;
%                     continue
%                 end
                
                % Low res image (not used; for debug)
                
%                 pxst2 = npx2*(fridixst2-1)+1;
%                 pxend2 = npx2*fridixend2;
%                 im2 = TDMS_readTDMSFile(fullfile(datadir, ...
%                     [thisimj2 '.tdms']), ...
%                     'SUBSET_GET', [pxst2 pxend2], ...
%                     'SUBSET_IS_LENGTH', false);
%                 im2 = reshape(im2.data{3}(1:npx2*len2)', [x_px2 y_px2 len2]);
                
                %%  Obtain paramecia centroids
                
                C_low_wide = NaN(len2,2);
                c = 1;
                for fridix = fridixst2:fridixend2
                    objid = tracks(fridix);
                    if objid == 0
                        C_low_wide(c,:) = [nan nan];
                    else
                        if objid <= length(objects(fridix).centroid)
                            C_low_wide(c,:) = objects(fridix).centroid(objid,:);
                        else
                            fprintf('Error in obj id @ frame %d\n',fridix)
                            C_low_wide(c,:) = [nan nan];
                        end
                    end
                    c = c+1;
                end
                
                % Correct centroids to small fov
                C_low_small = [C_low_wide(:,1)-(aq.datainterp(im2frames,3)-x_px/2), ...
                    C_low_wide(:,2)-(aq.datainterp(im2frames,2)-y_px/2)];
                
                mix = isnan(C_low_wide(:,1));
                
                % Requires at least 4 datapoints to interpolate
                if sum(~mix) >= 4
                    
                    % Interpolate any missing values
                    if any(mix)
                        C_low_wide_int = interp1(find(~mix), ...
                            C_low_wide(~mix,:),1:size(C_low_wide,1),'linear');
                    else
                        C_low_wide_int = C_low_wide;
                    end
                    
                    % Interpolate to hi res timebase
                    C_hi_wide_int = interp1(single(im2buffers), ...
                        C_low_wide_int,single(imbuffers), ...
                        'linear');
                    lastnans_hi = find(~isnan(C_hi_wide_int(:,1)));
                    lastnans_ix = false(len,1);
                    
                    if length(lastnans_hi) >= 4
                        if lastnans_hi(end) < len
                            lastnans_ix(lastnans_hi(end)+1:end) = 1;
                        end
                        
                        %%  Replace frames with tracking in hi resolution
                        
                        % Correct centroids to small fov
                        C_hi_small_int = [C_hi_wide_int(:,1)-(aq.datainterp(imframes,3)-x_px/2), ...
                            C_hi_wide_int(:,2)-(aq.datainterp(imframes,2)-y_px/2)];
                        
                        %                     % Create matrix to be filled out with updated positions
                        %                     C_hi_small2 = NaN(size(C_hi_small_int));
                        %                     % Fill with already known positions
                        %                     [Lia,Locb] = ismember(imbuffers,im2buffers);
                        %                     C_hi_small2(Lia,:) = C_low_small(Locb(Locb ~= 0),:);
                        
                        C_hi_small2 = C_hi_small_int;
                        
                        % Check which centroids are within small fov boudaries
                        C_hi_small_ix = all(C_hi_small_int > 0 & ...
                            C_hi_small_int < x_px | ...
                            lastnans_ix, ...
                            2);
                        
                        if all(C_hi_small_ix)
                            st = 1;
                            ed = length(C_hi_small_ix);
                        elseif all(~C_hi_small_ix)
                            st = []; ed = [];
                        else
                            st = find(diff(C_hi_small_ix) == 1);
                            ed = find(diff(C_hi_small_ix) == -1);
                            if isempty(st) && ~isempty(ed)
                                st = [1;st];
                            elseif isempty(ed) && ~isempty(st)
                                ed = [ed; length(C_hi_small_ix)];
                            end
                            if ed(1) < st(1)
                                st = [1;st];
                            end
                            if st(end) > ed(end)
                                ed = [ed; length(C_hi_small_ix)];
                            end
                            %                         rmvix = ed-st < 70;
                            %                         st(rmvix) = []; ed(rmvix) = [];
                        end
                        
                        if ~isempty(st) && hispeedtrack
                            for m = 1:length(st)
                                cframes = st(m):ed(m);
                                nframes = length(cframes);
                                objects_hi = cell(nframes,1);
                                for f1 = 1:nframes
                                    frame = cframes(f1);
                                    objects_hi{f1} = afap_paramecia_segment(im(:,:,frame),2,'Both');
                                    % new method
                                    if ~lastnans_ix(frame)
                                        DD = pdist2(C_hi_small_int(frame,:),single(objects_hi{f1}));
                                        DDv = DD < 10;
                                        if any(DDv)
                                            [~,DDvix] = min(DD);
                                            C_hi_small2(frame,:) = objects_hi{f1}(DDvix,:);
                                        end
                                    else
                                        lastvalue_ix = find(~isnan(C_hi_small2(:,1)));
                                        lastvalue_ix = lastvalue_ix(end);
                                        if ~isempty(lastvalue_ix)
                                            if ~isempty(objects_hi{f1})
                                                for i = 0:9
                                                    thislastvalue = single(C_hi_small2(lastvalue_ix-i,:));
                                                    DD = pdist2(thislastvalue,single(objects_hi{f1}));
                                                    DDv = DD < 10;
                                                    if any(DDv)
                                                        [~,DDvix] = min(DD);
                                                        C_hi_small2(frame,:) = objects_hi{f1}(DDvix,:);
                                                        break
                                                    end
                                                end
                                            end
                                        else
                                            break
                                        end
                                    end
                                end
                            end
                        end
                        
                        % interpolate new trajectory
                        mix2 = isnan(C_hi_small2(:,1));
                        C_hi_small2_int = interp1(single(imbuffers(~mix2)), ...
                            C_hi_small2(~mix2,:),single(imbuffers), ...
                            'linear');
                        
                        firstvalue = find(~isnan(C_hi_small2(:,1)));
                        firstvalue = firstvalue(1);
                        if firstvalue ~= 1
                            C_hi_small2_int(1:firstvalue-1,:) = C_hi_small_int(1:firstvalue-1,:);
                        end
                        
                        %%  Compute measurements on bout start and end
                        
                        bt = 1;
                        for btix = btsix
                            btstit = aq.btf(btix).itst;
                            btendit = aq.btf(btix).ited;
                            
                            % image frame
                            imixs = [btstit-framest+1 ...
                                btendit-framest];
                            
                            % Bout start and end
                            for b = 1:2
                                % Segment fish
                                O = afap_im_segment_fish(im(:,:,imixs(b)),eyethr,bodythr);
                                
                                if ~isempty(O)
                                    %%
                                    if b == 1
                                        % Mean paramecia trajectory preceding bout
                                        buffer_low = aq.datainterp(btstit,1);
                                        frame_low = afaptdms_GetFrLocation_ph(buffer_low, frid2);
                                        C_prev = NaN(nobjfrms,2);
                                        c = 1;
                                        for fridix = frame_low-nobjfrms+1:fridixst2
                                            objid = tracks(fridix);
                                            if objid == 0
                                                C_prev(c,:) = [nan nan];
                                            else
                                                if objid <= length(objects(fridix).centroid)
                                                    C_prev(c,:) = objects(fridix).centroid(objid,:);
                                                else
                                                    C_prev(c,:) = [nan nan];
                                                end
                                            end
                                            c = c+1;
                                        end
                                        
                                        o_mvec = [nanmean(diff(C_prev)) 0];
                                        o_mori = atan2d(o_mvec(2),o_mvec(1));
                                        
                                        % Compute parameters
                                        aq.Hc(hixx).frame(bt).st = btstit;  % frame
                                        target_cent_hi_small = C_hi_small2_int(imixs(b),:);
                                        cp_cent_hi_wide = [(aq.datainterp(btstit,3)-x_px/2)+O.cp(1), ...
                                            (aq.datainterp(btstit,2)-y_px/2)+O.cp(2)];  % eye center point corrected for low res FOV
                                        fish2target_vec = [target_cent_hi_small(1)-O.cp(1), ...
                                            target_cent_hi_small(2)-O.cp(2), 0];  % fish 2 paramecia vector
                                        fish2target_ori = atan2d(O.hvec(1)*fish2target_vec(2)-O.hvec(2)*fish2target_vec(1), ...
                                            O.hvec(1)*fish2target_vec(1)+O.hvec(2)*fish2target_vec(2)); % fish 2 paramecia orientation
                                        fish2target_traj_ori = atan2d(O.hvec(1)*o_mvec(2)-O.hvec(2)*o_mvec(1), ...
                                            O.hvec(1)*o_mvec(1)+O.hvec(2)*o_mvec(2));   % 180 deg == oposite directions; 0 == same directions
                                        fish2target_dist = sqrt((target_cent_hi_small(1)-O.cp(1))^2 + ...
                                            (target_cent_hi_small(2)-O.cp(2))^2);    % fish 2 paramecia distance (px)
                                        
                                        % Write to struct
                                        aq.Hc(hixx).target_cent_hi_small(bt).st = target_cent_hi_small;
                                        aq.Hc(hixx).cp_cent(bt).st = O.cp;
                                        aq.Hc(hixx).cp_cent_corr(bt).st = cp_cent_hi_wide;
                                        aq.Hc(hixx).fish2target_vec(bt).st = fish2target_vec;
                                        aq.Hc(hixx).fish2target_ori(bt).st = fish2target_ori;
                                        aq.Hc(hixx).fish2target_traj_ori(bt).st = fish2target_traj_ori;
                                        aq.Hc(hixx).fish2target_dist(bt).st = fish2target_dist;
                                        aq.Hc(hixx).target_mean_vec(bt).st = o_mvec;
                                        aq.Hc(hixx).target_mean_ori(bt).st = o_mori;
                                        aq.Hc(hixx).heading_vec(bt).st = O.hvec;
                                        aq.Hc(hixx).heading_ori(bt).st = O.h_ori;
                                        
                                        %                                 if length(btsix) >= 2 && btix == btsix(end-1)
                                        %                                     figure
                                        %                                     imagesc(im(:,:,imixs(b)))
                                        %                                     hold on
                                        %                                     plot(target_cent_hi_small(1),target_cent_hi_small(2),'ko')
                                        %                                     plot(O.cp(1),O.cp(2),'ro')
                                        %                                     plot(C_hi_small_int(imixs(b),1),C_hi_small_int(imixs(b),2),'r*')
                                        %                                     text(140,20,sprintf('f2td = %0.2g',fish2target_dist),'color','w')
                                        %                                     text(140,30,sprintf('f2to = %0.2g',fish2target_ori),'color','w')
                                        %                                     title(sprintf('hixx = %d (%d)',hixx,showfig))
                                        %                                     hold off
                                        %                                     showfig = 0;
                                        %                                 end
                                    else
                                        % Compute parameters
                                        aq.Hc(hixx).frame(bt).ed = btendit;  % frame
                                        target_cent_hi_small = C_hi_small2_int(imixs(b),:);
                                        cp_cent_hi_wide = [(aq.datainterp(btendit,3)-x_px/2)+O.cp(1), ...
                                            (aq.datainterp(btendit,2)-y_px/2)+O.cp(2)];  % eye center point corrected for low res FOV
                                        fish2target_vec = [target_cent_hi_small(1)-O.cp(1), ...
                                            target_cent_hi_small(2)-O.cp(2), 0];  % fish 2 paramecia vector
                                        fish2target_ori = atan2d(O.hvec(1)*fish2target_vec(2)-O.hvec(2)*fish2target_vec(1), ...
                                            O.hvec(1)*fish2target_vec(1)+O.hvec(2)*fish2target_vec(2)); % fish 2 paramecia orientation
                                        fish2target_dist = sqrt((target_cent_hi_small(1)-O.cp(1))^2 + ...
                                            (target_cent_hi_small(2)-O.cp(2))^2);
                                        
                                        % Write to struct
                                        aq.Hc(hixx).target_cent_hi_small(bt).ed = target_cent_hi_small;
                                        aq.Hc(hixx).cp_cent(bt).ed = O.cp;
                                        aq.Hc(hixx).cp_cent_corr(bt).ed = cp_cent_hi_wide;
                                        aq.Hc(hixx).fish2target_vec(bt).ed = fish2target_vec;
                                        aq.Hc(hixx).fish2target_ori(bt).ed = fish2target_ori;
                                        aq.Hc(hixx).fish2target_dist(bt).ed = fish2target_dist;
                                        aq.Hc(hixx).heading_vec(bt).ed = O.hvec;
                                        aq.Hc(hixx).heading_ori(bt).ed = O.h_ori;
                                        
                                    end
                                end
                            end
                            bt = bt+1;
                        end
                    end
                end
            end
        end
    end
    k = k+1;
end

%%
waitbar(1,h,'Saving...')
save(fullfile(datadir,'aq'),'aq','-v7.3');

close(h)
end