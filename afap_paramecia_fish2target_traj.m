function [tori,tdist,C_hi_small2_int] = afap_paramecia_fish2target_traj(datadir,tid,fst,fed,frid,x_px,y_px,frid2,aq)
% Reports fish-target id features for given targetID (tid) at input frame
% interval

h = waitbar(0,'Initializing');
eyethr = 0.3;   % eye thresolh for fish segmentation
bodythr = 0.01; % body threshold for fish segmentation
thisimj = 'im';
thisimj2 = 'imbkd';
    
if nargin < 5
    % Load image metadata
    [frid, x_px, y_px] = afaptdms_GetInfo(datadir, thisimj);
    [frid2, ~, ~] = afaptdms_GetInfo(datadir, thisimj2);    
    % Load aq
    waitbar(0,h,'Loading aq structure...')
    load(fullfile(datadir,'aq'),'aq');    
end
npx = x_px*y_px;

% Load objects file
objfile = dir(fullfile(datadir,'objects*.mat'));
if ~isempty(objfile)
    objpath = fullfile(objfile(1).folder,objfile(1).name);
    waitbar(0,h,'Loading objects file...')
    load(objpath)
else
    fprintf('No objects file found on %s\n',datadir);
    return
end

% Load object tracks file
if exist(fullfile(datadir,'tracks.mat'),'file')
    tracksfile = matfile(fullfile(datadir,'tracks'));
else
    fprintf('No tracks file found on %s\n',datadir);
    return
end

% load track
tracks = tracksfile.tracks(tid,:);   

% image buffers
bst = aq.datainterp(fst,1);
bed = aq.datainterp(fed,1);

% im
fridixst = afaptdms_GetFrLocation_ph(bst, frid);
fridixend = afaptdms_GetFrLocation_ph(bed, frid);

imbuffers = frid(fridixst:fridixend);
len = length(imbuffers);
imframes = find(ismember(aq.datainterp(:,1),imbuffers));

% imbkd
fridixst2 = findnearest(frid2,bst,1); fridixst2 = fridixst2(1);
if frid2(fridixst2) > bst && fridixst2 ~= 1
    fridixst2 = fridixst2-1; end

fridixend2 = findnearest(frid2,bed,-1); fridixend2 = fridixend2(1);
if frid2(fridixend2) < bed
    fridixend2 = fridixend2+1; end

im2buffers = frid2(fridixst2:fridixend2);
len2 = length(im2buffers);

% Load hi res image
waitbar(0,h,'Loading image frames...')
pxst = npx*(fridixst-1)+1;
pxend = npx*fridixend;
im = TDMS_readTDMSFile(fullfile(datadir, ...
    [thisimj '.tdms']), ...
    'SUBSET_GET', [pxst pxend], ...
    'SUBSET_IS_LENGTH', false);
im = reshape(im.data{3}(1:npx*len)', [x_px y_px len]);

% Get low frame rate centroids for track
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
        end
        
        if ~isempty(st)
            for m = 1:length(st)
                waitbar(m/length(st),h,'Improving centroid accuracy...')
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
        
        tori = NaN(len,1);
        tdist = NaN(len,1);
        for f = 1:len
            waitbar(f/len,h,'Computing distances...')
            O = afap_im_segment_fish(im(:,:,f),eyethr,bodythr);
            
            if ~isempty(O)
                %%
                % Compute parameters
                target_cent_hi_small = C_hi_small2_int(f,:);
                fish2target_vec = [target_cent_hi_small(1)-O.cp(1), ...
                    target_cent_hi_small(2)-O.cp(2), 0];  % fish 2 paramecia vector
                fish2target_ori = atan2d(O.hvec(1)*fish2target_vec(2)-O.hvec(2)*fish2target_vec(1), ...
                    O.hvec(1)*fish2target_vec(1)+O.hvec(2)*fish2target_vec(2)); % fish 2 paramecia orientation
                fish2target_dist = sqrt((target_cent_hi_small(1)-O.cp(1))^2 + ...
                    (target_cent_hi_small(2)-O.cp(2))^2);    % fish 2 paramecia distance (px)
                
                % Write to struct
                tori(f,1) = fish2target_ori;
                tdist(f,1) = fish2target_dist;
            end
        end
    else
        error('Not enough data points...')
    end
else
    error('Not enough data points...')
end
close(h);
end