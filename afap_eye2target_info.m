function afap_eye2target_info(datadir)
% Gets target parameters (orientation) with respect to eye position.
% Needs scored movies
%
% Pedro Henriques, November 2017

if nargin < 1
    datadir = uigetdir;
end

thisimj = 'im';
[frid, x_px, y_px] = afaptdms_GetInfo(datadir, thisimj);
npx = x_px*y_px;

load(fullfile(datadir,'aq'),'aq');

eyethr = 0.3;   % eye thresolh for fish segmentation
bodythr = 0.01; % body threshold for fish segmentation

% hixs = find(cellfun(@(x) x > 0,{aq.Hc.TargetID}));  % instances where target ID was identified
hixs = find(cellfun(@any,{aq.Hc.TargetID}));

k = 1;
h = waitbar(0,'Initializing');
for hixx = hixs
    waitbar(k/length(hixs),h,'Processing')
    for b = 1:size(aq.Hc(hixx).frame,2)
        
        framest = aq.Hc(hixx).frame(b).st;
        frameend = aq.Hc(hixx).frame(b).ed;
        
        if ~isempty(framest) && ~isempty(frameend)
            bufferst = aq.datainterp(framest,1);
            bufferend = aq.datainterp(frameend,1);
            
            % im
            fridixst = afaptdms_GetFrLocation_ph(bufferst, frid);
            fridixend = afaptdms_GetFrLocation_ph(bufferend, frid);
            len = fridixend-fridixst;
            
            pxst = npx*(fridixst-1)+1;
            pxend = npx*(fridixend-1)+1;
            im = TDMS_readTDMSFile(fullfile(datadir, ...
                [thisimj '.tdms']), ...
                'SUBSET_GET', [pxst pxend], ...
                'SUBSET_IS_LENGTH', false);
            im = reshape(im.data{3}(1:npx*len)', [x_px y_px len]);
            
            
            for i = 1:2
                if i == 1
                    frame = framest;
                    O = afap_im_segment_fish(im(:,:,1),eyethr,bodythr);
                else
                    frame = frameend;
                    O = afap_im_segment_fish(im(:,:,end),eyethr,bodythr);
                end
                
                lix = [O.stats.Side] == 0;
                rix = [O.stats.Side] == 1;
                lc = O.stats(lix).Centroid;
                rc = O.stats(rix).Centroid;
                
                tc = aq.Hc(hixx).target_cent(b).st;
                
                fc = aq.datainterp(frame,2:3);
                
                lc_lfov = [(fc(2)-x_px/2)+lc(1), ...
                    (fc(1)-y_px/2)+lc(2)];
                rc_lfov = [(fc(2)-x_px/2)+rc(1), ...
                    (fc(1)-y_px/2)+rc(2)];
                
                l2tv = [tc(1)-lc_lfov(1), tc(2)-lc_lfov(2), 0];
                r2tv = [tc(1)-rc_lfov(1), tc(2)-rc_lfov(2), 0];
                
                lev = O.stats(lix).eye_vector.* [1 -1 1];
                rev = O.stats(rix).eye_vector.* [1 -1 1];
                
                Leye2target_ori = atan2d(lev(1)*l2tv(2)-lev(2)*l2tv(1), ...
                    lev(1)*l2tv(1)+lev(2)*l2tv(2));
                
                Reye2target_ori = atan2d(rev(1)*r2tv(2)-rev(2)*r2tv(1), ...
                    rev(1)*r2tv(1)+rev(2)*r2tv(2));
                
                if i == 1
                    aq.Hc(hixx).Leye2target_ori(b).st = Leye2target_ori;
                    aq.Hc(hixx).Reye2target_ori(b).st = Reye2target_ori;
                else
                    aq.Hc(hixx).Leye2target_ori(b).ed = Leye2target_ori;
                    aq.Hc(hixx).Reye2target_ori(b).ed = Reye2target_ori;
                end
            end
        end
    end
    k = k+1;
end

waitbar(1,h,'Saving')
save(fullfile(datadir,'aq'),'aq');
close(h)
end