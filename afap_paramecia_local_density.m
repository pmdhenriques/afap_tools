function [] = afap_paramecia_local_density(datadir)
% Adds local paramecia density information to scored aq datasets
%
% Pedro Henriques

tau = 30;   % for exponential decay function
min_dist = 5;   % min distance (mm) for paramecia to be counted in local paramecia numbers
min_dist2 = 10;
min_ang = 120;  % min angle wr to fish as above
btlab = {'st','ed'};

h = waitbar(0,'Loading structures');

%% Load structures

load(fullfile(datadir,'aq'));

if isfield(aq.Hc,'Para_Local_n')
    aq.Hc = rmfield(aq.Hc,'Para_Local_n');
end
if isfield(aq.Hc,'Para_Local2_n')
    aq.Hc = rmfield(aq.Hc,'Para_Local2_n');
end

try
    load(fullfile(datadir,'objects'));
catch
    load(fullfile(datadir,'objects_v2.mat'));
end

tracksfile = matfile(fullfile(datadir,'tracks'));

[frid, x_px, y_px] = afaptdms_GetInfo(datadir,'im');
[frid2, ~, ~] = afaptdms_GetInfo(datadir,'imbkd');

ix = find(~cellfun(@isempty,{aq.Hc.frame}));    % valid indexes
for hix = 1:length(ix)
    waitbar(hix/length(ix),h,'Processing...')
    hixx = ix(hix);
    tid = aq.Hc(hixx).TargetID;
    tid = tid(1);
    if any(tid)
        tracks = tracksfile.tracks(tid,:);
        for f = 1:length(aq.Hc(hixx).frame)
            for i = 1:length(btlab)
                frame = aq.Hc(hixx).frame(f).(btlab{i});
                if ~isempty(frame)
                    buffer = aq.datainterp(frame,1);
                    fridix = findnearest(buffer,frid2,-1);
                    im2frame = find(aq.datainterp(:,1) == frid2(fridix));
                    objid = tracks(fridix);
                    
                    C_all = objects(fridix).centroid;   % All paramecia centroids
                    D = pdist2(C_all,C_all);    % Pairwise distances
                    Ld = nansum(exp(-D./50));   % Local density (exponential decay) for every paramecia
                    
                    %  Local paramecia numbers
                    
                    C_all_corr = [C_all(:,1)-(aq.datainterp(im2frame,3)-x_px/2), ...
                        C_all(:,2)-(aq.datainterp(im2frame,2)-y_px/2)]; % correct centroids to small fov
                    
                    im = afap_im_extract(datadir,im2frame,im2frame,'im',frid,x_px,y_px,aq);   % small fov image frame that corresponds to paramecia frame
                    O = afap_im_segment_fish(im,0.3,0.01);  % detect fish for this frame
                    
                    % Find missed paramecia with better segmentation
                    C_new = afap_paramecia_segment(im,2,'Both');
                    if ~isempty(C_new)
                        DD = pdist2(C_all_corr,C_new);
                        C_new_ix = ~any(DD < 5);
                        C_all_corr = [C_all_corr; C_new(C_new_ix,:)];
                    end
                    
                    if ~isempty(O)
                        cp = O.cp;  % Eye center point
                        hvec = O.hvec;  % Body vecror
                    else
                        cp = aq.Hc(hixx).cp_cent(f).(btlab{i});
                        hvec = aq.Hc(hixx).heading_vec(f).(btlab{i});
                    end
                    
                    fish2target_vec = [C_all_corr(:,1)-cp(1), ...
                        C_all_corr(:,2)-cp(2)];  % fish 2 paramecia vector
                    fish2target_ori = atan2d(hvec(1).*fish2target_vec(2)-hvec(2)*fish2target_vec(:,1), ...
                        hvec(1).*fish2target_vec(:,1) + hvec(2).*fish2target_vec(:,2)); % fish 2 paramecia orientation
                    fish2target_dist = sqrt((C_all_corr(:,1)-cp(1)).^2 + ...
                        (C_all_corr(:,2)-cp(2)).^2)./aq.settings.DISPLAYpxpermm;    % fish 2 paramecia distance (mm)
                    
                    % Number of paramecia in the near vicinity
                    Ln = sum(fish2target_dist <= min_dist & abs(fish2target_ori) <= min_ang);
                    Ln2 = sum(fish2target_dist <= min_dist2 & abs(fish2target_ori) <= min_ang);
                    
                    % Add to structure
                    
                    if objid ~= 0 && objid <= length(objects(fridix).centroid)
                        aq.Hc(hixx).Target_Ld(f).(btlab{i}) = Ld(objid);
                    else
                        aq.Hc(hixx).Target_Ld(f).(btlab{i}) = nan;
                    end
                    aq.Hc(hixx).All_Ld(f).(btlab{i}) = nansum(exp(-D./tau));
                    aq.Hc(hixx).Para_Local_n_5mm(f).(btlab{i}) = Ln;
                    aq.Hc(hixx).Para_Local_n_10mm(f).(btlab{i}) = Ln2;
                else
                    aq.Hc(hixx).Target_Ld(f).(btlab{i}) = nan;
                    aq.Hc(hixx).All_Ld(f).(btlab{i}) = nan;
                    aq.Hc(hixx).Para_Local_n_5mm(f).(btlab{i}) = nan;
                    aq.Hc(hixx).Para_Local_n_10mm(f).(btlab{i}) = nan;
                end
            end
        end
    end
end

waitbar(1,h,'Saving...')
save(fullfile(datadir,'aq'),'aq','-v7.3');
close(h)
end