function [frames] = afap_im_get_frames(datadir,vis,lat,ix,aq)
% Get frame start and end for input index of visual stimulus
%
% Example:
% vis -- 1 = loom; 2 = tap; 3 = grat; 4 = hunt
% lat -- 0 = Left; 1 = Right
% ix -- 24

% Load structures
if nargin < 5
    disp('Loading aq');
    load(fullfile(datadir,'aq'),'aq');
end

if vis == 1
    fst = aq.loomz(ix,1);
    fed = aq.loomz(ix,2);
elseif vis == 2
    fst = aq.tapz(ix,1);
    fed = aq.tapz(ix,2);
elseif vis == 3
    fixs = find(aq.visstim(:,1) == vis & aq.visstim(:,2) == lat);
    fst = aq.visstim(fixs(ix),3);
    fed = aq.visstim(fixs(ix),4);
elseif vis == 4
    fst = aq.huntep(ix,1);
    fed = aq.huntep(ix,2);
end

frames = [fst fed];

end