function Tpts = afap_im_segment_fish_tail(bodybw,cp,sbcp,npts,seglen)
% Finds tail point centroids from imput image from afap's 'im' tdms files.
% cp - Eye center point
% sbcp - Swim bladder centroid
% npts - number of points to find
% seglen - length of segment between each point
%
% Input parameters can be found by running the afap_im_segment_fish script.
% 
% DEPRICATED - afap_im_segment_fish now runs tail segmentation if prompt
%
% Pedro Henriques, 2018

bodydbw = bwdist(~bodybw);

Tpts = NaN(npts,2);
Tpts(1,:) = sbcp;  % First point is sb centroid
for i = 1:npts-1
    if i == 1
        cx = seglen*2*cos(0:0.1:2*pi)+Tpts(i,1);  % Draw a circle
        cy = seglen*2*sin(0:0.1:2*pi)+Tpts(i,2);
    else
        cx = seglen*cos(0:0.1:2*pi)+Tpts(i,1);
        cy = seglen*sin(0:0.1:2*pi)+Tpts(i,2);
    end
    [cxix,cyix,cI] = improfile(bodydbw,cx,cy,'bilinear');   % Build intencity profile of the circle perimeter
    cI2 = circshift(cI,round(length(cI)/4));    % shift intensity profile to correct for peaks that are in the extremities
    warning('off','all')
    [pks,locs] = findpeaks(cI,'MinPeakDistance',pi*seglen,'MinPeakHeight',1);  % Find peaks
    [pks2,~] = findpeaks(cI2,'MinPeakDistance',pi*seglen,'MinPeakHeight',1);   % Find peaks in the shifted vector
    warning('on','all')
    
    if length(pks) < 2 && length(pks2) >= 2     % Add a peak in the extremetry
        pks = pks2; locs = [1; locs];
    elseif length(pks) > 2
        % Use only two largest peaks
        [pks, ix] = sort(pks,'descend');
        pks = pks(1:2);
        locs = locs(ix(1:2));
    elseif length(pks) == 2
        if i == 1
            dd = abs(sqrt((cxix(locs)-cp(1)).^2 + (cyix(locs)-cp(2)).^2));  % Second point is defined by the maximum distance to the eye-eye center point
        else
            dd = NaN(length(pks),1);    % Next iterations look for points that maximize the distance to the previous point
            for j = 1:length(pks)
                dd(j) = abs(sqrt((cxix(locs(j))-Tpts(i-1,1))^2 + (cyix(locs(j))-Tpts(i-1,2))^2));
            end
        end
        [~, ix] = max(dd);
        Tpts(i+1,:) = [cxix(locs(ix)),cyix(locs(ix))];
    else
        % Terminate for loop if cant find more than one
        % point
        break
    end
end
end