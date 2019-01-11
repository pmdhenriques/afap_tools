function [C,varargout] = afap_paramecia_segment(im,pthr,method)
% Segments paramecia from afap movies and returns the centroid coordenates
%
% Pedro Henriques, Oct 2017

if nargin < 3
    method = 'Both';
end
if nargin < 2 || isempty(pthr)
    pthr = (max([min(max(im,[],1))  min(max(im,[],2))])); % Paramecia segmentation threshod        
end

gfilt = (fspecial('gaussian',7,1)); % Image filter
ndist = 5; % Minimum distance to find a neighbour between the different algorithms

switch method
    case 'LocalMaxima'
        % Use FastPeakFind to find local maxima
        [cent,~,imfilt] = FastPeakFind_ph(im,pthr,gfilt,2,1);
    case 'WeightedCentroid'
        % Slower method to find more difficult objects
        cent = FastPeakFind_ph(im,pthr,gfilt,2,2);
    case 'Both'
        [cent,~,imfilt] = FastPeakFind_ph(im,pthr,gfilt,2,1);
        cent = [cent; FastPeakFind_ph(im,pthr,gfilt,2,2)];
end

% Reshape
C = zeros(length(cent)/2,2);
C(:,1) = cent(1:2:end);
C(:,2) = cent(2:2:end);

switch method
    case {'LocalMaxima','Both'}
        % Delete centroids that colocalize with the fish body
        imbw = bwmorph(imfilt > 1,'open');        
        s = regionprops(imbw,{'Area','PixelList'});        
        [~,ix] = max([s.Area]);
        sf = s(ix);       
        sp = s(setdiff(1:size(s,1),ix));
        idx = ismember(C,sf.PixelList,'rows');
        C(idx,:) = [];
        
        
        % Compute pairwise distances
        
        dd = triu(pdist2(C,C));
        
        % Find centroids with no close neighbours and add them to final centroids
        
        [r,~] = find(dd > 0 & dd < ndist);
        if ~isempty(r)
            C(r,:) = [];
        end
end

if exist('imfilt','var')
    varargout{1} = imfilt;
end
varargout{2} = sp;
end