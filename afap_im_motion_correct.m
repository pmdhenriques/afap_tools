function [imt,offset] = afap_im_motion_correct(im,refpos)
% Translated input image sequence to either a reference position or the
% first weighted centroid of the input sequence
%
% Pedro Henriques, Fev 2018

thr = 0.1;  % Body threshold

if nargin < 2
    refpos = [90 90];
end

nfr = size(im,3);
offset = zeros(nfr,2);

imt = zeros(size(im),'uint8');
imt(:,:,1) = im(:,:,1);

h = waitbar(0,'Init');
for i = 1:size(im,3)
    waitbar(i/size(im,3),h,'motion correcting')
    imbw = imbinarize(im(:,:,i),thr);
    stats = regionprops(imbw,im(:,:,i),{'WeightedCentroid','Centroid','Area'});
    [~,mix] = max([stats.Area]);
    ccent = stats(mix).WeightedCentroid;
    if i == 1 && isempty(refpos)
        refpos = ccent;
        continue
    end
    offset(i,:) = refpos-ccent;
    imt(:,:,i) = imtranslate(im(:,:,i), ...
        offset(i,:),'cubic');
end
close(h)
end
