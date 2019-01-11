function imreg = afap_im_registration(imrot,imtemp)
% Does rigid registratin of afap im images to a template image.
% It is advised to run afap_im_rotate first
%
% Pedro Henriques, Oct 2017

[optimizer,metric] = imregconfig('monomodal');

if nargin < 2
    load('\\128.40.155.187\data2\Bianco_lab\Pedro\NI project\afap_analysis\imreg_mean', ...
        'imregmean');
    imtemp = imregmean;
end

imreg = zeros(size(imrot),'uint8');
for i = 1:size(imrot,3)
    imreg(:,:,i) = imregister(imrot(:,:,i),imtemp, ...
        'rigid',optimizer,metric,'PyramidLevels',4);
end
end