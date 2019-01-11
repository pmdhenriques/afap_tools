function [im] = afap_im_extract(datadir,itst,ited,imstr,frid,x_px,y_px,aq)
% Extracts image frames from the TDMS files for the iterations inputed.
%
% Pedro Henriques, Oct 2017

proceed = 'y';

if nargin < 8
    [frid, x_px, y_px] = afaptdms_GetInfo(datadir, imstr);
    load(fullfile(datadir,'aq'),'aq');
elseif nargin < 4
    imstr = 'im';
end

npx = x_px*y_px;

bufferst = aq.datainterp(itst,1);
bufferend = aq.datainterp(ited,1);

fridixst = findnearest(frid,bufferst,1);
fridixend = findnearest(frid,bufferend,1);

len = length(fridixst:fridixend);

pxst = npx*(fridixst-1)+1;
pxend = npx*fridixend;

if len > 15000
    proceed = input('Length > 15k frames. Proceed? [y/n] ','s');
end

switch proceed
    case 'y'
        im = TDMS_readTDMSFile(fullfile(datadir, ...
            [imstr '.tdms']), ...
            'SUBSET_GET', [pxst pxend], ...
            'SUBSET_IS_LENGTH', false);
        
        im = reshape(im.data{3}(1:npx*len)', [x_px y_px len]);
    otherwise
        return
end
end