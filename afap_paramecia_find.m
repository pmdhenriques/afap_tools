function objects = afap_paramecia_find(datadir,frames)
% Segments paramecia objects from the afap pipeline
%
% Uses new method from afap_paramecia_segment. Not tested fully on whole
% FOV datasets.
%
% Pedro Henriques, Oct 2017

if nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P', ...
        'Select fish to process');
end

imstr = 'imbkd';
[frid, x_px, y_px] = afaptdms_GetInfo(datadir, imstr);

if nargin < 2
    frames = 1:length(frid)-1;
end

objects = struct(...
    'frame', [], ...
    'frid', [], ...
    'centroid', [], ...
    'imout', []);

npx = x_px*y_px;

%%
parfor_progress(length(frames));
parfor frame = frames
    parfor_progress;
    
    bst = npx*frame;
    bend = bst+npx-1;
    im = TDMS_readTDMSFile(fullfile(datadir, ...
        [imstr '.tdms']), ...
        'SUBSET_GET', [bst bend], ...
        'SUBSET_IS_LENGTH', false);
    im = reshape(im.data{3}', [x_px y_px]);
    
    %% finding paramecia
    
    C = afap_paramecia_segment(im);
    
    objects(frame).frame = frame;
    objects(frame).frid = frid(frame);
    objects(frame).centroid = C;
end
parfor_progress(0);

disp('Saving...')
save(fullfile(datadir,'objects_v2.mat'),'objects');
       
end