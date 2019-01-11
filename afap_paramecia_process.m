function objects = afap_paramecia_process(datadir,fst,fed)
% Segments paramecia objects from the afap pipeline
%
% Uses old method, where all object files from NI ablation experiments were
% derived.
%
% Pedro Henriques

if nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\2P', ...
        'Select fish to process');
end

if ~exist(fullfile(datadir,'objects.mat'),'file')
    imstr = 'imbkd';
    [frid, x_px, y_px] = afaptdms_GetInfo(datadir, imstr);
    
    if nargin < 2
        fst = 1;
        if nargin < 3
            fed = length(frid)-1;
        end
    end
    
    thr = 0.01;
    ccthr = 80;
    immedfilter = [3 3];
    
    objects = struct(...
        'frame', [], ...
        'frid', [], ...
        'centroid', [], ...
        'area', [], ...
        'bbox', [], ...
        'imout', []);
    
    npx = x_px*y_px;
    
    %%
    parfor_progress(length(fst:fed));
    parfor frame = fst:fed
        parfor_progress;
        
        bst = npx*frame;
        bend = bst+npx-1;
        im = TDMS_readTDMSFile(fullfile(datadir, ...
            [imstr '.tdms']), ...
            'SUBSET_GET', [bst bend], ...
            'SUBSET_IS_LENGTH', false);
        im = reshape(im.data{3}', [x_px y_px]);
        
        %% Filtering and finding objects
        
        imF = medfilt2(im,immedfilter);
        imBW = imbinarize(imF,thr);
        imBW = bwmorph(bwmorph(bwmorph(imBW,'close',Inf),'fill'),'clean');
        cc = regionprops(imBW,'basic');
        cc = cc([cc.Area] <= ccthr);
        
        objects(frame).frame = frame;
        objects(frame).frid = frid(frame);
        objects(frame).bbox = cat(1,cc.BoundingBox);
        objects(frame).centroid = cat(1,cc.Centroid);
        objects(frame).area = cat(1,cc.Area);
    end
    parfor_progress(0);
    
    % Deal with empty frames
    ix = find(cellfun(@isempty,{objects.centroid}));
    for i = ix
        objects(i).centroid = [0 0];
    end
    
    disp('Saving...')
    save(fullfile(datadir,'objects.mat'),'objects');
else
    disp('objects file already exists!')
end
end