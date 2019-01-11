function [] = afap_paramecia_save_frames(datadir,imstr,fst,fed,frmskp,noutfr)
% Writes TIFF files of input tdms video from start (fst) and end (fed)
% frames, skiping 'frmskp' frames and outputting 'noutfr' number of frames
% for every frame skip.
% 
% Used to build small video frames across the whole experiment to manually
% count the number of paramecia.
%
% Pedro Henriques

if nargin < 6
    noutfr = 1;
    if nargin < 5
        frmskp = 1;
            if nargin < 3
                fst = 1;
                if nargin < 2
                    imstr = 'imbkd';
                    if nargin < 1
                        datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\NI project\Ablations', ...
                            'Select fish to process');
                    end; end; end; end; end

[frid, x_px, y_px] = afaptdms_GetInfo(datadir, imstr);
if nargin < 4
    fed = length(frid); end

h = waitbar(0,'Initializing');

npx = x_px*y_px;

%%

for frame = fst:frmskp:fed
    if frame+noutfr <= length(frid)
        fprogress = frame/length(frid);
        waitbar(fprogress,h, ...
            sprintf('Processing images... Frame (%d/%d)', ...
            frame,length(frid)));
        bst = npx*(frame);
        bend = (bst+npx*noutfr)-1;
        im = TDMS_readTDMSFile(fullfile(datadir, ...
            [imstr '.tdms']), ...
            'SUBSET_GET', [bst bend], ...
            'SUBSET_IS_LENGTH', false);
        im = reshape(im.data{3}', [x_px y_px noutfr]);
        
        % Save TIFF
        
        if ~exist(fullfile(datadir,'Paramecia_count'),'dir')
            mkdir(fullfile(datadir,'Paramecia_count'));
        end
        
        tf = Tiff(fullfile(datadir,'Paramecia_count', ...
            sprintf('%s_f%d_f%d.tif',imstr,frame,frame+noutfr-1)),'w');
        
        tags.Photometric = Tiff.Photometric.MinIsBlack;
        tags.Compression = Tiff.Compression.None;
        tags.ImageLength = y_px;
        tags.ImageWidth = x_px;
        tags.SampleFormat = Tiff.SampleFormat.UInt;
        tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tags.Orientation = Tiff.Orientation.TopLeft;
        tags.BitsPerSample = 8;
        
        for fr = 1:noutfr
            tf.setTag(tags)
            tf.write(im(:,:,fr));
            if ~(fr == noutfr)
                tf.writeDirectory()
            end
        end
        tf.close()
    end
end
close(h)
end