function [] = afap_im_save_movie(datadir,fst,fed,motioncorr)
% Saves compressed movie from input frame interval to disk
%
% Pedro Henriques

chunksz = 10000;    % Chunk size in frames to process movie

[frid, x_px, y_px] = afaptdms_GetInfo(datadir,'im');
load(fullfile(datadir,'aq'),'aq');

writerObj = VideoWriter( ...
    fullfile(datadir,sprintf('im_movie_f%d_f%d',fst,fed')), ...
    'Motion JPEG AVI');
writerObj.FrameRate = 700;
writerObj.Quality = 100;
open(writerObj)

chks = fst:chunksz:fed;
if mod((fst+fed+1),chunksz) ~= 0
    chks = [chks, fed];
end
nchks = length(chks);

h = waitbar(0,'Init');
for i = 1:nchks-1
    waitbar(i/nchks-1,h,'processing')
    itst = chks(i);
    ited = chks(i+1);
    
    im = afap_im_extract(datadir,itst,ited,'im',frid,x_px,y_px,aq);
    
    if motioncorr
        [im,~] = afap_im_motion_correct(im,[105 90]);
    end
    
    im = permute(repmat(im,1,1,1,3),[1 2 4 3]);
    writeVideo(writerObj,im);
end
close(writerObj);
close(h);
end
