function afap_paramecia_tracks_reparse(datadir,n_tracks2p)
% Computes the tracks matrix when the input data from simpletracker_ph
% doesn't fit into memory.

if nargin < 2
    n_tracks2p = 50000;
end

disp('Loading structures')

objfile = dir(fullfile(datadir,'objects*.mat'));
if ~isempty(objfile)
    objpath = fullfile(objfile(1).folder,objfile(1).name);
    load(objpath)
else
    disp('No objects file found!')
    return
end

load(fullfile(datadir,'tracks'))
adjacency_tracks = tracks;
clear tracks

points = {objects.centroid};
maxpoints = max(cellfun(@length,points));
n_cells = cellfun(@(x) size(x, 1), points);

if maxpoints <= 255
    datatype = 'uint8';
elseif maxpoints > 255 && maxpoints <= 65535
    datatype = 'uint16';
else
    datatype = 'uint32';
end

n_slices = numel(points);
n_tracks = length(adjacency_tracks);

%%

tfilename = fullfile(datadir,'tracks_rep.mat');
m = matfile(tfilename,'Writable',true);
h = waitbar(0,'Initializing');
for i = 1:n_tracks2p:n_tracks
    if n_tracks-i < n_tracks2p
        tlen = n_tracks-i;
    else
        tlen = n_tracks2p;
    end
    
    tracks = zeros(tlen, n_slices,datatype);
    for t = 1:tlen
        waitbar((i+t-1)/n_tracks,h,'Reparsing adjacency')
        
        adjacency_track = adjacency_tracks{i+t-1};
        for j = 1 : numel(adjacency_track)
            
            cell_index = adjacency_track(j);
            
            % We must determine the frame this index belong to
            tmp = cell_index;
            frame_index = 1;
            while tmp > 0
                tmp = tmp - n_cells(frame_index);
                frame_index = frame_index + 1;
            end
            frame_index = frame_index - 1;
            in_frame_cell_index = tmp + n_cells(frame_index);
            
            tracks(t,frame_index) = in_frame_cell_index;            
        end        
    end
    waitbar((i+t-1)/n_tracks,h,'Saving')
    m.tracks(i:tlen+i-1,1:n_slices) = tracks;
end
clear m
close(h)

end