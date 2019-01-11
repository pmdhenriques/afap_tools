function [] = afap_paramecia_tracks(datadir)
% Script to connect paramecia objects into tracks from the afap pipeline.
% Needs to have run afap_paramecia_process before and have a objects.mat
% file in the input directory
%
% Pedro Henriques


max_linking_distance = 15;
max_gap_closing = 10;
debugmd = false;
linkmethod = 'Hungarian';

if nargin < 1
    datadir = uigetdir('\\128.40.155.187\data2\Bianco_lab\Pedro\NI project\Ablations\NI', ...
        'Select objects directory');
end

if ~exist(fullfile(datadir,'tracks.mat'),'file')
    objfile = dir(fullfile(datadir,'objects*.mat'));
    if ~isempty(objfile)
        objpath = fullfile(objfile(1).folder,objfile(1).name);
        load(objpath)
        
        [tracks,err] = simpletracker_ph({objects.centroid}, ...
            'MaxLinkingDistance',max_linking_distance, ...
            'MaxGapClosing',max_gap_closing, ...
            'Debug',debugmd, ...
            'Method',linkmethod);
        
        disp('Saving...')
        save(fullfile(datadir,'tracks.mat'),'tracks','-v7.3')
        
        if err{1} == 1
            afap_paramecia_tracks_reparse(datadir);
        end
                
    else
        errordlg('Objects file doesn''t exist')
        return
    end
else
    disp('tracks file already exists!')
end

end