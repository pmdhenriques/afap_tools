function [] = afap_paramecia_tracks_join(datadir,tracks2j)
% Joins two or more paramecia tracks together. Usually called from within
% afap_paramecia_score.m script.
%
% Pedro Henriques


m = matfile(fullfile(datadir,'tracks.mat'),'Writable',true);

vclass = class(m.tracks(1,1));

t = zeros(size(tracks2j,2),size(m.tracks(1,:),2),vclass);
for i = 1:size(tracks2j,2)
    t(i,:) = m.tracks(tracks2j(i),:);
end
tj = sum(t,1,'native');

for i = 1:size(tracks2j,2)
    if i == 1
        m.tracks(tracks2j(i),:) = tj;
    else
        m.tracks(tracks2j(i),:) = zeros(1,size(t,2),vclass);
    end
end
clear m

disp('Great success!')

end
