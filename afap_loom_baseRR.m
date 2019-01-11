function [] = afap_loom_baseRR(datadir)
% Calculates the baseline escape response rate (measured using the same
% parameters as for escape responses during looming stimuli
%
% Pedro Henriques, Jul 2018

load(fullfile(datadir,'aq'),'aq')

wnd = mode(aq.loomz(:,3));  % Loom window size
nbts = size(aq.bouts,1);    % Number of bouts
isst = false(nbts,1);   % Bouts where stimulus being presented
for i = 1:nbts
    isst(i) = any(aq.bouts(i,1) >= aq.visstim(:,3) & ...
        aq.bouts(i,1) <= aq.visstim(:,4));
end

sttmsum = sum(aq.visstim(:,4)-aq.visstim(:,3)); % Total frames where stim was present
tottm = size(aq.data,1);    % Total number of frames
remtm = tottm-sttmsum;  % Total number of frames without stim

ix = aq.bouts(~isst,7) >= aq.settings.escapethr;    % Bouts where fish escapes
lmbasert = sum(ix)/(remtm/wnd); % Loom base response rate

aq.loom_base_rr = lmbasert;

try
    save(fullfile(datadir,'aq'),'aq')
catch
    save(fullfile(datadir,'aq'),'aq','-v7.3')
end
end