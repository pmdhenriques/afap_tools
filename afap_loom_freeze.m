function [] = afap_loom_freeze(datadir)
% Finds freezing responses for looming stimuli presentation
%
% 1 - start it
% 2 - freeze ibi length (frames)
% 3 - freeze ibi length (s)
% 4 - latency (frames
%
% Pedro Henriques, Jan 2018

%%

load(fullfile(datadir,'aq'),'aq');
load('Z:\Pedro\NI_project\Ablations\NI\Bilateral\analysis\freeze_thr.mat','frzthr')

for i = 1:size(aq.loomz,1)
    loomst = aq.loomz(i,1);
    btix = find([aq.btf.Loombout] == i);
    if ~isempty(btix)
        % check if end of last bout is within loom boundary
        if btix(1)-1 > 0
            lsted = aq.btf(btix(1)-1).ited;
            if lsted(1) > loomst
                btix = [btix(1)-1 btix];
            end
        end
        st = [aq.btf(btix).st];
        ed = [aq.btf(btix).ed];
        stfrm = [aq.btf(btix).itst];
        edfrm = [aq.btf(btix).ited];
        
        if length(st) >= 2
            ibi = st(2:end)-ed(1:end-1);
            ibifrm = stfrm(2:end)-edfrm(1:end-1);
            ibit = log2(ibi);
            frzix = find(ibit >= frzthr);
            if ~isempty(frzix)
                aq.loom_freeze(i,1) = edfrm(frzix(1));
                aq.loom_freeze(i,2) = ibifrm(frzix(1));
                aq.loom_freeze(i,3) = ibi(frzix(1));
            else
                aq.loom_freeze(i,1) = nan;
                aq.loom_freeze(i,2) = nan;
                aq.loom_freeze(i,3) = nan;
            end
        else
            aq.loom_freeze(i,1) = nan;
            aq.loom_freeze(i,2) = nan;
            aq.loom_freeze(i,3) = nan;
        end
    else
        aq.loom_freeze(i,1) = nan;
        aq.loom_freeze(i,2) = nan;
        aq.loom_freeze(i,3) = nan;
    end
end
aq.loom_freeze(:,4) = aq.loom_freeze(:,1)-aq.loomz(:,1);

% Save
try
    save(fullfile(datadir,'aq'),'aq');
catch
    save(fullfile(datadir,'aq'),'aq','-v7.3');
end
end