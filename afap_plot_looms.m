function afap_plot_looms(aq,idxs)
% plots figures with tail, eye and loom size features for aq structures

offset = 700;

if nargin < 2
    idxs = 1:size(aq.loomz,1);
elseif size(idxs,1) > size(idxs,2)
    idxs = idxs';
end

for i = idxs
    stix = aq.loomz(i,1)-offset;
    edix = aq.loomz(i,2)+offset;
    L = sgolayfilt(aq.Langl(stix:edix),2,21);
    R = sgolayfilt(aq.Rangl(stix:edix),2,21);
    
    figure
    plot(L)
    hold on
    plot(R)
    plot(aq.Verg(stix:edix))
    plot(aq.hunton(stix:edix).*aq.pthrON,'k-')
    plot(aq.hunton(stix:edix).*aq.pthrON,'k-')
    line([0 length(L)],[55 55],'LineStyle','--','Color','g')
    ylim([-50 80])
    yyaxis right
    plot(aq.fvelV(stix:edix))
    
    thisLv = aq.loomz(i,5);
    [F, ~] = afap_loomingstim(thisLv, aq.settings.loomST, aq.settings.loomED, aq.settings.DISPLAYpxpermm, aq.settings.DVIEWmm, aq.settings.datarate);
    plot(offset:length(F)+offset-1, 2.*atand(F./(aq.settings.DISPLAYpxpermm.*aq.settings.DVIEWmm)), 'LineWidth', 2, 'Color','k');
    plot(aq.loomz(i,6)+offset,aq.loomz(i,7),'k*','MarkerSize',10)
    ylim([0 120])
    
    if aq.settings.analysisversion >= 180119
        title(mat2str(aq.loomz(i,[5,12,14:16])))
    else
        title(mat2str(aq.loomz(i,[5,14])))
    end
end
end