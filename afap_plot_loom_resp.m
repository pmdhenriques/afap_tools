function afap_plot_loom_resp(aq)
%% Plots response rate and size (latency), separating between hunt-on and off conditions
%
% Pedro Henriques, Jan 2018

ori = unique(aq.loomz(:,14));
RR = NaN(2,length(ori));

honix = any(aq.loomz(:,8:9),2);
rix = ~isnan(aq.loomz(:,6));

figure
for o = 1:length(ori)
    oix = aq.loomz(:,14) == ori(o);
    RR(1,o) = sum(~honix & rix & oix) / ...
        sum(~honix & oix);  % Not hunting
    RR(2,o) = sum(honix & rix & oix) / ...
        sum(honix & oix);   % Hunting
end
bar(RR')
ylim([0 1])
set(gca,'XTickLabel',string(ori))
ylabel('Response rate')
xlabel('Loom laterality')
legend({'~Hunting','Hunting'})

figure
for o = 1:length(ori)
    oix = aq.loomz(:,14) == ori(o);
        
    ix1 = ~honix & rix & oix;
    ix2 = honix & rix & oix;
    L1 = aq.loomz(ix1,7); if isempty(L1); L1 = nan; end
    L2 = aq.loomz(ix2,7); if isempty(L2); L2 = nan; end    
    
    subplot(1,length(ori),o)
    boxplot(padcat(L1,L2))
    hold on
    if any(ix1)
        plot(1,aq.loomz(ix1,7),'bo')
    end
    if any(ix2)
        plot(2,aq.loomz(ix2,7),'ro')
    end
    ylim([0 120])
    set(gca,'XTickLabel',{'~Hunt','Hunt'})
    title(sprintf('Ori = %d',ori(o)))
    if o == 1
        ylabel('Loom size @ response (deg)')
    end
end
end