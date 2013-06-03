function [] = visualize_bin_probs(t, bins, probs, tnow, maxtime)

hold on
axis([0, maxtime, 0, numel(bins)])
for i = 1:numel(bins)
    maxprob = 1.1*max([probs{i,1}, probs{i,2}]);
    plot(t,probs{i,1}/maxprob-i+6,'r')
    plot(t,probs{i,2}/maxprob-i+6,'b')
    plot(t,zeros(1,numel(t))-i+6,'k')
end

ylabels = {};
for i = 1:numel(bins)
    ylabels{i} = sprintf('Bin %d', bins(numel(bins)-i+1));
end
plot([tnow, tnow], [0, numel(bins)], 'g');

AX = gca;
set(AX,'YTick',(1:numel(bins))-0.5);
set(AX,'YTickLabel',ylabels);
set(AX,'YGrid','off');
set(AX,'XGrid','on');
