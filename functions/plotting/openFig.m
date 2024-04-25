function openFig(fig_title)
figs = findall(groot,'Type','Figure');
for i=1:length(figs)
    if strcmp(figs(i).Name,fig_title)
       set(0,'CurrentFigure',figs(i))
       shg
    end
end
end

