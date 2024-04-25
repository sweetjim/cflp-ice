function cmap = cmapCheck(cmap)
f = figure('Visible','off');
switch colormapclasses(cmap)
    case 'cmocean'
        cmap = cmocean(cmap);
    otherwise
        cmap = colormap(cmap);
        return
end
delete(f)
end