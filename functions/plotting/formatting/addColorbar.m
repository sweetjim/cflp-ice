function varargout = addColorbar(varargin)
%% Add Colorbar
% Adds a colorbar to the specified axes (or current axes) with a title
% formatted in LaTeX or in plain script.
% -------------------------------------------------------------------------
% %  Singular arguments (Optional):
% -------------------------------------------------------------------------
%   lineplot: [char] (Optional)
%   If the axes is a lineplot with sequential colorodering.
%   Must additionally specify the c-axis limits and colormap (cmocean
%   supported). (Default is false, [0,1], 'gray').
%
%   cmap:   [char]
%   Specified colormap (see 'colormap' and 'cmocean' documentation) to be
%   used in conjunction with 'lineplot' (default is 'gray').
%
%   flip:   [char]
%   Flips the colormap (default is false).
%
%   ivert: [char]
%   Inverts the colormap (default is false).
%
%   reverse: [char]
%   Invert the order of lines (default is false).
%
%   invisible: [char]
%   Set the colorbar visibility to 'off' to leverage the custom colormap
%   when handling multiple axes (default is false) - note, the colorbar
%   will still exist!
%
%   colorbar: [char]
%   Do not create the colorbar - only the colormap (default is false).
%
%   log:    [char]
%   Set the colormap scale to logarithmic (defualt is linear).
%
%   latex:  [bool] (Optional)
%   Latex formatting (default is false).
%
%   location: [char]
%   Location of the colorbar (default is empty).
%       - 'north', 'south', 'east', 'west'
%       - 'northoutside', 'southoutside', 'eastoutside', 'westouside'
%
%
% -------------------------------------------------------------------------
% % Name-value arguments (Optional):
% -------------------------------------------------------------------------
%   linemap: [vector] 
%   Colormapper for "lineplot" (default is empty).
%
%   ax:     [axes] (Optional)
%   Axes to add colorbar to (default is gca).
%
%   label:  [char] (Optional)
%   Colorbar label (default is '').
%
%   fs:     [int] (Optional)
%   Fontsize (default is 10).
%
%   limits: [double, vec]
%   Two value argument that sets the limits of the colorbar axes. To be
%   used in conjunction with 'lineplot' (default is [0,1]).
%
%   pivot: [double]
%   Colormap pivot argument for 'cmocean' colormaps (default is []).
%
%   levels: [int]
%   Number of levels in the colormap (default is 256).
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%
type        = [];
invisible   = false;
cbar        = true;
flipcmap    = false;
invert      = false;
reverse     = false;
logscale    = false;
ax          = [];
label       = '';
fs          = 10;
latex       = false;
lineplot    = false;
linemap     = [];
cmap        = 'gray';
location    = [];
type        = 'default';
limits      = [];
pivot       = [];
levels      = [];
tile        = 1;
tilespan    = [1 1];

cmoceanmaps = {'thermal','balance','haline','delta','solar','curl'...
    'ice','diff','gray','tarn',...
    'oxy','deep','dense','phase',...
    'algae','matter','turbid','topo',...
    'speed','amp','tempo','rain','imfuse'};
diverging = {'balance','curl','diff','tarn'};
brewermaps = { 'BrBG','Accent','Blues','PuBuGn',...
    'PiYG','Dark2','BuGn','PuRd','PRGn','Paired','BuPu','Purples',...
    'PuOr','Pastel1','GnBu','RdPu',...
    'RdBu','Pastel2','Greens','Reds',...
    'RdGy','Set1','Greys','YlGn',...
    'RdYlBu','Set2','OrRd','YlGnBu',...
    'RdYlGn','Set3','Oranges','YlOrBr',...
    'Spectral','PuBu','YlOrRd'};
colormaps = {'parula','jet','hsv','hot','cool','spring',...
    'spring','summer','autumn','winter',...
    'gray','bone','copper','pink','lines',...
    'colorcube','prism','flag','white'};
parseInput(varargin);

if numel(ax)>1
    isTL = cellfun(@(x) isa(x,'matlab.graphics.layout.TiledChartLayout'),varargin);
    isAX = cellfun(@(x) isa(x,'matlab.graphics.axis.Axes'),varargin);
    varargin = varargin(~any([isTL;isAX]));
    arrayfun(@(x) addColorbar(x,varargin{:}),ax)
    return
end
if ~isempty(pivot)&&any(~contains(cmap,diverging))
    pivot = [];
end

if isempty(ax)
    ax = gca;
else
    if numel(ax)>1
        o1 = cellfun(@(x) isequal(x,'ax'),varargin);
        o2 = cellfun(@(x) isa(x,'matlab.graphics.axis.Axes'),varargin);
        varargin = varargin(~logical(o1+o2));
        arrayfun(@(x) addColorbar('ax',x,varargin{:}),ax)
        return
    end
end

if isempty(location)
    c = colorbar(ax);
else
    c = colorbar(ax,'Location',location);
    if strcmp(location,'layout')
        delete(findall(ax.Parent,'Type','Colorbar'))
        c = colorbar(ax,'Location',location);
        c.Layout.Tile = tile;
        c.Layout.TileSpan = tilespan;
    end
end
warning off
c.Label.String          = label;
warning on
if latex
    c.Label.Interpreter     = 'latex';
    c.TickLabelInterpreter  = 'latex';
end
c.FontSize = fs;
if isempty(levels)
    levels = 256;
end

if ~isempty(limits)
    caxis(ax,limits)
end
switch type
    case 'brewer'
        cmap = brewermap(levels,cmap);
    case 'cmocean'
        f = figure('Visible','off');
        if ~isempty(pivot)
            cmap = cmocean(cmap,'pivot',pivot,ax);
        else
            cmap = cmocean(cmap,ax);
        end
        delete(f)
    otherwise
        f = figure('Visible','off');
        cmap = colormap(cmap);
        delete(f)
end
cmap = interp1(1:size(cmap,1), cmap, linspace(1,size(cmap,1),levels),'linear');
if flipcmap
    cmap = flipud(cmap);
end
if invert
    cmap = 1-cmap;
end
if invisible
    c.Visible = 'off';
end
if lineplot
    % h = ax.Children;

    h = findall(ax,'Type','Line');
    if reverse
        h = flipud(h);
    end
    C = imresize(cmap,[length(h) 3],'nearest');

    if ~isempty(linemap)
        if numel(linemap)~=length(h)
            error('Incompatiable assignment. Linemap does not have the same dimensions as the line data.')
        end
        if any(gradient(linemap)<0)&&any(gradient(-linemap)>0)
            error('Linemap must be a monotonic increasing vector')
        end
        C    = interp1(linspace(linemap(1),linemap(end),size(C,1)),C,linemap,'linear');
        if isempty(limits)
            limits  = [linemap(1) linemap(end)];
        end
    end

    for i=1:length(h)
        h(i).Color =  C(i,:);
    end
    if isempty(limits)
        limits = [0,1];
    end
    caxis(ax,limits)
else
    if reverse
        cmap = flipud(cmap);
    end
end

assignin('base','cbar',c)
if ~cbar
    delete(c)
end
if logscale
    set(ax,'ColorScale','log')
end
set(ax,'Colormap',cmap)

if nargout>0
    varargout{1} = c;
end

%% Input parser
    function parseInput(varargin)

        m = 1;
        items = varargin{:};
        [ax,items] = axesCheck(items);
        for k=1:length(items)
            switch items{m}
                %% Name arguments
                case 'latex'
                    latex       = true;
                case {'line','lineplot'}
                    lineplot    = true;
                case 'linemap'
                    linemap     = namevalue;
                case cmoceanmaps
                    cmap        = items{m};
                    type        = 'cmocean';
                case brewermaps
                    cmap        = items{m};
                    type        = 'brewer';
                case colormaps
                    cmap        = items{m};
                case {'colorbar','colormap','nocolorbar'}
                    cbar        = false;
                case 'flip'
                    flipcmap    = true;
                case 'invert'
                    invert     = true;
                case {'reverse'}
                    reverse    = true;
                case {'log','logscale'}
                    logscale    = true;
                case 'linear'
                    logscale    = false;
                case {'invisible','tile','subplot'}
                    invisible   = true;
                    %% Name-value arguments
                case 'ax'
                    ax          = namevalue;
                case {'label','title'}
                    label       = namevalue;
                case 'fs'
                    fs          = namevalue;
                case {'lim','limit','limits'}
                    limits      = namevalue;
                case 'pivot'
                    pivot       = namevalue;
                case 'levels'
                    levels      = namevalue;
                case 'location'
                    location    = namevalue;
                    if strcmp(location,'layout')
                        m = m+1;
                        tile     = items{m};
                        m = m+1;
                        tilespan = items{m};
                    end
                        
            end
            m = m+1;
            if m>length(items);break;end
        end
        function out = namevalue
            out = items{m+1};
            m   = m+1;
        end
    end
end

