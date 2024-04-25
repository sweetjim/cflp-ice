function addlabels(varargin)
%% ADDLABELS is a function for efficient labelling of axes
% Valid NAME-VALUE input arguments:
%   'ax'        [mxn axes/tiledlayout array]
%   'x'         [1xn cell array, string]
%   'y'         [1xn cell array, string]
%   'z'         [1xn cell array, string]
%   'c'         [1xn cell array, string]
%   'title'     [1xn cell array, string]
%   'fs'        [double]
% Valid VALUE arguments
%   'latex'     sets axes and colorbar interpreter to latex
%   'append'    appends / replaces labels existing labels
% ________________________________________________________________________________
% ________________________________________________________________________________
% EXAMPLES:
% ________________________________________________________________________________
%       x = linspace(0,1);
%       y = x.^2;
%       z = x;
%       plot3(x,y,z)
%
%  CASE 1: Label current axes
%       addlabels('x','Time (s)','y','Velocity (m/s)','fs',20)
%
%  CASE 2: Append current axes
%       addlabels('y','$U$ (ms$^{-1}$)','latex')
% ________________________________________________________________________________
%  CASE 3: Label 1x3 array of axes with unique titles
%       TL=tiledlayout(1,3);
%           for i=1:3
%           nexttile;
%           plot(x,x.^i)
%       end
%       addlabels(TL,'x','$x$','y','$f(x)$','title',{'$x^1$','$x^2$','$x^3$'},'latex')
%
%  CASE 4: Label 2x3 array of axes with unique titles and array of x-labels
%  and y-labels
%       TL=tiledlayout(2,3);
%           for i=1:6
%           nexttile;
%           plot(x,x.^i)
%           addlabels(...
%           'x','$x$',...
%           'y','$f(x)$',...
%           'title',sprintf('$%s^%i$','x',i),...
%           'latex')
%       end
%       addlabels(TL,...
%           'x',{'Row 1','Row 2'},...
%           'y',{'Column 1','Column 2','Column 3'},...
%           'array',...
%           'latex')
%
%
% ________________________________________________________________________________
% ________________________________________________________________________________

%%
appending   = false;
ax          = [];
x_str       = '';
y_str       = '';
z_str       = '';
c_str       = '';
title_str   = '';
fs          = 10;
latex       = false;
useArray    = false;
parseInput(varargin);
if useArray
    if ~isa(ax,'matlab.graphics.layout.TiledChartLayout')
        error('Expected "ax" input is not a TiledLayout class')
    end
    parent  = ax;
    dims    = ax.GridSize; % [row col]
    ax      = ax.Children;
    if numel(ax)==1
        error('Expected input for an "array" is missing an array of axes.')
    end

    if ~iscellstr(x_str)
        x_str = {x_str};
    end
    if ~iscellstr(y_str)
        y_str = {y_str};
    end
    dim_str = [numel(y_str) numel(x_str)]; % [row col]

    if any(dims~=dim_str)
        incorrect   = {'y-label','x-label'};
        fixStr      = [incorrect{(dims~=dim_str)} ' array'];
        if all(dims~=dim_str)
            fixStr = 'both x-label and y-label arrays';
        end
        error('Incompatible array size with TiledLayout grid.\n Check size of %s',fixStr)
    end

    layout = arrayfun(@(x) parent.Children(x).Layout.Tile,1:numel(ax));
    ax      = flipud(ax);
    x_str   = fliplr(x_str);

    tileX = max(layout)-(0:dims(1)-1);
    tileY = dims(1).*(0:dims(2))+1;
    tileY = tileY(tileY<max(layout));

    if latex
        arrayfun(@(x) addlabels('ax',ax(tileX(x)),'x',x_str(x),'latex'),1:numel(tileX))
        arrayfun(@(x) addlabels('ax',ax(tileY(x)),'y',y_str(x),'append','latex'),1:numel(tileY))
    else
        arrayfun(@(x) addlabels('ax',ax(tileX(x)),'x',x_str(x)),1:numel(tileX))
        arrayfun(@(x) addlabels('ax',ax(tileY(x)),'y',y_str(x),'append'),1:numel(tileY))
    end
    return
end

if isempty(ax)
    ax = gca;
end
set(findall(ancestor(ax(1),'Figure'),'Type','Colorbar'),'UserData',[])

if isa(ax,'matlab.graphics.layout.TiledChartLayout')
    ax = findall(ax,'Type','Axes');
end

skipCbar = false;
if ~iscell(c_str)
    if isempty(c_str)
        skipCbar = true;
    end
    c_str = {c_str};
end

% check for latex script
if any(contains({x_str,y_str,z_str},'$'))
    latex = true;
end


for i=1:numel(ax)
    axi = ax(i);
    if ~skipCbar
        ci  = findAssociatedColorbar(axi);
        c_str0 = c_str;
        if numel(c_str)==numel(ax)
            c_str0 = c_str{i};
        end
        if ~strcmp(c_str0,'')
            ci.Label.String = c_str0;
            ci.Label.Interpreter    = 'tex';
            ci.TickLabelInterpreter = 'tex';
            if latex
                ci.Label.Interpreter    = 'latex';
                ci.TickLabelInterpreter = 'latex';
            end
            ci.FontSize = fs;
        end
    end

    if numel(title_str)>1 && isa(title_str,'cell')
        title_stri = title_str{i};
    else
        title_stri = title_str;
    end

    if appending
        fs = axi.FontSize;
        if ~isempty(x_str),if latex,try xlabel(axi,x_str,'Interpreter','latex'),end,else,xlabel(axi,x_str),end,end
        if ~isempty(y_str),if latex,try ylabel(axi,y_str,'Interpreter','latex'),end,else,ylabel(axi,y_str),end,end
        if ~isempty(z_str),if latex,try zlabel(axi,z_str,'Interpreter','latex'),end,else,zlabel(axi,z_str),end,end
        if ~isempty(title_str),if latex,try title(axi,title_str,'Interpreter','latex'),end,else,title(axi,title_str),end,end
        set(axi,'FontSize',fs)
        continue
    end

    if latex
        try %#ok<*TRYNC>
            xlabel(axi,x_str,'Interpreter','latex')
            ylabel(axi,y_str,'Interpreter','latex')
            title(axi,title_stri,'Interpreter','latex')
            set(axi,'TickLabelInterpreter','latex')
            zlabel(axi,z_str,'Interpreter','latex')
            ci.Label.Interpreter     = 'latex';
            ci.TickLabelInterpreter  = 'latex';
        end
    else
        try
            xlabel(axi,x_str,'Interpreter','tex')
            ylabel(axi,y_str,'Interpreter','tex')
            title(axi,title_stri,'Interpreter','tex')
            set(axi,'TickLabelInterpreter','tex')
            zlabel(axi,z_str,'Interpreter','tex')
            ci.Label.Interpreter     = 'tex';
            ci.TickLabelInterpreter  = 'tex';
%             xlabel(axi,x_str)
%             ylabel(axi,y_str)
%             title(axi,title_stri)
%             zlabel(axi,z_str)
        end
    end
end
try
    set(ax,'FontSize',fs)
catch
    if isa(ax,'matlab.graphics.layout.TiledChartLayout')
        tiletitle = get(ax,'Title'); set(tiletitle,'FontSize',fs);
        tilex = get(ax,'Xlabel'); set(tilex,'FontSize',fs);
        tiley = get(ax,'Ylabel'); set(tiley,'FontSize',fs);
    end
end


%% Input parser
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        [ax,items]= axesCheck(items);
        for k=1:length(items)
            switch items{m}
                %% Name arguments
                case 'latex'
                    latex   = true;
                    %% Name-value arguments
                case 'ax'
                    ax      = namevalue;
                case 'title'
                    title_str   = namevalue;
                case {'x','x_str'}
                    x_str   = namevalue;
                case {'y','y_str'}
                    y_str   = namevalue;
                case {'z','z_str'}
                    z_str   = namevalue;
                case {'c','c_str'}
                    c_str   = namevalue;
                case 'fs'
                    fs      = namevalue;
                case 'append'
                    appending = true;
                case 'array'
                    useArray = true;
            end
            m = m+1;
            if m>length(items);break;end
        end
        function out = namevalue
            out = items{m+1};
            m   = m+1;
        end
    end
%% Additional functions
    function c = findAssociatedColorbar(ax)
        ap = ax.Position;
        cb = findall(ancestor(ax,'Figure'),'Type','Colorbar');
        cb = setxor(cb,findall(cb,'UserData','tagged'));
        sameXOrigin = arrayfun(@(x) ( ...
            within(ap(1),padding(x.Position(1),.1))) ...
            ,cb);
        sameYOrigin = arrayfun(@(x) ( ...
            within(ap(2),padding(x.Position(2),.1))) ...
            ,cb);
        sameOrigin = any([sameXOrigin sameYOrigin],2);

        xr = padding(cumsum(ap([1 3])),.1);
        yr = padding(cumsum(ap([2 4])),.1);
        withinAxes = arrayfun(@(x) ...
            any(within(cumsum(x.Position([1 3])),xr))&& ...
            any(within(cumsum(x.Position([2 4])),yr)) ...
            ,cb);

        
        if any(sameOrigin)&&~all(sameOrigin)
            c = cb(sameOrigin);
        elseif any(withinAxes)
            c = cb(withinAxes);
        end
        c.UserData = 'tagged';
        function values = padding(values,percent)
            values = values.*[1-percent 1+percent];
        end
    end
    function out = within(x,xrange,type)
        %% Within
        % Returns a logical argument that depends on whether 'x' falls within the
        % range of 'xrange'. Can also specify inclusive (DEFAULT i.e., <=) or
        % exclusive bounds (i.e., <)
        %
        % Examples:
        % >> within(0,[-1 1])
        % >> ans =
        % >>    logical
        % >>        1
        % --------------------
        % >> within(5.1,0:5)
        % >> ans =
        % >>    logical
        % >>        0
        % --------------------
        % >> within(1,1:2,'exclusive')
        % >> ans =
        % >>    logical
        % >>        0
        %%
        x1 = xrange(1);
        x2 = xrange(end);
        if nargin<3
            type = 'inclusive';
        end
        switch type
            case 'inclusive'
                out = (x>=x1&x<=x2);
                if diff(xrange)<0
                    out = (x<=x1&x>=x2);
                end
            case 'exclusive'
                out = (x>x1&x<x2);
                if diff(xrange)<0
                    out = (x<x1&x>x2);
                end
        end
    end
    function [ax,items] = axesCheck(items)
        %% AXESCHECK
        % Script used within parseInput to determine if first argument is an axes
        hasAxArg = cell2mat(cellfun(@(x) contains(x,'ax'),items,'ErrorHandler',@(~,~) false,'UniformOutput',false));
        hasAxObj = cellfun(@(x) isa(x,'matlab.graphics.axis.Axes'),items);
        hasTLObj = cellfun(@(x) isa(x,'matlab.graphics.layout.TiledChartLayout'),items);
        ax = [];
        if ~any(hasAxArg)
            if any(hasAxObj)
                ax      = items{hasAxObj};
                items   = items(~hasAxObj);
            elseif any(hasTLObj)
                ax      = items{hasTLObj};
                items   = items(~hasTLObj);

            end
            return
        end
    end
end