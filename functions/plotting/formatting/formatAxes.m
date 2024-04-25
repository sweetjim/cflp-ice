classdef formatAxes<handle
    %FORMATAXES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fig
        ax
        tiled
        interpreter
        fontsize
        font
        savepath
        view
    end
    
    methods
        %% CONSTRUCTOR
        function obj = formatAxes(axes,varargin)
            % Constructor
            %
            %%
            obj.savepath    = cd;
            obj.interpreter = 'none';
            obj.fontsize    = 15;
            obj.font        = 'Calibri';
            obj.view        = [];
            parseInput(varargin)
            
           
            
            % If treenode selection (i.e. app_analysis_comparison.mlapp)
            if contains(class(axes),'Node')
%                idxclass = contains(arrayfun(@(x) class(x.NodeData),axes,'UniformOutput',false),{'Tile','Axes'});
               axes = axes.NodeData;
            end
            
            if contains(class(axes),'Figure')
               obj.fig      = axes;
               axes         = axes.Children;
            end
            
            obj.tiled       = [];
            obj.ax          = axes;
            
            % Check if axes is a tiledlayout object
            if contains(class(obj.ax),'Tile')
                obj.tiled    = axes;
                obj.ax       = axes.Children;
                
                % Check for colorbar
            end
            
            
            
            %% Input parser
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case {'default','latex'}
                            obj.interpreter = items{m};
                        case listfonts
                            obj.font = items{m};
                        case {'sp','savepath'}
                            obj.savepath   = namevalue;
                        case {'fontsize','fs','FontSize'}
                            obj.fontsize = namevalue;
                        case {'tl','tile','TiledLayout','tiled','tiledlayout'}
                            obj.tiled   = namevalue;
                        case 'view'
                            obj.view = namevalue;
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
        %% MAIN FUNCTIONS
        function parse(obj)
            %%
            items = {'ax','tile'};
            if isempty(obj.tiled)
               items = {'ax'};
            end
            
            %% FONT CHANGES
            for i=1:numel(items)
                switch items{i}
                    case 'ax'
                        active = obj.ax;
                    case 'tile'
                        active = obj.tiled;
                end
                
                % Apply font changes
                if isa(active,'matlab.graphics.layout.TiledChartLayout')
                    tiletitle = get(active,'Title'); set(tiletitle,'FontSize',obj.fontsize);
                    tilex = get(active,'Xlabel'); set(tilex,'FontSize',obj.fontsize);
                    tiley = get(active,'Ylabel'); set(tiley,'FontSize',obj.fontsize);
                else
                    set(active,...
                    'FontName',obj.font,...
                    'FontSize',obj.fontsize)
                    active = active(arrayfun(@(x) contains(class(x),'Axes'),active)); % get axes
                    obj.setview;
                end
                
                
                
                    % Change strings to accomodate fonts
                    %%
                    for j=1:numel(active)
                        %%
                        label = {'XLabel','YLabel','ZLabel','Title'};
                        for k=1:numel(label)
                            try
                                tmp = eval(sprintf('active(j).%s.String;',label{k}));
                            catch
                                continue
                            end
                            if isempty(tmp)
                                continue
                            end
                            if contains(tmp,'\textsf')
                               tmp = extractAfter(tmp,'\textsf');
                               tmp(1) = '';
                               tmp(end) = '';
                            end
                            if numel(tmp)==1
                                tmp = tmp{1};
                            end
                            if ~contains(obj.font,'Helvetica')
                               tmp = sprintf('%s{%s}','\textsf',tmp);
                            end
                            eval(sprintf('active(j).%s.String="%s";',label{k},tmp))
                            eval(sprintf('active(j).%s.Interpreter="%s";',label{k},obj.interpreter))
                        end
                    end
                
            end
            
        end
        function save(obj)
            figure(obj.fig)
            out2latex(obj.savepath)
        end
        function setview(obj)
            if ~isempty(obj.view)
                try %#ok<TRYNC>
                    set(findobj(obj.ax,'Type','axes'),'View',obj.view)
                end
            end
        end
    end
    methods (Hidden)
        
    end
end

