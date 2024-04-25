function addtick(varargin)
%% Add Tick
% Adds a colorbar to the specified axes (or current axes) with a title
% formatted in LaTeX or in plain script.
% -------------------------------------------------------------------------
% Example
% -------------------------------------------------------------------------
% 
% addtick('x',0.1,'color','r')
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%%

ticklabel = [];
[cmap,ax,dir,color,newticks] = parseInput(varargin);
if numel(ticklabel)==1;ticklabel={ticklabel};end
tick_axes = fieldnames(newticks);

for i=1:length(tick_axes)
   switch tick_axes{i}
       case 'x'
           ticks        = newticks.x;
           tick_dir     = 'Xtick';
           tick_color   = 'XColor';
       case 'y'
           ticks        = newticks.y;
           tick_dir     = 'Ytick';
           tick_color   = 'YColor';
       case 'z'
           ticks        = newticks.z;
           tick_dir     = 'Ztick';
           tick_color   = 'ZColor';
   end
   ticklabel0 = get(ax,[tick_dir 'label']);
   [update_ticks,idx] = unique(cat(2,get(ax,tick_dir),ticks));
   ticklabel0 = unique(cat(1,ticklabel0,ticklabel));
   ticklabel0 = ticklabel0(idx);
   set(ax,tick_dir,update_ticks,[tick_dir 'label'],ticklabel0)
   
   if ~isempty(color)
      if numel(color)==1
          for j=1:length(ticks)
          pos(j) = find(update_ticks==ticks(j),1);
          end
      else
          
      end
   end
   

   switch tick_axes{i}
       case 'x'
           ax.XTickLabel{3} = 2;
       case 'y'
           
       case 'z'
           
   end
end


function [cmap,ax,dir,color,newticks] = parseInput(varargin)
        cmap        = 'amp';
        ax          = gca;
        dir         = 'x';
        color       = [];
        newticks    = struct([]);
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                %% Name arguments
                case colormapclasses
                    cmap    = items{m};               
                case isa(gca,'matlab.graphics.axis.Axes')
                    ax      = items{m};
                %% Name-value arguments
                case {'x','xtick'}
                    dir          = 'x';
                    newticks(1).x   = namevalue;
                case {'y','ytick'}
                    dir          = 'y';
                    newticks(1).y   = namevalue;
                case {'z','ztick'}
                    dir          = 'z';
                    newticks(1).z   = namevalue;
                case 'label'
                    ticklabel = namevalue;
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

