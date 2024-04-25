function latexformat(varargin)
%% Latex-format
% Formats the size of the figure to be a fraction of the standard linewidth (A4)
% 
% -------------------------------------------------------------------------
% %  Name-value arguments (Optional):
% -------------------------------------------------------------------------
%   'fig'   :   [figure]    - desired figure (defualt is gcf)
%   'h'     :   [double]    - height fraction of the linewidth
%   'w'     :   [double]    - width fraction of the linewidth
%   'square':   [char]      - sets the height to equal the width
%   'ratio' :   [vec]       - sets the height to width ratio
%       [w h]               - vectorised ratio
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%
h       = [];
w       = [];
ratio   = [];
fig     = gcf;
parseInput(varargin)

linewidth = 5.6188; %[in]

set(fig,'Units','inches')
pos = get(fig,'Position');
pos([1 2]) = 1;
if sum([h w ratio])==0
    %% Default - set the figure width to be the linewidth
    set(fig,'Position',[pos(1) pos(2) linewidth pos(4)]);
elseif isempty(ratio)
    if isempty(h)
        h =  1;
    elseif isempty(w)
        w = 1;
    end
    set(fig,'Position',[pos(1) pos(2) linewidth.*[w h]]);
else
    set(fig,'Position',[pos(1) pos(2) linewidth.*[ratio(1) ratio(2)]]);
end
    
%% Input parser
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'h','y'}
                    h = namevalue;
                case {'w','x'}
                    w = namevalue;
                case 'square'
                    ratio = 1;
                case 'ratio'
                    ratio = namevalue;
                case 'fig'
                    fig  = namevalue;
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

