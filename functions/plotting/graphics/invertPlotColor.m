function invertPlotColor(varargin)
%%
persistent state
if nargin==1
    switch varargin{1}
        case {'k','black'}
            colordef black
            set(gcf,'Color','k')
        case {'w','white','default'}
            colordef white
            set(gcf,'Color','w')
    end
else
    if isempty(state)||~state
        state = 1;
        colordef black
        set(gcf,'Color','k')
    elseif state
        state = 0;
        set(gcf,'Color','w')
        colordef white
    end
end
end

