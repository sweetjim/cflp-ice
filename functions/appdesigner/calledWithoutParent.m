function [state,parent] = calledWithoutParent(comp,name,varargin)
p = inputParser;
addOptional(p,'Icon','')
parse(p,varargin{:})
icon        = p.Results.Icon;

state = false;
parent = comp.Parent;
if isa(comp.Parent,'matlab.ui.Figure')
    ss = get(0,'ScreenSize');
    comp.Parent.Position = [[ss([3 4])/2-comp.Position([3 4])/2] comp.Position([3 4])]; %#ok<NBRAK2>
    comp.Parent.Name = name;
    if ~isempty(icon)
        if isfile(icon)
            comp.Parent.Icon = icon;
        end
    end
    state = true;
end
if ~nargout
    clear state
end

