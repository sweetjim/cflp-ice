function newtick(ax,value,string,axis)
%% New tick marker
% Inserts a new tick at the prescribed value for the appropriate axis.
% -------------------------------------------------------------------------
% %  Singular arguments:
% -------------------------------------------------------------------------
%   ax  :   [axes]      -   specified axes
%   value : [double]    -   value of new tick
%   string: [char]      -   new tick label
%   axis  : [char]      -   specific axis
%       'x','y','z'
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%
switch axis
    case 'x'
        tick        = xticks;
        ticklabels  = xticklabels;
    case 'y'
        tick        = yticks;
        ticklabels  = yticklabels;
    case 'z'
        tick        = zticks;
        ticklabels  = zticklabels;
end

tick0       = tick;
tick(end+1) = value;
[tick,idx]  = sort(tick);

if length(unique(tick))==length(tick0)
    pos     = find(diff(tick)==0,1,'first');
    tick    = tick0;
    ticklabels = {ticklabels{1:pos-1},string,ticklabels{pos+1:end}}';
else
    pos         = find(idx==length(tick));
    ticklabels = {ticklabels{1:pos-1},string,ticklabels{pos:end}}';
end


switch axis
    case 'x'
        set(ax,'XTick',tick,'XTickLabel',ticklabels)
    case 'y'
        set(ax,'YTick',tick,'YTickLabel',ticklabels)
    case 'z'
        set(ax,'ZTick',tick,'ZTickLabel',ticklabels)
end
end

