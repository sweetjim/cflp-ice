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
        ax      = ax.Children;
        items   = items(~hasTLObj);
    end
    return
end


