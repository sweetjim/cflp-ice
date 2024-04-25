function ax = subAxes(parent,index,varargin)
%SUBAXES takes a TiledLayout object and creates a nested array of axes
%using the same arguments as TiledLayout.
%
%Usage:
%tile   = tiledlayout(1,3);
%a1     = subAxes(tile,1,2,2,'TileSpacing','none');
%a2     = subAxes(tile,2,1,2,'TileSpacing','none');
%a2     = subAxes(tile,3,3,2,'TileSpacing','none');


if ~isa(parent,'matlab.graphics.layout.TiledChartLayout')
    error('Parent must be a TiledLayout object')
end
tile = tiledlayout(parent,varargin{:});
tile.Layout.Tile = index;
ax = arrayfun(@(idx) nexttile(tile,idx),1:prod(tile.GridSize));
end

