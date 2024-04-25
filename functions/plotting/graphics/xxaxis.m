function [AXES] = xxaxis(x1,y1,x2,y2,ax)
%XXAXIS Summary of this function goes here
%   Detailed explanation goes here
if nargin<5
    ax = gca;
end
ax1 = axes('Position',ax.Position);
plot(ax1,x1,y1,'-r')
ax1.XColor = 'r';
ax1.YColor = 'r';
ax2 = axes('Position',ax1.Position);
plot(ax2,x2,y2,'-k')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
ax2.Position = ax1.Position;
if nargout>0
    AXES = [ax1 ax2];
end
end

