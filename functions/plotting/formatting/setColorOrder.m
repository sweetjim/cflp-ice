function setColorOrder(ax,cmap)
%SETCOLORORDER Summary of this function goes here
%   Detailed explanation goes here
h = ax.Children;
% if reverse
%     h = flipud(h);
% end
for i=1:length(h)
    isline(i) = isa(h(i),'matlab.graphics.chart.primitive.Line');  %#ok<AGROW>
end
h = h(isline);
C = imresize(cmapCheck(cmap),[length(h) 3],'nearest');
for i=1:length(h)
    h(i).Color =  C(i,:);
end
end

