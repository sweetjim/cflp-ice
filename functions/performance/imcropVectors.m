function [x,y] = imcropVectors(x,y,roi)
%IMCROPVECTORS applies the 2D ROI to the corresponding vectors
roix = roi(1):(roi(3)+roi(1));
roiy = roi(2):(roi(4)+roi(2));
if max(roix)>numel(x)
    roix = roix(1):numel(x);
end
if max(roiy)>numel(y)
    roiy = roiy(1):numel(y);
end
x = x(roix);
y = y(roiy);
end

