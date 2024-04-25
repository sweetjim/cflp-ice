function [range,range1d] = array2range(vector,N)
%ARRAY2RANGE converts a 1D vector of indicies (from an array of length N)
%into a 2D range array


if isrow(vector)
    vector = vector';
end
if nargin==2
    vector = sort([1 N vector'])';
end 
range        = repmat(vector,[1 2]);
range(:,2)   = circshift(range(:,2),-1)-1;
range(end,:) = [];

range1d = sort(reshape(range,1,[]));
end

