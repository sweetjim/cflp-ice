function mask = polyshape2mask(pshape,m,n,ds)
%POLYSHAPE2MASK converts a polyshape object to a mask
% Inputs:
%   pshape  -   [polyshape or vertices]
%   m       -   [int] number of rows
%   n       -   [int] number of columns
%   ds      -   [double] conversion (OPTIONAL, default = 1)
if nargin<4
    ds = 1;
end
if isa(pshape,'polyshape')
    mask = poly2mask(pshape.Vertices(:,1)*ds,pshape.Vertices(:,2)*ds,m,n);
    return
end
    mask = poly2mask(pshape(:,1)*ds,pshape(:,2)*ds,m,n);
end

