function [map,values] = findintegers(array)
%% FINDINTEGERS returns the integer values and indexing map from an array
% Outputs: 
%   map     - index locations
%   values  - list of integer values in array
%%
[values,map] = unique(floor(array));

end

