function [map,values] = roundto(array,value)
%% ROUNDTO returns the values and indexing map from an array that is rounded to the set value
% Outputs: 
%   map     - index locations
%   values  - list of integer values in array
%
% Example:
%   x = linspace(0,1);
%   [map,vals] = roundto(x,0.5)
%   map     = [1 26 76]
%   vals    = [0 0.5 1]
%%
[values,map]=unique(floor(array/value)*value);
end

