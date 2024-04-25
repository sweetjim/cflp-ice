function output = array(input)
% ARRAY converts a [1 2] index vector into a [1 n] index vector such that 
% numel(output) = input(1):input(2)

output = input(1):1:input(2);