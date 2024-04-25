function output = rescaleMax(input,absolute)
%RESCALEMAX simply divides the input by its unweighted matrix maximum
if nargin<2
    absolute = false;
end
maxval = max(input,[],'all','omitnan');
if absolute
    maxval = max(abs(input),[],'all','omitnan');
end
output = input./maxval;

end

