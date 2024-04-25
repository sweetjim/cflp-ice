function flatarray = flatten(array)
%FLATTEN reduces a 2D array into a 1D array
flatarray     = reshape(array,[1 size(array,1)*size(array,2)]);
end

