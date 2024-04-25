function output = rmse(sample,real)
%RMSE calculates the root-mean-square-error between two arrays
%%
output = sqrt((sample-real).^2);
end

