function [out,outX,outY] = imcrop2(in,X,Y,roi)
%% IMCROP2 is an extension of <a href="matlab:help('imcrop')">imcrop</a>
% that also crops the axis vectors.
%%
out = imcrop(in,roi);

DX = roi(1):roi(3)+1;
DY = roi(2):roi(4)+1;

outX = X(DX);
outY = Y(DY);
end

