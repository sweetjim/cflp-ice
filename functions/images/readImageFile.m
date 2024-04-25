function image = readImageFile(imagepath,isBinary)
%READIMAGE Summary of this function goes here

%%
if nargin<3
    isBinary = false;    
end
image = double(imread(imagepath));

if length(size(image))>2 && isBinary
    image = imread(imagepath);
    image = rgb2gray(image);
end
end

