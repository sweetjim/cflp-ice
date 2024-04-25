function out = longestCurve(in,n)

if ischar(n)&&strcmp(n,'all')
    biggest_n   = 1e3;
elseif ~ischar(n)
    biggest_n   = n;
else
    biggest_n   = 1;
end

CC          = bwconncomp(in,8);
numPixels   = cellfun(@numel,CC.PixelIdxList);
[~,idx]     = sort(numPixels,'descend');


if biggest_n>length(idx)
    biggest_n=length(idx);
end

list    = idx(1:biggest_n);


for i=1:CC.NumObjects
    if ~sum(i==list)>0
        in(CC.PixelIdxList{i}) = 0;
    end
end

out = in;
end

