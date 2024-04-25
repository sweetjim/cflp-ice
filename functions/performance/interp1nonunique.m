function yq = interp1nonunique(x,y,xq)
%INTERP1NONUNIQUE interpolates non-unique valued x-coordinates by
%introducing a floating point argument.
yq = interp1(cumsum(ones(size(x)))*eps + x,y,xq);
end

