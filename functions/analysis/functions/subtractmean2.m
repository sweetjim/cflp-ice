function hout = subtractmean2(h)
hx = mean(h,1);
hy = mean(h,2);

Hx = repmat(hx,[size(h,1) 1]);
Hy = repmat(hy,[1 size(h,2)]);

hout = h-(Hx+Hy);
end

