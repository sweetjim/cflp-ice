function [map0,xcoords,ycoords] = uncertainty2D(x,y,xu,yu,number_std)
% UNCERTAINTY2D computes the 2D uncertainty on a set of related vectors x
% and y with uncertainties xu and yu respectively.
if nargin<5
    number_std  = 3;
end
span = 1;
hasmatrix = false;
if ismatrix(x)&&ismatrix(y)
    [~,dim] = max(size(x));
    span    = min(size(x));
    xnumel = size(x,dim);
    ynumel = size(y,dim);
    hasmatrix = true;
else
    xnumel = numel(x);
    ynumel = numel(yu);
end

% coordinates
xcoords = linspace(min(x,[],'all'),max(x,[],'all'),xnumel);
ycoords = linspace(min(y,[],'all'),max(y,[],'all'),ynumel);

xpos = arrayfun(@(i) find(xcoords>=i,1,'first'),x);
ypos = arrayfun(@(i) find(ycoords>=i,1,'first'),y);

% 2d uncertainty
un = sqrt(xu.^2+yu.^2);

% make 1:3 sigma clouds
map0    = zeros(xnumel,ynumel,span);


for k=1:span
    for j=1:number_std
        map = zeros(xnumel,ynumel);
        for i=1:xnumel
            dilatefactor = round(mean([find(xcoords>un(i,k)*j,1,'first') find(ycoords>un(i,k)*j,1,'first')]));
            if isnan(dilatefactor)
                continue
            end
            xrange      = xpos(i,k)+round(dilatefactor/2).*[-1 1];
            yrange      = ypos(i,k)+round(dilatefactor/2).*[-1 1];

            xrange(xrange<1)=1;xrange(xrange>xnumel)=xnumel;
            yrange(yrange<1)=1;yrange(yrange>ynumel)=ynumel;
            map(xrange(1):xrange(2),yrange(1):yrange(2)) = 1;
        end
        map0(:,:,k) = map0(:,:,k)+map;
    end
end

if hasmatrix
    tmp0 = zeros(size(map));
    for j=fliplr(1:number_std)
        tmp = double(logical(sum(map0>j-1,3)));
        tmp0=tmp0+tmp;
    end
    map0 = tmp0;
end


map0 = map0';
if ~nargout
    clf
    contourf(xcoords,ycoords,map0)
    hold on
    plot(x,y,'k','LineWidth',2)
end

return