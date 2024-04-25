function hFound = pointer2object(figurehandle,object)
%POINTER2OBJECT determines if the pointer is positioned on the desired
%object
pos         = figurehandle.CurrentPoint;
hFound      = [];
tmp         = findobj(figurehandle);

for n=1:length(object)
    for j=1:length(tmp)
        try %#ok<TRYNC>
            if strcmp(tmp(j).Tag,object{n})
                idx = j;
                break
            end
        end
    end

    hTag = tmp(idx);
    hPos = getpixelposition(hTag,true);  % recursive to find absoulute position
    if within(pos(1),hPos(1).*[1 1]+[0 hPos(3)])&&within(pos(2),hPos(2).*[1 1]+[0 hPos(4)])
        hFound = hTag;
        break
    end
end
    function out = within(x,xrange,type)
        %% Within
        % Returns a logical argument that depends on whether 'x' falls within the
        % range of 'xrange'. Can also specify inclusive (DEFAULT i.e., <=) or
        % exclusive bounds (i.e., <)
        %
        % Examples:
        % >> within(0,[-1 1])
        % >> ans =
        % >>    logical
        % >>        1
        % --------------------
        % >> within(5.1,0:5)
        % >> ans =
        % >>    logical
        % >>        0
        % --------------------
        % >> within(1,1:2,'exclusive')
        % >> ans =
        % >>    logical
        % >>        0
        %%
        x1 = xrange(1);
        x2 = xrange(end);
        if nargin<3
            type = 'inclusive';
        end
        switch type
            case 'inclusive'
                out = (x>=x1&x<=x2);
                if diff(xrange)<0
                    out = (x<=x1&x>=x2);
                end
            case 'exclusive'
                out = (x>x1&x<x2);
                if diff(xrange)<0
                    out = (x<x1&x>x2);
                end
        end


    end
end
