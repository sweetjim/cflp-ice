function [objs,sorting] = sortobjs(objs,heirarchy)
%SORTOBJS Summary of this function goes here
if nargin<2
    heirarchy = 'simple';
end

switch heirarchy
    case 'simple'
        classarr = arrayfun(@(x) class(x),objs,'UniformOutput',false);
        if any(contains(classarr,'Image'))
            [~,sorting]     = sort(classarr(~contains(classarr,'Image')));
            sorting(end+1)  = find(contains(classarr,'Image'));
        end
        objs = objs(sorting);
    otherwise
        % to add
end
Ia = sorting;
Ib = 1:numel(objs);
end

