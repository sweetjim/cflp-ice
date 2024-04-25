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
if isempty(xrange)
    out = false;
    return
end
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

