function pshape = shadedErrorBarPoly(x,y,varargin)
%SHADEDERRORBARPOLY generates a shaded error bar from a patch object.
% Arguments:
%   x               - [1xn or nx1 vector] independent variable
%   y               - [2xn or nx2 vector] dependent variable error bounds
%   'auto'          - [char] automatically determine errorbounds from raw
%                     data
%   'direction'     - [char] direction of the independent variable: 'x' (default) or 'y'
%   'color'         - [char, RBG triplet] color of the patch (default is 'r')
%   'facealpha'     - [double] alpha transparancy of the patch face (default is 0.1)
%   'edgealpha'     - [double] alpha transparancy of the patch edge (default is 0.1)
%   'ax'            - [axes] axes to plot (default is gca)
direction   = 'x';
color       = 'r';
facealpha   = .1;
edgealpha   = .1;
ax          = [];
auto        = false;
parseInput(varargin)
if isempty(ax)
    ax = gca;
end

if auto
    dir = find(any(size(x)'==size(y),2));
    y = mean(y,dir,'omitnan')+[-1 1].*std(y,[],dir,'omitnan'); 
end


[d1,d2]=size(y);
if d1>=d2
    ypos=y(:,1);
    yneg=y(:,2);
else
    ypos=y(1,:);
    yneg=y(2,:);
end


if any(isnan(y))
    omit = max([find(~isnan(ypos),1,'first') find(~isnan(yneg),1,'first')]);
    yneg(1:omit) = [];
    ypos(1:omit) = [];
    x(1:omit) = [];
    omit = min([find(~isnan(ypos),1,'last') find(~isnan(yneg),1,'last')]);
    yneg(omit:end) = [];
    ypos(omit:end) = [];
    x(omit:end) = [];
end

if isrow(x);x = x';end
if isrow(yneg);yneg=yneg';end
if isrow(ypos);ypos=ypos';end
warning off

omit = ~any(isnan(yneg),2);
x=x(omit);
yneg = yneg(omit);
ypos = ypos(omit);
switch direction
    case 'x'
        pshape = patch(ax,[x;flipud(x)],[yneg;flipud(ypos)],color);
        %pshape = polyshape([x;flipud(x)],[yneg;flipud(ypos)]);
    case 'y'
        pshape = patch(ax,[yneg;flipud(ypos)],[x;flipud(x)],color);
        %pshape = polyshape([yneg;flipud(ypos)],[x;flipud(x)]);
end
set(pshape,'FaceAlpha',facealpha,'EdgeAlpha',edgealpha)
warning on
% holding = ishold(ax);
% if ~holding;hold(ax,'on');end
% plot(pshape)
% if ~holding;hold(ax,'off');end

if ~nargout
    clear pshape
end
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case 'ax'
                    ax = namevalue;
                case {'k','b','r','g','c','m','w'}
                    color = items{m};
                case 'color'
                    color = namevalue;
                case {'direction'}
                    direction = namevalue;
                case 'edgealpha'
                    edgealpha = namevalue;
                case 'facealpha'
                    facealpha = namevalue;
                case {'auto','automatic'}
                    auto = true;
                otherwise
                    if isa(items{m},'double')
                        if (numel(items{m})==3)
                            color = items{m};
                        end
                    end
            end
            m = m+1;
            if m>length(items);break;end
        end
        function out = namevalue
            out = items{m+1};
            m   = m+1;
        end
    end
end

