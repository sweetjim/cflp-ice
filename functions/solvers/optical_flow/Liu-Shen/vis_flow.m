function [Ox,Oy] = vis_flow (ux,uy,varargin)% gx, offset, mag, col,ds,ax)
% Quiver, with subsampling
%
% Usage:
% gx      = 25;     % samples
% offset  = 0;      % arrow offset
% mag     = 1;      % arrow magnitude
% col     = 'b';    % arrow color
% ds      = 1;      
% ax      = [];
% X       = [];
% Y       = [];
%% Inputs
gx      = 25;
offset  = 0;
mag     = 1;
col     = 'b';
ds      = 1;
ax      = [];
X       = [];
Y       = [];
parseInput(varargin)
if isempty(ax)
    ax = gca;
end
%% Quiver
[sy sx] = size(ux);
if (gx==0)
    jmp = 1;
else
    jmp = floor(sx/gx);
    jmp = jmp + (jmp==0);
end

indx = (offset+1):jmp:sx;
c = 1;
CX = [];
CY = [];
for j=(1+offset):jmp:sy
    Vx(c,:) = ux(j,indx);
    Vy(c,:) = uy(j,indx);
    CX(c,:) = indx;
    %CY(c,:) = ones(size(indx)).*(sy-j+1);
    CY(c,:) = ones(size(indx)).*j;
    c = c+1;
end

if (isnan(Vx(1,1)))
    Vx(1,1) = 1;
    Vy(1,1) = 0;
    CX(1,1) = 1;
    CY(1,1) = 1;
end

%%
M = ~isnan(Vx) & ~isnan(Vy);
if sign(mean(diff(Y)))<1
   CY = flipud(CY);
end
xx = CX(M)/ds;
yy = CY(M)/ds;
if all([~isempty(X) ~isempty(Y)])
    xx = rescale(xx,min(X),max(X));
    yy = rescale(yy,min(Y),max(Y));
end



H = quiver (ax,xx,yy, Vx(M), Vy(M), mag);
s = size(ux);
% axis (ax,[0 s(2)/ds 0 s(1)/ds]);
set (H, 'Color', col);

%% Output
switch nargout
    case 0
        clear Ox;
        clear Oy;
    case 1
        Ox = H;
    otherwise
        Ox = Vx;
        Oy = Vy;
end
%% Function
    function parseInput(varargin)
    m = 1;
    items = varargin{:};
    for k=1:length(items)
        switch items{m}
            case 'ax'
                ax = namevalue;
            case {'points','gx'}
                gx = namevalue;
            case {'mag','magnitude'}
                mag = namevalue;
            case {'offset'}
                offset = namevalue;
            case {'color','col','c'}
                col = namevalue;
            case {'ds'}
                ds = namevalue;
            case {'X','x'}
                X = namevalue;
            case {'Z','z','Y','y'}
                Y = namevalue;
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

