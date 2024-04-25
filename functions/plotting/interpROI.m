function varargout = interpROI(roi1,roi2,steps,currentstep,type)
%INTERPROI interpolates between two ROIs over an interval t.
[x1,y1] = checkROI(roi1);

if nargin<=2
    verts = 1e2;
    if nargin==2
        if numel(roi2)==1
            verts = roi2;
        elseif numel(roi2)>1
            error('Cannot interpolate single ROI with specified vertices')
        end
    end
    [x1,y1] = interpolation(x1,y1,verts,'linear');
    varargout{1} = [x1;y1]';
    return
end
[x2,y2] = checkROI(roi2);
if nargin<5
    type = 'polyshape';
end

if numel(x1)==numel(x2)
    verts = numel(x1);
else
    verts = max([numel(x1) numel(x2)]);
end

% p = [0 cumsum(hypot(diff(x1),diff(y1)))'];
% p1 = linspace(p(1),p(end),verts);
% x1 = interp1(p,x1,p1,'previous');
% y1 = interp1(p,y1,p1,'previous');
% 
% p = [0 cumsum(hypot(diff(x2),diff(y2)))'];
% p1 = linspace(p(1),p(end),verts);
% x2 = interp1(p,x2,p1,'previous');
% y2 = interp1(p,y2,p1,'previous');

[x1,y1] = interpolation(x1,y1,verts);
[x2,y2] = interpolation(x2,y2,verts);

if currentstep<0;currentstep=0;elseif currentstep>steps;currentstep=steps;end
switch currentstep
    case 0
        P = polyshape(x1,y1,'Simplify',true);
        if strcmpi(type,'vertices')
            varargout{1} = P.Vertices;
            return
        end
        varargout{1} = P;
        return
    case steps
        warning off
        P = polyshape(x2,y2,'Simplify',true);
        warning on
        if strcmpi(type,'vertices')
            varargout{1} = P.Vertices;
            return
        end
        varargout{1} = P;
        return
    case 'all'
        currentstep = [];
end
% create surface
X = [x1;x2];
Y = [y1;y2];

if nargin<3
    steps = 2;
end
if nargin<4
    currentstep = step/2;
end
Z = [x1*0;x2*0+steps];%imresize(,[t numel(x1)],'bilinear');

fig     = figure('Visible','off');
if isempty(currentstep)
    currentstep = 1:steps;
end
S       = contour3(X,Y,Z,'LevelList',currentstep,'linewidth',3);
[x,y]   = C2xyz(S);

delete(fig)
% delete(findall(0,'Type','Figure','Visible','off'))

if ~nargout
    contour3(X,Y,Z,'LevelList',currentstep,'linewidth',3)
    %     surface(X,Y,Z,'edgecolor','none');
    light
    view(0,90)
    clear X Y
elseif nargout==1
    if numel(currentstep)>1
        x=reshape(x,steps,[]);
        y=reshape(y,steps,[]);
        tmp = images.roi.Polygon;
    else
        warning off all
        P = polyshape(x,y,'Simplify',true,'KeepCollinearPoints',false);
        warning on
    end
    if P.overlaps
        [x,y]       = P.boundary;
        [~,sorting] = sort(y);
        vertices    = fillmissing([x(sorting) y(sorting)],'linear');
        if strcmpi(type,'vertices')
            varargout{1} = vertices;
            return
        end
        P.Vertices  = vertices;
    end
    %     P.Vertices(isnan(P.Vertices))=[];
    varargout{1} = P;
elseif nargout==2
    varargout{1} = X;
    varargout{2} = Y;
end

%% Nested functions
    function [x,y] = interpolation(x,y,verts,method)
        if nargin<4
            method = 'previous';
        end
        p = [0 cumsum(hypot(diff(x),diff(y)))'];
        p1 = linspace(p(1),p(end),verts);
        x = interp1(p,x,p1,method);
        y = interp1(p,y,p1,method);
    end
    function [x,y] = checkROI(roi)
        status = cellfun(@(x) isa(roi,['images.roi.' x]),{'Rectangle','Polygon'});
        if any(status)
            switch class(roi)
                case 'images.roi.Rectangle'
                    x = (roi.Position(1).*[1 1]+[0 roi.Position(3)])';
                    y = (roi.Position(3).*[1 1]+[0 roi.Position(4)])';
                otherwise
                    x = roi.Position(:,1);
                    y = roi.Position(:,2);
            end
            return
        end
        x = roi(:,1);
        y = roi(:,2);
        [y,sorted] = unique(y);
        x = x(sorted);
    end
    function pshape = createROI(X,Y)
        pgon        = polyshape(X,Y,'Simplify',true);
        [~,order]   = sort(pgon.Vertices(:,2));
        pos         = pgon.Vertices(order,:);
        fig         = figure('Visible','off');
        pshape      = drawpolygon('Position',pos);
        delete(fig)
    end
end