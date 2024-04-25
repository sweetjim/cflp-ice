function [map,crop,warpfield] = warpCoordinates(imref,pgon,nomap,warpfield)
%WARPCOORDINATES Summary of this function goes here
%   Detailed explanation goes here
%% refine pos
if nargin<3
    nomap = false;
end
ims = imref;
if size(imref,3)>1
    imref   = imref(:,:,1);
end
posref     = pgon.Position;
warning off
% posref  = interpROI(pos,100);
% posref  = smooth2a(posref,5,0);
% pgon.Position = posref;
[xa,ya] = boundingbox(polyshape(posref));
crop    = [xa(1) ya(1) diff(xa) diff(ya)];
%% get mask (field 1)
mask    = createMask(pgon);
map     = imgaussfilt(bwdist(~mask),10).*mask;
if nomap
    return
end
imref   = imcrop(imref,crop);
MASK    = imcrop(mask,crop);
field0  = map;
U       = meshgrid(1:size(imref,2),1:size(imref,1));

%% warping field (v1)
% find coordinates of "streamlines"
hasField = true;
if nargin<4
    hasField = false;
end

if ~hasField
    fprintf('Calculating warping field\n')
    steps       = linspace(0,1,1e2);
    [x1,y1]     = C2xyz(contour(rescaleMax(field0),'LevelList',steps));
    idx         = repmat(find(x1==max(x1))',1,2);
    idx(:,2)    = circshift(idx(:,2),-1)-1;
    idx(end,:)  = [];
    flag        = logical(idx(:,1)-idx(:,2));
    idx(~flag,:)=[];

    dim = size(idx,1);
    X1 = cell(1,dim);
    Y1 = X1;

    for j=1:dim
        X1{j} = x1(idx(j,1):idx(j,2));
        Y1{j} = y1(idx(j,1):idx(j,2));
    end
    pts = max(cellfun(@numel,X1));
    X1  = cellfun(@(x) imresize(x,[1 pts]),X1,'UniformOutput',false);
    Y1  = cellfun(@(x) imresize(x,[1 pts]),Y1,'UniformOutput',false);
    x1  = reshape(cell2mat(X1),pts,[]);
    y1  = reshape(cell2mat(Y1),pts,[]);
    warpfield = cat(3,x1,y1);
else
    x1 = warpfield(:,:,1);
    y1 = warpfield(:,:,2);
end

% ims = preloaded-movmean(preloaded,10,3);
imdim = size(ims,3);
TMP   = zeros(pts,dim,imdim);

fprintf('Applying warping field\n')
for k=1:imdim
    TMP(:,:,k)=warpingFieldLoop(ims(:,:,k),crop,pts,dim,warpfield);
end

%% warping field (v2)
% map     = imregdemons(field0,U, ...
%     'DisplayWaitbar',false, ...
%     'AccumulatedFieldSmoothing',1);

%% Nested functions
    function out = warpingFieldLoop(ims,crop,pts,dim,wf)
        im  = imcrop(ims,crop);
        out = zeros(pts,dim);
        %%
        for jj=1:pts
            for ii=1:dim
                ypos = floor(wf(jj,ii,2));
                xpos = floor(wf(jj,ii,1));
                xpos(xpos<1)=1;
                ypos(ypos<1)=1;
%                 displayProgress('Warping',jj+ii-1,1,pts*dim)
                
                try
                    out(jj,ii)=im(ypos,xpos);
                catch
                    continue
                end
            end
        end
    
    end
end

