% function lagrangian_of(images)
% %LAGRANGIAN_OF attempts to create a simplistic Lagrangian reference frame
% %(from a Euclidean reference frame) by shifting the image time-series by
% %pixel quantities per frame and perform an Optical Flow analysis.
% 
% [fig,ax]=previewImages(images,'gray',[-10 10]);
% 
% 
% end
%%
shift = -0;
k = 1;

clear vf
K = size(images,3)-1;
for k=1:K
im1 = images(:,:,k);
im2 = images(:,:,k+1);%fraccircshift(images(:,:,k+1),shift,2);

im1(~mask)=0;
im2(~mask)=0;

[xa,ya]=boundingbox(polyshape(roi(1).Position)); 
crop = [xa(1) ya(1) diff(xa) diff(ya)];
mask=createMask(roi(1));
im1 = imcrop(im1,crop);
im2 = imcrop(im2,crop);
vf = velocity_field(im1,im2,'fast',25,200);
if k==1
    U = zeros([size(vf.u) K]);
    V = U;
end
U(:,:,k) = vf.u;
V(:,:,k) = vf.v;
displayProgress('Iterating',k,1,K,'delete')
pause(1e-4)
end
tmp = imcrop(mask.*(1+images(:,:,1)*0),crop);

U = mean(U,3);
V = mean(V,3);

%%

[xa,ya]=boundingbox(polyshape(roi(1).Position)); 
crop = [xa(1) ya(1) diff(xa) diff(ya)];
mask=createMask(roi(1));
fac = std2(hypot(u,v));
u = imcrop(U.*mask,crop)/fac;
v = imcrop(V.*mask,crop)/fac;


figure(6)
tiledlayout(3,1)
nexttile
imagesc(imcrop(im1,crop)-imcrop(im2,crop))
colormap gray
streamslice(u,v,5)
nexttile
imagesc(u)
clim = caxis;
cmocean('balance','pivot',0)
hold all
streamslice(u,-v,1)
contour(flowfun(u,v),'k','LineWidth',2)
caxis(clim)
nexttile
imagesc(v)
clim = caxis;
cmocean('balance','pivot',0)
hold all
contour(flowfun(u,v),'k')
caxis(clim)
