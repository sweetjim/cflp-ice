function compare_of_implementations
ims = imageDatastore("functions\solvers\optical_flow\", ...
    "IncludeSubfolders",true, ...
    'FileExtensions',{'.tif'},...
    'LabelSource','foldernames');

pairs = {'oval','twin','wall','vortex'};
idx = find(contains(ims.Files,pairs{2},'IgnoreCase',1));

im1 = readimage(ims,idx(end-1));
im2 = readimage(ims,idx(end));
clc

figure(1)
clf
imagesc(im1)%,colormap gray,out2latex('renders\process\examples\case2_a','r',600)
imagesc(im2)%,out2latex('renders\process\examples\case2_b','r',600)
imagesc(imfuse(im1,im2))

clf
out=velocity_field(double(im1),double(im2),'slow',20,2e2,2000);
u = out.u;
v = out.v;
psi = flowfun(u,v,'+');
clf
imagesc(im1)
hold on,vis_flow(gather(u),gather(v),'gx',200)
streamslice(u,v,4)
return
%%
% imagesc(double(im1)-double(im2))
%% Calculation
% parameters
method = 'fast';
clc
lambda = 1;
alpha = [1 50 100 200 400];
its   = [1 50 100 200 500];


vfcalc = @(method,alpha,its,lambda) velocity_field(double(im1),double(im2),method,alpha,its,lambda);
vf = vfcalc(method,1,1,1);
%
vf=preallocateStruct(vf,[numel(alpha),numel(its)]);
clear n
k = 1;
c_a = numel(alpha);
c_i = numel(its);
q = parallel.pool.DataQueue;
parfevalOnAll(@clear, 0,'all');
parfevalOnAll(@warning, 0,'off','all');
afterEach(q,@displayProgressparfor)
displayProgressparfor([],c_a*c_i)
parfor i=1:c_a
    for j=1:c_i
        out = vfcalc(method,alpha(i),its(j),lambda);
        vf(i,j) = out;
        %         k=k+1;
        %         displayProgress('Calculating',k,1,numel(alpha)*numel(its)+1,'delete')
        send(q,[])
    end
end
% vf = reshape(vf,c_a*c_i,1);
%% Plotting
figure(2)
tile = tiledlayout(5,5,'TileSpacing','compact','Padding','compact');

for i=1:c_a
    for j=1:c_i
        nexttile;
        %     imagesc(gather(vf(i).u./(2*std2(vf(i).u))))
        imagesc(gather(vf(i,j).v./(2*std(vf(i,j).v,[],'all','omitnan'))))
        %         imagesc(vf.im2)
        %         caxis([0 1])
        %         colormap gray
        %         hold on
        %         vis_flow(vf.u,vf.v,'gx',80,'mag',2)
        %         streamslice(vf.u,vf.v)
    end
    displayProgress('Building',i,1,numel(vf))
end
addlabels('ax',tile,'array','x',alpha,'y',its)
addlabels('ax',tile,'x','Diffusion (\alpha)','y','Iterations')
set(tile.Children,'XTick',[],'YTick',[])
arrayfun(@(x) cmocean('balance','pivot',0,x),tile.Children)
%%
figure(1)
hold on
delete(findall(gcf,'Type','Quiver'))
delete(findall(gcf,'Type','Line'))
vis_flow(vf(end,end).u,vf(end,end).v,'gx',100,'mag',2)