function loopImages(files,x,z,loop,rotation)
%%
dt = 0.3;
clf

V = VideoWriter(folder,'mpeg-4');
open(V);
set(gcf,'color','w')

for i=loop
    %%
    imin = loadimage(files,i,rotation);
    imagesc(x,z,gradient(imin));
    cmocean('gray')
    caxis([-10 10])
    addlabels(gca,sprintf('t = %.2f (s)',dt*i)) %getdt(files,i,'min')
    axis xy
    writeVideo(V,getframe(gcf))
    pause(.01)
end
close(V)
end

