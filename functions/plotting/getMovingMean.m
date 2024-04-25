function movint = getMovingMean(files,loop,rotation,bin,trail,pos)
im0     = loadimage(files,1,rotation);
int     = zeros(size(im0,1),size(im0,2),bin);
for i=pos:pos+bin-1
   int(:,:,i) = loadimage(files,loop(i),rotation);
end
movint = movmean(int,trail,3);
end

