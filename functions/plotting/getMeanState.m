function int = getMeanState(files,loop,rotation)
im0     = loadimage(files,1,rotation);
int     = zeros(size(im0));
for i=loop
   im  =  loadimage(files,i,rotation); 
   int = int+double(im);
   displayProgress('Integrating',i,loop(1),loop(end))
end
int = int/length(loop);
end

