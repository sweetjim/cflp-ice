function fft1D(in,var,mask_range)
[Fk,k]=ffts(in,mean(diff(var)),1,length(var));

if isempty(mask_range)
    mask = 1;
else
    mask  = double(var>mask_range(1)&var<mask_range(2));
end
fin = iffts(Fk.*mask,length(k),1,length(k));
tiledlayout(2,1)

nexttile
plot(k,log10(abs(Fk).^2))

nexttile
plot(var,in,'k',var,fin,'--r')
end

