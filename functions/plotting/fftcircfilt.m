function fout = fftcircfilt(in,radius)
%%
Ft      = fft2(in);
[nx,ny] = size(in);
[X,Y]   = meshgrid(-ny/2+1:ny/2,-nx/2+1:nx/2);
R2      = X.^2+Y.^2;
mask    = R2<radius^2;
fout    = ifft2((Ft.*mask));
fout    = fout-min(fout,[],'all');
fout    = fout/max(fout,[],'all');
fout    = real(fout);
% imagesc((fout))
end

