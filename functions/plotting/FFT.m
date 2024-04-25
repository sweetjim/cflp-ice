function [Fhat,F,k] = FFT(f,var)
n       = length(var);
F       = fft(f,n);
freq    = 1/(mean(diff(var))*n)*(-n/2:n/2-1);
L       = 1:floor(n/2+1);
PDS     = F.*conj(F)/n;

k       = freq;%(L);
Fhat    = PDS;%L);
end

