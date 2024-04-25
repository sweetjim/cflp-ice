function [Fk,k,l] = ffts2(array2d,vec1,vec2)
%% FFTS2
% Two-dimensional FFT with dimensional vector arguments. Returns the
% transformed matrix and associated wavevectors.
% Parameters: 
%   array2d : 2-dimensional array to perform FFT2 on
%   vec1    : first dimension vector
%   vec2    : second dimension vector
% Outputs:
%   Fk      : Fourier transform of array2d
%   k       : wavevector of first vector
%   l       : wavevector of second vector
%%
if nargin<2
    vec1 = 1:size(array2d,1);
end
if nargin<3
    vec2 = 1:size(array2d,2);
end

[Fk,k] = ffts(array2d,mean(diff(vec1)),1,numel(vec1),1);
[Fk,l] = ffts(Fk,mean(diff(vec2)),2,numel(vec2),1);

% [F_time,omega_space] = ffts(anom,dt,2,length(t),1);                     % FFT in time
% [F_space,k_space]    = ffts(F_time,dx,1,size(anom,1),1);                % FFT in space
% F_space              = fftshift(fftshift(F_space,1),2);                 % FFT centering

k   = (k-mean(k)-(k(2)-k(1))/2);
l   = -(l-mean(l)-(l(2)-l(1))/2);

end

