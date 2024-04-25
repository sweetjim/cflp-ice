function ddc_plume(varargin)
%DDC_PLUME Summary of this function goes here
%   Detailed explanation goes here
q       = -50/4000; % heat flux (dT/cp)
Ti      = 0;        % wall temperature
kappa   = 0.6/4000; % W / m2
m       = 2*pi/0.3; % lengthscale of DDC layers (m)
i       = complex(0,1);
wbar    = 0.3;      % vert. velocity amplitude (m/s)
D       = 0.04;     % plume width (m)
k       = sqrt(i*m*wbar./kappa-m^2);
z       = linspace(-2*pi/m*5,0,100);
x       = linspace(0,D,50);
[X,Z]   = ndgrid(x,z);
T       = Ti-q/kappa*(X+imag(-i./k./cos(k*D)*sin(k*X).*exp(i*m*Z)));
dTdx    = -q/kappa*(1+imag(-i./cos(k*D)*cos(k*X).*exp(i*m*Z)));
subplot(1,3,[1 2]);
pcolor(x,z,T');shading flat;
colorbar
title('Temp w.r.t ice face at x=0');
subplot(1,3,3);
plot(dTdx(1,:),z,dTdx(end,:),z);
title('dT/dx at edges');
end

