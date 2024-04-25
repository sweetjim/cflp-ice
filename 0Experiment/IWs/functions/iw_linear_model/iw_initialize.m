function out = iw_initialize(H,z,x,h,N,u,ut,omega,f,nu,n,alpha0)
%% iw_initialize 
% Models an bottom-forced internal wave using Bell's (1975) linearized
% equations, Lewellyn & Young's (2001) finite depth reformulation, and
% Shakespeare's (2021) biharmonic viscous dissipative regime. 
% 
% Parameters:
%   H: domain depth / height
%   z: vertical vector
%   x: horizontal vector
%   h: z(x) profile of forcing feature (i.e., topography)
%   N: stratification (rad/s)
%   u: linear velocity (m/s)
%   ut: oscillating velocity magnitude (m/s)
%   omega: oscillation frequency (rad/s)
%   f: Coriolis frequency (rad/s)
%   nu: viscosity (m2/s)
%   n: number of harmonics to compute
%   
%%
if (N/omega<max(n)+1)&&u==0
    n = 0:N/omega;
end

% complex numbers
i=complex(0,1);

% topo transform
[hhat,k]    = ffts(h,x(2)-x(1),2, length(x));

% grids
[K,Z,nn]    = ndgrid(k,z,n);
[Hhat,~,~]  = ndgrid(hhat,z,n);
% nth-horizontal Laplacian dissipation
alpha0      = nu.*K.^alpha0;

% downstream waves i.e. -n*omega+k*u  (n<=0)
OMEGA0      = -nn*omega+K*u;
% the topo boundary condition
w0          = Hhat.*i.*OMEGA0.*besselj(-nn,(K*ut)/omega);
% the vertical wavenumber
m0          = real(abs(K).*sqrt((N.^2-(OMEGA0).^2)./(OMEGA0.^2-f.^2)).*sign(OMEGA0));
% make zero evanescent solutions
m0((N.^2-OMEGA0.^2)./(OMEGA0.^2-f.^2)<0)=0;
% the decay rate
gamma0      = real(alpha0.*m0.*OMEGA0.*(2*N^2-OMEGA0.^2-f^2)./(N^2-OMEGA0.^2)./(f^2-OMEGA0.^2)/2);


% upstream waves  i.e. n*omega+k*u  (n>=0)
OMEGA1      = nn*omega+K*u;
% the topo boundary condition
w1          = Hhat.*i.*OMEGA1.*besselj(nn,(K*ut)/omega);
% the vertical wavenumber
m1          = real(abs(K).*sqrt((N.^2-OMEGA1.^2)./(OMEGA1.^2-f.^2)).*sign(OMEGA1));
% make zero evanescent solutions
m0((N.^2-OMEGA1.^2)./(OMEGA1.^2-f.^2)<0)=0;
% the decay rate
gamma1      = real(alpha0.*m1.*OMEGA1.*(2*N^2-OMEGA1.^2-f^2)./(N^2-OMEGA1.^2)./(f^2-OMEGA1.^2)/2);

% depth dependence of flow (what=1 at Z=-H and what=0 at Z=0)
what0       = (exp(-i*(m0.*Z)-gamma0.*(Z-H))-exp(i*(m0.*Z)+gamma0.*(Z+H)))./(exp(i*m0.*H+2*gamma0.*H)-exp(-i*m0.*H));
what1       = (exp(-i*(m1.*Z)-gamma1.*(Z-H))-exp(i*(m1.*Z)+gamma1.*(Z+H)))./(exp(i*m1.*H+2*gamma1.*H)-exp(-i*m1.*H));
what1(isnan(what1)) = 0;
what0(isnan(what0)) = 0;

% sum the upstream (n<=0) and  downstream (n=>0) to give net
what        = (w0.*what0+w1.*what1);
% correct n=0 harmonic component - it is a duplicated when summing n>=0 and n<=0
what(nn==0) = w0(nn==0).*what0(nn==0);
what(isnan(what))   = 0;
% calculate IFFT
w           = squeeze(sum(iffts(what,length(x),1,length(x)),3));

%% Buoyancy

b0          = w0.*N^2./(i.*(OMEGA0));          % Buoyancy (continuity eq in fourier space)
b1          = w1.*N^2./(i.*(OMEGA1));
bhat        = (b0.*what0+b1.*what1);
% bhat(isnan(bhat))=0;

bhat(nn==0) = bhat(nn==0)/2; % remove duplicate n=zeros

% bhat(nn==0)=w0(nn==0).*bhat0(nn==0);

bhat(isnan(bhat))       = 0;
bhat(isinf(abs(bhat)))  = 0;

b_nn    = iffts(bhat,length(x),1,length(x));
b       = squeeze(sum(b_nn,3));
%% Output
out = struct(...
    'k',K,...
    'nn',nn,...
    'w0',w0,...
    'w1',w1,...
    'what0',what0,...
    'what1',what1,...
    'OMEGA0',OMEGA0,...
    'OMEGA1',OMEGA1,...
    'w',w,...
    'u',-cumsum(gradient(w')'),...%-w*tan(atan(omega/N)),...
    'b',b,...
    'bnn',b_nn);
end