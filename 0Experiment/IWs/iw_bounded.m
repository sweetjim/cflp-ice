function out = iw_bounded
%%
i = 1;
j = 1;
aij = 1;
res = 5e2;
omega = 1;
x = linspace(0,1,res);
z = linspace(0,2,res);

n = 10;
m=20;
k0 = n*pi/max(x);
m0 = m*pi/max(z);
[X,Z,nn]= ndgrid(x,z,1:2);
t = 0;

psi = sin(nn.*X*k0).*sin(nn.*Z*m0).*exp(-complex(0,1)*omega*t);


imagesc(sum(psi,3))

%%
res = 5e2;
H = 1;
L = 2;
x = linspace(0,H,res);
z = linspace(0,L,res);

k0 = 5*pi/H;
m0 = 5*pi/L;
mode = 10;
[X,Z,nn,mm]= ndgrid(x,z,mode,mode);

N       = 1;
omega   = N*(1+((1*L)./(1*H)).^2).^(-1/2);
t       = pi*omega;
psi     = real(sin(nn.*X*k0).*sin(nn.*Z*m0).*exp(-complex(0,1)*omega*t));
imagesc(x,z,sum(psi,[3 4]))
end