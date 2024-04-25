function [eta,rho0,rho] = calculateLayers(Tinf,Sinf,Z,rho,discretize,rho0_temp)
%CALCULATELAYERS 
% From Huppert & Turner (1980). 
%   eta = 0.66pm0.06*(density(0,Sinf)-density(Tinf,Sinf))/dzRhoinf
% From Hewitt (2020)
%   eta = 2.85*St/(E+St)*c/sin(phi)/dzRhoInf*DTa0*Drhoief/Ltilde
% From Magorrian & Wells (2016)
%   eta = Drho_u/drzRhoinf
%   where Drho_u = M/(E+M)*Drhoi0ef
%       assume Drhoi0ef=0.024*rho0
%       let     M = mdot/U
%   eta = 0.024*rho0/(dzRhoinf*(1+EU/mdot))
%% Input parsing
if nargin<4
    rho     = density(Sinf,Tinf);
end
if isempty(rho)
    rho     = density(Sinf,Tinf);
end
if nargin<5
    discretize = false;
end
%% HT80
Tw      = 0; % 
L       = property('L');
c       = property('water','cp','T',Tinf);
Sw      = L*Sinf/(L+c*(Tinf-Tw));

if isempty(rho)
    rho     = density(Sinf,Tinf);
end
if nargin<6
    rho0_temp = 0;
end
rho0    = density(ones(size(Sinf)).*Sw,rho0_temp);

if mean(diff(Z))<0
    rho     = fliplr(rho);
    rho0    = fliplr(rho0);
end
eta     = abs(.66*(1+.06*[-1 1]').*(rho0-rho)/mean(gradient(rho,Z)));
%% H20
phi     = pi/2;           % slope angle [rad]
Cd      = 2.5e-3;           % drag coefficient              | !!UNCERTAINTY!!
E0      = 3.6e-2;           % entrainment coefficient       |
St      = 1.1e-3;           % heat transfer coeff           |
Stm     = 3.1e-5;           % salt transfer coeff           |

E       = E0*sin(phi);      % entrainment rate param
alpha   = 3.87e-5;          % thermal expansion coeff
beta    = 7.86e-4;          % haline contraction coeff
cS      = 2.009e3;          % heat capacity ice
cL      = 3.974e3;          % heat capacity liquid
L       = 3.35e5;           % latent heat
Gamma   = 5.73e-2;          % liquidus freezing offset (linear approx.)
g       = 9.81;             % gravity
rho0    = 999.8;            % 0degC freshwater reference density
nu      = 1e-6;             % kinematic viscosity
kappaS  = property('water','kappaS','T',0); % molecular solutal diffusivity

Drhoief = 0.024*rho;
Ti      = 1e-3;
Ltilde  = L+cS*(liquidus('T',0)-Ti);
DTa0    = Tinf(1);
dZrhoA  = abs(mean(gradient(rho,Z)));
ETA     = 2.85*St/(E+St)*cL/sin(phi)/dZrhoA*Tinf*Drhoief/Ltilde;
%% Discretisation / layer steps
if discretize
    %%
    ETA     = abs(mean(eta));
    [Z,map] = sort(Z);
    ETA     = ETA(map);
    %%
    newidx  = find(Z==min(Z));
    newz    = min(Z);
    i = 1;
    clear idxs
    
    while newz<max(Z)
        newz    = ETA(newidx)+newz;
        newidx  = find(Z>newz,1,'first');
        if isempty(newidx)
            break
        end
        idxs(i) = newidx;
        i = i+1;
    end
    idxs = [1 idxs];
    
    vals = diff(Z(idxs));
    for i=1:numel(idxs)-1
        ETA(idxs(i):idxs(i+1))=vals(i);
    end
    eta = ETA(map);
end

