function rho = density(S,T,type)
%% Density
% Calculates the quasi-linear density of water from haline and thermal
% concentrations (in [g/kg] and [degC]) referenced from 
% (S,T) = (0 psu,0 degC) at 0 bar.
%
% Parameters: REQUIRED
%   S   -   [1xn vector, double] salinity      (psu/perthousand)
%   T   -   [1xm vector, double] temperature   (degC)
% Parameters: OPTIONAL
%   type -  [integer,string] 
%           If S(1xn) and T(1xm) are both vectors
%           Use '2d' or 2 to output nxm density (DEFAULT)
%           Use '1d' or 1 to output 1xn density (if n==m)
%% TEOS
if nargin<3
    type = '2d';
end


areVectors = [isvector(S) isvector(T)];

if all(areVectors)
    areRows = [isrow(S) isrow(T)];
    areCols = [iscolumn(S) iscolumn(T)];

    if and(any(areRows),any(areCols))
        [s,t]   = meshgrid(S,T);
        rho     = sw_dens0(s,t);
        return
    end
    if numel(S)==numel(T)
        rho = sw_dens(S,T,0);
        return
    end
    if numel(T)==1
        rho = arrayfun(@(x) sw_dens(x,T,0),S);
        return
    end
    if numel(S)==1
        rho = arrayfun(@(x) sw_dens(S,x,0),T);
        return
    end
    [s,t]   = meshgrid(S,T);
    rho     = sw_dens(s,t,0);
else
    rho = sw_dens(S,T,0);
    switch type
        case {'1d',1}
        rho = diag(rho)';
    end
end
return

%% Ruddick
rho0 = .99708; % g/cm3
a0 = -2.539e-4;
a1 = -4.968e-6;
a2 = 2.7e-8;
a3 = .69976;
a4 = .14042;
a5 = .337;
b0 = 1.6803e-3;
b1 = -3.551e-5;
b2 = 3.52e-7;
c0 = 2.714e-3;
c1 = -8.11e-5;
c2 = 9e-7;

t0      = @(T) T-25;
s0      = @(S) S*1e-3;
drho    = @(S,T) ...
            (a0*t0(T)+a1*t0(T).^2+a2*t0(T).^3)+...
            (a3*s0(S)+a4*s0(S).^2+a5*s0(S).^3)+...
        -   (b0*t0(T)+b1*t0(T).^2+b2*t0(T).^3).*s0(S)+...
            (c0*t0(T)+c1*t0(T).^2+c2*t0(T).^3).*s0(S).^2;
rho     = 1e3*(drho(S,T)+rho0);
%% Roquet et al. (2015)
% Defining a simplified yet “realistic” equation of state for seawater
rho = @(T,S) -.011/2*(T-4+.25*S).^2+.77*S+999.8;
end

