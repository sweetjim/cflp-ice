function Tb = Tb_heatflux(Sb,Sinf,Ti,Tinf,varargin)
%% Interface temperature (from heat-flux equation)
% Calculate the interface temperature from the interface and ambient
% salinity, ice and ambient temperature, and thickness of the ice 
% (default is H_x = 0.1 [m]).
% 
% Parameters: Sb, Sinf, Ti, Tinf, (Optional) H_x
% 
%%

H_x     = .1;               % thickness of ice [m]

c       = struct(...        % [J kg / degC]
    'i',[2009],...          % specific heat cap. of ice
    'w',[3974]);            % specific heat cap. of water

gamma   = struct(...        % []
    'T',[2.2e-2],...        % heat transfer coefficient
    'Ti',property('kappa','ice','T',Ti)/H_x,...  % heat in ice transfer coefficient 
    'S',[6.2e-4]);          % salinity transfer coefficient

rho     = struct(...        % [kg / m3]
    'i',[918.7],...         % ice density
    'w',density(Sinf,0,Tinf,0));  % water density

L       = 3.35e5;           % [J / kg]


A       = L/c.w*gamma.S/gamma.T;  % ~ 2.4
B       = rho.i./rho.w*c.i/c.w*gamma.Ti/gamma.T;

DeltaSinf   = Sb-Sinf;
DeltaS      = -Sb;


Tb = 1./(1+B)*(-DeltaSinf./DeltaS*A + Tinf+B*Ti);

end

