function [rho,S,T,z] = icesetup(Tamb,Samb,n)
%% Set spaces
% Generates structures for rho, S, and T where each structure contains the
% following variables:
%   'inf'   -   far-field   [double]
%   'b'     -   interface   [vec]
%   'i'     -   ice         [double]
%   'w'     -   far-field   [double]
% Calculates the interface conditions from a 4th-order polynomial
% thermal-driving-liquidus equation.
%%

addpath(genpath('functions'))

if nargin<3
    n = 100;
end
if Samb==0||Samb<0
    Samb = 1e-5;
end
x   = linspace(0,1.1,n);
z   = linspace(0,1,n);

S   = struct(...            % [g / kg]
    'inf',Samb,...            % ambient salinity
    'b',[],...              % interface salinity
    'i',0);                 % ice salinity

T       = struct(...        % [degC]
    'inf',Tamb,...             % ambient temperature
    'b',0,...              % interface temperature
    'i',-5);                % ice temperature

S.w     = linspace(S.i,S.inf,n);      % initial salinty profile
T.w     = T.inf;    % initial temperature profile

theta   = linspace(0,90,n);   % [deg]
g       = 9.8;              % [m / s2]
D       = 1e-2*0.146e-6;    % [m2 / s]

V       = @(theta,deltaS,rho) cos(deg2rad(theta)).*deltaS.*(rho.w./rho.i)*Gamma.S;

rho = struct(...            % [kg / m3]
    'w',density(S.w,T.w),...
    'b',[],...
    'i',property('ice','rho','T',T.i),...
    'ref',density(0,0));

L       = 3.35e5;           % [J / kg]
Gamma_S = 5.73e-2;          % [degC]
Gamma   = gammas('h',.1,'Ti',T.i);

try
    T.b     = Tb4thOrder(rho,S,T);
    S.b     = liquidus('S',T.b);
    rho.b   = density(S.b,T.b);
catch
    disp('Sort out Tb4thorder matrix dims')
end
end

