function out = interface(varargin)
%% Interface, Liquidus relation
% Calculates the linear density of water due to haline and thermal
% concentrations and the liquidus line.
% 
% Parameters: (Optional) z
% 
% Outputs a structure containing the density field (as a function of T and
% S, assuming T0~999 [kg / m3] and S0=0 [g/kg]), the density at the
% liquidus curve (at the presribed depth; default is z = 0 [m]), S and T
% arrays and vectors, as well as the liquidus freezing temperature at the
% prescribed depth (Tf).

%% Defaults
if nargin==1
    z = varargin{1};
else
    z = 0:1e3:8e3;
end
%% Script
n = 1e3;
s = linspace(0,300,n);
t = linspace(-30,100,n);


[S,T]   = meshgrid(s,t);
Tf      = liquidus('T',s);
liquid  = T>Tf(1,:);
rho     = sw_dens(S,T,0);
% rho(rho==0)=property('ice','rho','T',-1);

rho_i = density(s,Tf(1,:));

out = struct(...
    'rho',rho,...
    'rho_b',rho_i,...
    'liquid',liquid,...
    'rho_ice',property('ice','rho','T',-1),...
    'S',S,...
    'T',T,...
    's',s,...
    't',t,...
    'Tf',Tf,...
    'Tf0',Tf(1,:));
end

