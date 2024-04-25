function h = DDL_h(rho_i,rho_inf,z)
%% Double diffusive layer thickness
% Calculates the DDL thickness according to Huppert & Turner's (1980)
% empirical theory given as
% 
% h = (0.65pm0.06) (rho(T0,Sinf) - rho(Tinf,Sinf)) / (drho/dz)
% 
% Parameters: rho_i, rho_inf, z

drhodz = gradient(rho_inf,z);
h = 0.65 * (rho_i - rho_inf)./drhodz; 

end

