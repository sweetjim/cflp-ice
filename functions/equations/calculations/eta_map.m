function deltaT = eta_map(S,T,z)
%% Delta_T map
% Generates a map of the boundary layer thickness from Tanny & Tsinober
% (1988)'s expression for 'eta', which is the boundary layer thickness in
% the Rayleigh number expression:
% 
% Ra = g*alpha*DeltaT*eta^3/(kappaT*nu)
% 
% where eta = Delta_T = alpha*(Tb-Tw)/(-beta*ds0/ds)
%%
[s,t]   = meshgrid(S.w,T.w);
dsdz    = s/(max(z)-min(z));
alpha   = sw_alpha(s,t,0);
beta    = sw_beta(s,t,0);
deltaT  = alpha.*(ones(size(t)).*T.b-t)./(-beta.*dsdz);
end