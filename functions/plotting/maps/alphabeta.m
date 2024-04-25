function [alpha,beta] = alphabeta(ds)
%% Alpha-beta map
% Generates a map of the alpha (thermal expansion coefficient) and beta
% (haline contraction coefficient) relationship for differing fluid
% temperatures and salinities at atmospheric pressure.
%%

if nargin==0
ds = 1e3;
end
S       = linspace(0,50,ds);
T       = linspace(-10,50,ds);
[s,t]   = meshgrid(S,T);
alpha   = sw_alpha(s,t,0);
beta    = sw_beta(s,t,0);

if nargout>0
return
end
clf
imagesc(S,T,(alpha./beta)')
hold on
contour(S,T,(alpha./beta)','k','ShowText','on')
hold off
axis xy
addColorbar('latex','fs',15,'title','$\alpha/\beta$','cmap','balance','pivot',0)
addlabels('latex','y','$S$ (g/kg)','x','$T$ ($^\circ$ C)','fs',15)
end

