function [rho,comp] = liquidusDensity(state,val)
%% LIQUIDUSDENSITY returns the density of seawater at the liquidus at reference sea level
% 
% Example: 
% 
% >> [rho,S] = liquidusDensity('T',-1)
% rho =
%    1.0145e+03 (kg / m3)
% S =
%    18.1818 (psu)

%%

persistent Sliquidus
getIdx = @(t) find(Sliquidus<t,1,'first');
s = linspace(0,50,1e4);
switch state
    case 'T'
        if isempty(Sliquidus)
            Sliquidus = zeros(size(s));
            for i=1:length(Sliquidus)
                Sliquidus(i) = liquidus('T',s(i));
            end
        end
        comp    = s(getIdx(val));
        rho     = density(comp,val);
    case 'S'
        comp    = liquidus('T',val);
        rho     = density(val,comp);
end

