function n = refractiveIndex(S,T,lambda)
%REFRACTIVEINDEX 
% Calculates the refractive index of seawater using the empirical formula derived from <a href="matlab:web('https://doi.org/10.1364/AO.34.003477','-browser')">Quan and Fry (1995)</a>. 
% Validity ranges:
%   T = [0 30] degC
%   S = [0 35] g/kg
%   lambda = [400 700] nm
%
% Parameters: S (required),T (required),lambda (optional, default 589.3nm)
% 
%%
if nargin<3
    lambda = 589.3;
end

% Empirical equation coefficient terms
n0 =  1.31405;
n1 =  1.779e-4;
n2 = -1.05e-6;
n3 =  1.6e-8;
n4 = -2.02e-6;
n5 =  15.868;
n6 =  0.01155;
n7 = -0.00423;
n8 = -4382;
n9 =  1.1455e6;

% Equation
n = n0+(n1+n2.*T+n3.*T.^2).*S + n4*T.^2+1./lambda.*(n5+n6*S+n7.*T)+n8./lambda.^2+n9./lambda.^3;

end

