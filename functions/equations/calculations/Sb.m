function SB = Sb(rho,T)
%% Interface salinity
% Calculates the interface salinity from the three-equation balance
%  Parameters:
%   rho  :   [struct]   -   [kg / m3]
%       'i'             -   ice density
%       'b'             -   interface density
%       'w' or 'inf'    -   ambient density
% 
%   T   :   [struct]    -   [degC]
%       'i'             -   ice temperature
%       'b'             -   interface temperature
%       'w' or 'inf'    -   ambient temperature
%% Equation
% Given as
% 
% Sb = 1/(2 GammaS) * [-(T0 + q_q + q_q / q_LST * (qR Ti - Tw)) pm sqrt((T0 + q_q + q_q /
% q_LST* (qR Ti - Tw))^2 + 4 GammaS q_q)]
% 
% where:
% q_q = q_LST/(1 + q_R)
% 
% q_LST = L/cw * gammaS/ gammaT
% 
% q_R   = rho_i/rho_w * ci/cw * gammaTice / gammaT
% 
%%

T0      = 8.32e-2;          % [degC]
lambda  = 7.61e-4;          % [degC / m]
Gamma_S = 5.73e-2;          % [degC]
Gamma   = gammas('Ti',T.i);

c       = struct(...
    'i',property('ice','cp','T',T.i),...
    'w',property('water','cp','T',T.w));

qR      = rho.i./rho.w * c.i/c.w * Gamma.Ti/Gamma.T;
qLST    = property('L')./c.w*Gamma.S/Gamma.T;
qq      = qLST./(1+qR);

A       = Gamma_S;
B       = 1./(1+qR).*(qR.*T.i+T.w)-T0-qq;
C       = qq.*S.w;

SB = 1./(2*A).*(-B+sqrt(B.^2 - 4*A.*C).*[-1 1]');

end

