function Tef = Teff(T)
%% Effective temperature of released meltwater (sensible heating)
% Calculates the "apparent" temperature of ice if it were a fluid by
% factoring in the latent heat compoonent from Jenkins (2011) & Magorrian &
% Wells (2016) as
% 
% T_ef = Tb - L/cw + cs/cw * (Ti-Tb)
% -------------------------------------------------------------------------
% %  Singular arguments:
% -------------------------------------------------------------------------
%
%   T   :   [struct]    -   [degC]
%       'i'             -   ice temperature         [vec, double]
%       'b'             -   interface temperature   [vec, double]
%       'w' or 'inf'    -   ambient temperature     [vec, double]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% 
L   = property('L');
cw  = property('water','cp','T',T.w);
ci  = property('ice','cp','T',T.w);
Tef = T.b - L./cw+ci./cw.*(T.i-T.b);
end

