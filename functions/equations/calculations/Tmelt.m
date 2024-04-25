function Tef = Tmelt(T)
%% Meltwater temperature (effective)
% Calculates the "effective" temperature of released meltwater (see
% Magorrian & Wells, 2016, p. 6 or 4675)
% 
% Given as
% 
% -------------------------------------------------------------------------
% %  Singular arguments:
% -------------------------------------------------------------------------
% T_ef = Tb - L/cw + cs/cw*(Ti-Tb)
% 
%   T   :   [struct]    -   [degC]
%       'i'             -   ice temperature
%       'b'             -   interface temperature
%       'inf'           -   ambient temperature
%%

L  = property('L');
cw = property('water','cp','T',T.inf);
ci = property('ice','cp','T',T.i);
Tef = T.b-L/cw+ci/cw.*(T.i-T.b);
    
end

