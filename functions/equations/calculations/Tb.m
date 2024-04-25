function TB = Tb(rho,S,T,varargin)
%% Interface temperature
% Calculates the interface temperature from the three-equation balance
% -------------------------------------------------------------------------
% %  Singular arguments:
% -------------------------------------------------------------------------
%   rho  :   [struct]   -   [kg / m3]
%       'i'             -   ice density
%       'b'             -   interface density
%       'w' or 'inf'    -   ambient density
%
%   S   :   [struct]    -   [degC]
%       'i'             -   ice temperature
%       'b'             -   interface temperature
%       'w' or 'inf'    -   ambient temperature
%
%   T   :   [struct]    -   [degC]
%       'i'             -   ice temperature
%       'b'             -   interface temperature
%       'w' or 'inf'    -   ambient temperature
% 
% -------------------------------------------------------------------------
% %  Singular arguments (Optional):
% -------------------------------------------------------------------------
%   'negative'          - returns only the negative sign component of the
%                           quadratic (default is both signs)
%   'positive'          - returns only the positive sign component of the
%                           quadratic
%% Equation
% Given as
%
% Tb = 1/(1+qR)*(-Twiq + T0*(1+qR) pm sqrt((Twiq - T0*(1+qR))^2 - 4*(1+qR)*(qLST GammaS Sw - T0 * Twiq)))
%
%
% where:
% qR    = (rhoi/rhow)*(ci/cw)*(gammaTi/gammaT)
%
% qLST  = (rhoi/rhow)*(gammaS/gammaT)*(L/cw)
%
% Twiq  = -Tw-qR*Ti-qLST
%%
output  = 'both';
parseInput(varargin);

T0      = 0;%8.32e-2;          % [degC]
Gamma_S = 5.73e-2;          % [degC]
Gamma   = gammas('Ti',T.i);

c       = struct(...
    'i',property('ice','cp','T',T.i),...
    'w',property('water','cp','T',T.w));

qi      = rho.i.*c.i*Gamma.Ti;
qw      = rho.w.*c.w*Gamma.T;
qL      = rho.w.*Gamma.S*property('L');
qiw     = qw+qi;
FTq     = T.w.*qw+T.i.*qi+qL;

A       = qiw;
B       = -(T0.*qiw+FTq);
C       = T0.*FTq+qL.*Gamma_S.*S.w;

TB      = 1./(2*A).*(-B+sqrt(B.^2 - 4*A.*C).*[-1 1]');
switch output
    case 'negative'
        TB = TB(1,:);
    case 'positive'
        TB = TB(2,:);
    otherwise
        return
end

%%
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'pos','+','+ve'}
                    output = 'positive';
                case {'neg','-','-ve'}
                    output = 'negative';
            end
            m = m+1;
            if m>length(items);break;end
        end
        function out = namevalue
            out = items{m+1};
            m   = m+1;
        end
    end
end

