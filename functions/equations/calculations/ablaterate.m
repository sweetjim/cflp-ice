function V = ablaterate(S,T,rho,varargin)
%% Ablation rate
% Calculates the ablation rate from a specified equation.
% 
% -------------------------------------------------------------------------
% %  Singular arguments:
% -------------------------------------------------------------------------
%   S   :   [struct]    -   [g/kg]
%       'i'             -   ice salinity
%       'b'             -   interface salinity
%       'w' or 'inf'    -   ambient salinity
% 
%   T   :   [struct]    -   [degC]
%       'i'             -   ice temperature
%       'b'             -   interface temperature
%       'w' or 'inf'    -   ambient temperature
% 
%   rho  :   [struct]   -   [kg / m3]
%       'i'             -   ice density
%       'b'             -   interface density
%       'w' or 'inf'    -   ambient density
% 
% -------------------------------------------------------------------------
% %  Singular arguments (Optional);
% -------------------------------------------------------------------------
%   'sf'  : [char]      - uses the salt-flux equation 
%       optional name-value argument 'theta' (default is theta = 0 deg)
% -------------------------------------------------------------------------
% % Name-value arguments (Optional):
% -------------------------------------------------------------------------
%       'theta' : [vec/double]  - slope parameter [degrees]
%       
%       'h_i'   : [double]      - ice thickness [m] (can be ignored unless
%                                   h~O(1e-6)) 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%
equation    = 'sf';
theta       = [];
h_i         = .1;

parseInput(varargin);
Gamma       = gammas('Ti',T.i,'h',h_i);
deltaS      = (S.b-S.w)./(S.i-S.b);
switch equation
    case {'sf','sfslope'}
        if isempty(theta)
            theta = 0;
        end
        V = cos(deg2rad(theta)).*deltaS.*(rho.w./rho.i)*Gamma.S;
end
%%
  function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'salt-flux','salt flux','SF','sf'}
                    equation = 'sf';
                case {'theta','slope'}
                    theta = namevalue;
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

