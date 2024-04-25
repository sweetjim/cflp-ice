function OUT = gammas(varargin)
%% Gamma (transfer coefficients)
% Calculates / states the transfer coefficients.
% 
% Parameters: (Optional, Name-value) 
%   'h'     : [double]  -   the ice thickness (default is 0.1 [m])
%   'Ti'    : [double]  -   the ice temperature (default is -5 [degC])
% 
% Outputs:
%   'gamma' : [struct]
%       'T'             - heat coeffecient
%       'Ti'            - heat in ice coefficient
%       'S'             - salt coefficient

%%
h   = .1; % [m]
Ti  = -5; % [degC]
parseInput(varargin);

OUT   = struct(...        % []
    'T',[2.2e-2],...        % heat transfer coefficient
    'Ti',property('ice','kappaT','T',Ti)/h,...             % heat in ice transfer coefficient (kappaT / h_ice (thickness))
    'S',[6.2e-4]);          % salinity transfer coefficient

%%
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'h_ice','thickness','H'}
                    h   = namevalue;
                case {'T','Ti'}
                    Ti  = namevalue;
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

