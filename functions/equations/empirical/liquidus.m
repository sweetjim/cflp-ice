function OUT = liquidus(varargin)
%% Liquidus curve
% Determines the liquidus temperature (interface temperature / freezing
% point) for fluid with a salt concentration (in [g/kg]) at a given depth
% (in [m]).
%
% Examples:
%  liquidus('T',300) = -49.0625     
%   is the interface temperature at 300 [g/kg] at 0 [m]
% 
% liquidus('T',300,'z',1e4) = -41.4525
%   is the interface temperature at 300 [g/kg] at -10 [km]
% -------------------------------------------------------------------------
% % Name-value arguments:
% -------------------------------------------------------------------------
%   'S'     :   [char]  - calculates the interface salinity
%       T [vec]         - temperature at which to evaluate
% 
%   'T'     :   [char]  - calculates the interface temperature
%       S [vec]         - salinity at which to evaluate
% 
% -------------------------------------------------------------------------
% % Name-value arguments (Optional):
% -------------------------------------------------------------------------
%   'z'     :   [char]  - depth of liquidus
%       z [vec]         - depth in meters
% 
%   'coarg' :   [char]  - coargument to input
%       S || T [vec]    - coargument vector (must be same size as input
%                           argument)
% -------------------------------------------------------------------------
% %  Singular arguments (Optional):
% -------------------------------------------------------------------------
%   'rho'   : [char]    - calculates the interface variable and outputs the
%                           seawater density.
% 
%   'linear': [char]    - calculates the interface variable from the linear
%                           approximation (default is polynomial).
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%
DENSITY = false;
method  = 'S';
equation= 'poly';
var     = 0;
z       = 0;
coarg   = [];
parseInput(varargin);
z       = abs(z);

%% Equation
lambda  = 7.61e-4;                  % [degC / m]
switch equation
    case 'poly'
        Sb  = @(T) -18.7*T-0.519*T.^2-0.00535*T.^3;     % salinity-temperature relation (Vancoppenolle1 et al, 2018)
    case 'linear'
        T0      = 0;%8.32e-2;          % [degC]
        Gamma_S = -5.73e-2;         % [degC]
        Sb  = @(T) 1/Gamma_S.*(T-T0);
end
t       = linspace(0,-50,1e5);
sb      = Sb(t);

% TL = @(x) (x.*(374000 + 10380*x + 107*x.^2))/20000
%% Calculation
OUT = zeros(1,length(var));
for i=1:length(var)
    switch method
        case {'S','Sb'}
            OUT(i) = Sb(var(i))-lambda*z;
        case {'T','Tb'}
            OUT(i) = t(find(sb>var(i),1,'first')) + lambda*z;
    end
    if DENSITY
        switch method
            case 'T'
                T = OUT(i);
                S = var(i);
            case 'S'
                S = OUT(i);
                T = var(i);
        end
        OUT(i) = density(S,T);
    end
end

if ~isempty(coarg)
    if numel(var)~=numel(coarg)
        error('Coargument vector must be the same size as the input vector')
    end
    if strcmp(method,'T')
        map = coarg>6;
        OUT(map) = 0;
    end
end

%% Input parser
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'S','Sb'}
                    method  = 'S';
                    var     = namevalue;
                case {'T','Tb'}
                    method  = 'T';
                    var     = namevalue;
                case {'z'}
                    z       = namevalue;
                case {'rho','density'}
                    DENSITY = true;
                case {'linear','poly','polynomial'}
                    equation = items{m};
                case 'coarg'
                    coarg = namevalue;
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
