function OUT = property(varargin)
%% Properties of ice and water
% 
% Contains information about the temperature dependence (when applicable)
% on density, thermal conductivity, specific heat capacity, thermal
% diffusivity, haline diffusivity.
% -------------------------------------------------------------------------
% %  Singular arguments:
% -------------------------------------------------------------------------
%   type: [char]    (default is 'ice')
%       'ice'       - selects ice properties
%       'water'     - selects water properties
%       'L'         - latent heating                            [J / kg]
%       'LL'        - generic modified latent heating           [-84.5 degC] 
%       'StT'       - Themal Stanton number                     [1.1e-3]
%       'StS'       - Haline Stanton number                     [3.1e-5]
% 
%   var: [char]     (default is 'kappaT')
%       'rho'       - density                                   [kg / m3]
%       'K'         - thermal conductivity                      [W / m K]
%       'cp'        - specific heat capacity                    [J / kg K]  
%       'kappaT'    - thermal diffusivity                       [m / s2]
%       'kappaTf'   - dynamic fluid thermal diffusivity (water) [m / s2]
%       'kappaS'    - haline diffusivity                        [m / s2]
%       'nu'        - kinematic viscosity (water only)          [m / s2]
% -------------------------------------------------------------------------
% %  Name-value arguments (Optional):
% -------------------------------------------------------------------------
%   T: [double,vec]  - specific temperature associated with variable of
%                      interest. 
% 
%   S: [double,vec]  - specific salinity associated with kappaTf as
%       kappaTf = kappaT(T)/(cp(T) * rho(S,T))
%       MUST HAVE 'T' specified.
% -------------------------------------------------------------------------
% % Outputs:
% -------------------------------------------------------------------------
% 
%   OUT: [struct]   - if T is unspecified
%       T           - vector of temperatures
%       var         - vector of var(T)
% 
%   OUT: [double]   - if T is specified
%       var         - double of var(T)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%
type = 'ice';
var = 'kappaT';
Tin = [];
S   = [];
parseInput(varargin);

if strcmp(var,'kappaTf')&&~isempty(S)&&~isempty(Tin)
    kappaT  = property('water','kappaT','T',Tin);
    c       = property('water','cp','T',Tin);
    
    if numel(Tin)>1&&numel(S)==1
        %% T-array
        rho = arrayfun(@(i) density(S,Tin(i)),1:length(Tin));
%         for i=1:length(Tin)
%             rho(i)  = density(S,Tin(i));
%         end
        OUT     = (kappaT./c).*1./rho;
    elseif numel(Tin)==1&&numel(S)>1
        %% S-array
        rho = arrayfun(@(i) density(S(i),Tin),1:length(S));
%         for i=1:length(S)
%             rho(i)   = density(S(i),Tin);
%         end
        OUT     = (kappaT./c).*1./rho;
    else
        %% ST-matrix or value
        rho     = density(S,Tin);
        OUT     = (kappaT./c)'.*1./rho;
    end
    
    return
elseif strcmp(var,'kappaTf')&&(isempty(S)||isempty(Tin))
    error('Specify T and S.')
end
if numel(Tin)>1
    t = Tin;
else
    t = linspace(-100,100);
end
kappaS = @(t)10.^(-5.144+0.0127.*t);    % Wasburn (1926). International Critical Tables of Numerical Data
kappaT = @(k,rho,cp) k./(rho.*cp);      % [m2 / s]
L       = 3.35e5;                       % [J / kg]

switch var
    case 'L'
        OUT = L;
        return
    case 'LL'
        OUT = -84.5;
        return
    case 'StT'
        OUT = 1.1e-3;
        return
    case 'StS'
        OUT = 3.1e-5;
        return
end



%% Empirical data
switch type
    case 'ice'
        T       = [0,-5,-10,-15,-20,-25,-30,-35,-40,-50,-60,-70,-80,-90,-100];
        %[deg C]
        rho     = [916.2,917.5,918.9,919.4,919.4,919.6,920,920.4,920.8,921.6,922.4,923.3,924.1,924.9,925.7];
        % [kg / m3]
        K       = [2.220,2.250,2.300,2.340,2.390,2.450,2.500,2.570,2.630,2.760,2.900,3.050,3.190,3.340,3.480];
        % [W / m K]
        cp      = 1e3*[2.050,2.027,2,1.972,1.943,1.913,1.882,1.851,1.818,1.751,1.681,1.609,1.536,1.463,1.389];
        % [J / kg K]
        KAPPAT   = kappaT(K,rho,cp);
    case 'water'
        Tnu     = [0.01,10,20,25,30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,240,260,280,300,320,340,360];
        % [deg C]
        nu      = 1e-6*[1.7918,1.3065,1.0035,0.8927,0.8007,0.6579,0.5531,0.4740,0.4127,0.3643,0.3255,0.2938,0.2677,0.2460,0.2123,0.1878,0.1695,0.1556,0.1449,0.1365,0.1299,0.1247,0.1206,0.1174,0.1152,0.1143];
        % [m2 / s]
        Tcp     = [0.01,10,20,25,30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,240,260,280,300,320,340,360];
        % [deg C]
        cp      = 1e3*[4.2199,4.1955,4.1844,4.1816,4.1801,4.1796,4.1815,4.1851,4.1902,4.1969,4.2053,4.2157,4.2283,4.2435,4.2826,4.3354,4.4050,4.4958,4.6146,4.7719,4.9856,5.2889,5.7504,6.5373,8.2080,15.004];
        % [J / kg K]
        TK      = [0.01,10,20,30,40,50,60,70,80,90,99.6];
        % [deg C]
        K       = 1e-3*[555.75,578.64,598.03,614.50,628.56,640.60,650.91,659.69,667.02,672.88,677.03];
        % [W / m K]
        TKAPPAT = [0,5,10,20,25,30,50,75,100];
        % [deg C]
        KAPPAT  = 1e-6*[0.13200,0.13500,0.13800,0.14300,0.14600,0.14800,0.15500,0.16200,0.16800];
end

%% Fitting

switch type
    case 'ice'
        Targin = T;
        switch var
            case 'kappaT'
                argin = KAPPAT;
            case 'kappaS'
                error('No data for haline diffusivity of ice')
            case 'kappaTf'
                error('Why are you calculating this?!')
            case {'cp','Cp'}
                argin  = cp;
            case {'k','K'}
                argin  = K;
            case 'rho'
                argin  = rho;
            case 'nu'
                error('No data for kinematic viscosity of ice')
        end
    case 'water'
        switch var
            case 'kappaT'
                Targin  = TKAPPAT;
                argin   = KAPPAT;
            case 'kappaS'
                Targin  = t;
                argin   = kappaS(t);
            case 'kappaTf'
                Targin  = TKAPPAT;
                argin   = KAPPAT;
            case {'cp','Cp'}
                Targin  = Tcp;
                argin   = cp;
            case {'k','K'}
                Targin = TK;
                argin   = K;
            case 'nu'
                Targin  = Tnu;
                argin   = nu;
            case 'rho'
                error('Use density function')
        end
end


[range,argout] = fitpoly(Targin,argin,3);

%% Output parser
if isempty(Tin)||numel(Tin)>1
    OUT = struct('T',t);
    switch var
        case 'kappaT'
            OUT.kappaT = argout;
        case 'kappaS'
            OUT.kappaS = argout;
        case {'cp','Cp'}
            OUT.cp      = argout;
        case {'k','K'}
            OUT.k       = argout;
        case 'rho'
            OUT.rho     = argout;
        case 'nu'
            OUT.nu      = argout;
    end
    OUT.validity = sprintf('T = %.3f to %.3f degC',range(1),range(2));
    
    if numel(Tin)>1
       OUT  = argout; 
    end
else
    idx = find(t>=Tin,1,'first');
    OUT     = argout(idx);
end

%% Polynomial fit

    function [range,out] = fitpoly(Targ,arg,n)
        range   = [min(Targ) max(Targ)];
        fit     = polyfitn(Targ,arg,n);
        out     = polyval(fit.Coefficients,t);
    end


%% Input parser
    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case {'ice','water'}
                    type = items{m};
                case {'rho','kappaT','kappaTf','kappaS',...
                        'cp','Cp','k','K','L','LL','nu',...
                        'StT','StS'}
                    var = items{m};
                case 'T'
                    Tin = namevalue;
                case 'S'
                    S   = namevalue;
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

