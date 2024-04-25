function TSL(varargin)
%% Temperature-Salinity-Liquidus
% Generates a plot of the density of water due to haline and thermal
% concentrations with the liquidus curve (freezing point).
% 
% -------------------------------------------------------------------------
% %  Singular arguments (Optional):
% -------------------------------------------------------------------------
%   'ice'   :   [char]  -    draws the liquidus slope and sets it
%   to an ice density at -5 degC
% 
%   'density' : [char]  - map of seawater density
%  
%   'Teff'  :   [char]  - recalculates the density of the ice to account
%   for latent heating and the temperature dependence on ice's density.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%
if nargin==0
    type = 'ice';
elseif nargin==1
   switch varargin{1}
       case {'Teff','Tef'}
           type = 'Teff';
       otherwise
           type = 'ice';
   end
end
%%
clf
tsl =    interface(0);
tsl.rho = smooth2a(tsl.rho,10,1);
% Plotting
rhostep = 980:20:max(tsl.rho,[],'all');
%%
% imAlpha=ones(size(tsl.S));
% imAlpha(isnan(tsl.rho))=0;

map0 = tsl.rho.*tsl.liquid;
switch type
    case 'ice'
        map = map0;
        map(map==0) = tsl.rho_ice;
    case 'Teff'
        
        T   = struct('i',-5,...
            'w',tsl.t);
        for i=1:length(tsl.t)
           T.b(i) = tsl.t(find(tsl.liquid(:,i)>0,1,'first')); 
        end
        TB          = (repmat(Teff(T),[length(tsl.s) 1]))';
        rho_i       = tsl.rho_ice.*(sw_BETA(tsl.S,tsl.T-TB)-sw_ALPHA(tsl.S,tsl.T-TB));
        map_Teff    = (tsl.rho_ice+rho_i).*~tsl.liquid;
        map         = map0+map_Teff;
    case 'density'
        map = tsl.rho;
end

imagesc(tsl.s,tsl.t,map);axis xy %,'AlphaData',imAlpha

set(gca,'color',1*[1 1 1]);
hold on;
% plot(tsl.s,tsl.Tf(1,:),'-c','LineWidth',1)
% plot(tsl.s,tsl.Tf(2:end,:),'--c','LineWidth',1)
switch type
    case 'ice'
        rho = tsl.rho.*tsl.liquid;
        rho(rho==0) = NaN;
        [C,h] = contour(tsl.s,tsl.t,rho,rhostep,'w','ShowText','on','LabelSpacing',2e2);
        clabel(C,h,'Color','w')
    case 'Teff'
        rho = tsl.rho.*tsl.liquid;
        rho(rho==0) = NaN;
        map_Teff(map_Teff==0) = NaN;
        minmax = [min(map_Teff,[],'all') max(map_Teff,[],'all')];
        ice_step = round(minmax(1):diff(minmax)/5:minmax(2),1,'decimal');
        
       
        [C,h] = contour(tsl.s,tsl.t,rho,rhostep,'w','ShowText','on','LabelSpacing',2e2);
        clabel(C,h,'Color','w')
        [C,h] = contour(tsl.s,tsl.t,map_Teff,ice_step,'r','ShowText','on');
        clabel(C,h,'Color','r')
    otherwise
        [C,h] = contour(tsl.s,tsl.t,tsl.rho,rhostep,'w','ShowText','on','LabelSpacing',2e2);
        clabel(C,h,'Color','w')
end
hold off;

addlabels(gca,'latex','x','S (g/kg)','y','T ($^\circ$C)','fs',12)
addColorbar(gca,'latex','label','$\rho(S,T)$ (kg m$^{-3}$)',...
    'cmap','ice','flip',...
    'levels',round(diff(caxis)/2)+1,'limits',[960 Inf],'fs',12)
latexformat('ratio',[1 .5])
out2latex('rho_liquidus')
end

