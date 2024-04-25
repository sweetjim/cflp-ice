
function getLayerVelocity(ddlayer,varargin)
%% Get Double Diffusive Layers
% Uses the results from 'getDDLayers' and calculates the vertical velocity
% of each layer as a function of time.
% -------------------------------------------------------------------------
% %  Parameters:
% -------------------------------------------------------------------------
%  analysis: [struct] (Required)
%   The output from program 'getEdgeAnalysis'.
% 
% -------------------------------------------------------------------------
% %  Singular arguments (Optional):
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% % Name-value arguments (Optional):
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% % Outputs:
% -------------------------------------------------------------------------
%   ddlayer: [struct]
%       fields:
%           method: [char]
%           The detection method used ('gradient' or 'peak').
% 
%           mode: [char]
%           The detection type used ('max' or 'min').
%  
%           resize: [int, vector]
%           The resizing coefficients in either [t] or [z,t].
%  
%           z,t: [double, vec]
%           Depth and time vectors.
%           
%           h: [double, matrix]
%           Ice thickness data.
% 
%           eta: [double, matrix]
%           Method-mode-based double diffusive layer data.
% 
%           layers: [double, vec]
%           Spatially averaged layer thicknesses (time series)
% 
%           spikes: [logical, vec]
%           Indexes of 'layers' that correspond to the time-gradient of
%           'layers' exceeding one standard deviation from the mean.
% 
%           root, savepath: [char]
%           Root folder (experiment) and the location of the saved
%           workspace (i.e. analysis, args, and others to come).
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%%
%% Initialization and parsers           
fig_title = 'Double-Diffusive Layer Velocity Analysis';
persistent fig 
type = 'cmocean';
[cmap,rerun,pivot,invert,smooth_t] = parseInput(varargin);
openFig(fig_title)
if rerun;fig = getfig;end
clf

if rerun
    if isempty(fig)
        fig = getfig;
    elseif ~isvalid(fig)
        fig = getfig;else
        clf(fig)
    end
else
    if isempty(fig)
        fig = getfig;
    elseif ~isvalid(fig)
        fig = getfig;
    else
        clf(fig)
    end
end

bmap    = ddlayer.binaryeta;
z       = ddlayer.z;
t       = ddlayer.t;
ds      = max(z)/length(z);
labels  = bwlabel(bmap')';
lblmax  = max(labels,[],'all');
velocity_label  = sprintf('%s (cm/%s)','$\partial_z\eta$',ddlayer.dt_str);
meanvel_label   = sprintf('%s (cm/%s)','$\langle\partial_z\eta\rangle$',ddlayer.dt_str);
time_label      = sprintf('t (%s)',ddlayer.dt_str);


%% Velocity calculation (distance mapping)
displacement    = zeros(length(t),lblmax);
meanz           = zeros(lblmax,1);
velocity        = displacement;
for i =1:lblmax
    VAR                 = labels==i;
    LT                  = getLayerThickness(islocalmax(VAR));
    doubles             = sum(LT==0,1)>1;
    newLT               = islocalmax(getLayerThickness(LT(:,doubles)==0));
    LT(:,doubles)       = getLayerThickness(newLT);
    
    dist                = max(LT,[],1)*ds;
    if ~isempty(smooth_t)
       dist             = smooth(dist,smooth_t); 
    end
    dist(dist==Inf)     = NaN;
    dist0               = dist(find(~isnan(dist),1));
    meanz(i)            = mean(z(mean(VAR,2)>0));
    displacement(:,i)   = dist-dist0;
    velocity(:,i)       = gradient(displacement(:,i),mean(diff(t)));   
end
for i=1:length(meanz)
    finaldisp(i)   = displacement(find(~isnan(displacement(:,i)),1,'last'),i);
end
clims           = [min(meanz) max(meanz)];
meanvelocity    = nanmean(velocity,1);
meanvelfit      = polyfitn(meanz,meanvelocity,1);
displacfit      = polyfitn(meanz,finaldisp,1);
%% Plotting
switch type
    case 'default'
        CMAP = colormap(cmap);
    case 'cmocean'
        if ~isempty(pivot)
            CMAP = cmocean(cmap,'pivot',pivot);
        else
            CMAP = cmocean(cmap);
        end
end
CMAP = imresize(CMAP,[length(meanz) 3]);
if strcmp(invert,'invert')
    CMAP = flipud(CMAP);
end

tile = tiledlayout(2,3,'TileSpacing','compact');
title(tile,'Vertical velocity and displacement of each ($k$) double diffusive layer',...
    'interpreter','latex','FontSize',15)
% -----------------------------------------------------------------------------
% DISPLACEMENT MAP
% -----------------------------------------------------------------------------
clear n
n(1) = nexttile(tile);
plot(n(1),t,displacement)
addlabels(n(1),...
    'x',time_label,...
    'y','$(z-z_{0,k})$ (cm)',...
    'title','Displacement',...
    'latex','fs',15)
ticks = xticks;
ticks(end+1) = round(mean(ddlayer.t(1)),2,'decimals');
set(n(1),'XTick',sort(unique(ticks)))
grid(n(1),'on')
addColorbar('ax',n(1),'label','$\langle z\rangle$','latex','fs',15,...
    'lineplot','limits',clims,cmap,'levels',length(meanz),...
    invert,'pivot',pivot,'invisible')
h = get(n(1),'Children');
legend(n(1),h(1),'$\langle z_k \rangle$',...
    'Location','northwest','interpreter','latex')

% -----------------------------------------------------------------------------
% VELOCITY MAP
% -----------------------------------------------------------------------------

n(2) = nexttile(tile,4);
plot(n(2),t,velocity)
addlabels(n(2),...
    'x',time_label,...
    'y',velocity_label,...
    'title','Velocity',...
    'latex','fs',15)
grid(n(2),'on')
addColorbar('ax',n(2),'label','$\langle z\rangle$','latex','fs',15,...
    'lineplot','limits',clims,cmap,'levels',length(meanz),'invisible',...
    invert,'pivot',pivot)


% -----------------------------------------------------------------------------
% MEAN VELOCITY 
% -----------------------------------------------------------------------------
n(3) = nexttile(tile,5);
scatter(n(3),meanz,meanvelocity,36,CMAP,'filled')
hold(n(3),'on');
plot(n(3),meanz,polyval(meanvelfit.Coefficients,meanz),'--k')
hold(n(3),'off');
xlim(n(3),[min(z) max(z)])
grid(n(3),'on')
addlabels(n(3),...
    'x','$z$ (cm)',...
    'y',meanvel_label,...
    'title','Mean-velocity',...
    'latex','fs',15)
legend(n(3),...
    {'$\langle z_k \rangle$',...
    sprintf('Fit (%s = %.1f)','$R^2$',meanvelfit.R2)},...
    'Location','northwest','interpreter','latex')

% -----------------------------------------------------------------------------
% FINAL DISPLACEMENT
% -----------------------------------------------------------------------------

n(4)        = nexttile(tile,2);
scatter(n(4),meanz,finaldisp,36,CMAP,'filled','s')
hold(n(4),'on');
plot(n(4),meanz,polyval(displacfit.Coefficients,meanz),'--k')
hold(n(4),'off');
addlabels(n(4),...
    'x','$z$ (cm)',...
    'y','$z_f$ (cm)',...
    'title','Final Displacement',...
    'latex','fs',15)
xlim(n(4),[min(z) max(z)])
grid(n(4),'on')
legend(n(4),...
    {'$\langle z_k \rangle$',...
    sprintf('Fit (%s = %.1f)','$R^2$',displacfit.R2)},...
    'Location','northwest','interpreter','latex')

% -----------------------------------------------------------------------------
% PHASE PORTRAIT 
% -----------------------------------------------------------------------------

n(5) = nexttile(tile,[2 1]);
numbins     = 20;
X           = reshape(displacement,[1 numel(displacement)]);
Y           = reshape(velocity,[1 numel(velocity)]);
remove      = imbinarize(isnan(X)+isnan(Y));
X(remove)   = [];
Y(remove)   = [];
HEAT        = heatCdata(X,Y,numbins);
x           = linspace(min(X),max(X),numbins);
y           = linspace(min(Y),max(Y),numbins);
% [values,~]  = hist3([X;Y]',numbins.*[1 1]);
% values(values==0) = NaN;
scatter(n(5),X,Y,36,HEAT,'filled')
% imagesc(x,y,values),shading flat
addColorbar('ax',n(5),'cmap','amp','limits',[min(HEAT) max(HEAT)],...
    'label',sprintf('Neighborhood count (%.4f cm %s %.4f %s)',...
    diff(x([1 2]))*1e2,'$\times$',...
    diff(y([1 2])),ddlayer.dt_str),...
    'latex','fs',15)
set(n(5),'ColorScale','log')

addlabels(n(5),...
    'x','$(z-z_{0,k})$ (cm)',...
    'y',velocity_label,...
    'title','Phase portrait',...
    'latex','fs',15)
grid(n(5),'on')
% h = get(n(5),'Children');
% legend(n(5),h([2 1]),...
%     {'$\langle z_k \rangle$',...
%     sprintf('Fit (%s = %.2f)','$R^2$',phasefit.R2)},...
%     'Location','northwest','interpreter','latex')

for i=1:length(n)
box(n(i),'on')
end
%% Outputs                              
layervelocity = struct(...
    'z',z,...
    't',t,...
    'layers',length(meanz),...
    'meanz',meanz,...
    'velocity',velocity,...
    'displacement',displacement,...
    'meanvelocity',meanvelocity,...
    'root',ddlayer.root,...
    'savepath',ddlayer.savepath);
save(layervelocity.savepath,'ddlayer','-append')
assignin('base','layervelocity',layervelocity)

%% Functions

    
%% Input parser                         
    function [cmap,rerun,pivot,invert,smooth] = parseInput(varargin)
        cmap        = 'amp';
        rerun       = false;
        pivot       = [];
        smooth      = 5;
        invert      = '';
        
        %%
        cmoceanmaps = {'thermal','balance','haline','delta','solar','curl'...
            'ice','diff','gray','tarn',...
            'oxy','deep','dense','phase',...
            'algae','matter','turbid','topo',...
            'speed','amp','tempo','rain'};
        colormaps = {'parula','jet','hsv','hot','cool','spring',...
            'spring','summer','autumn','winter',...
            'gray','bone','copper','pink','lines',...
            'colorcube','prism','flag','white'};
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                %% Name arguments
                case 'rerun'
                    rerun = true;
                case cmoceanmaps
                    cmap = items{m};
                    type = 'cmocean';
                case colormaps
                    cmap = items{m};
                    type = 'default';
                case {'flip','invert'}
                    invert  = 'invert';
                case 'pivot'
                    pivot   = namevalue;
            end
            m = m+1;
            if m>length(items);break;end
        end
        
        function out = namevalue
            out = items{m+1};
            m   = m+1;
        end
    end
    function fig = getfig
        fig = figure('Name',fig_title,'WindowStyle','docked','NumberTitle','off');
    end
end