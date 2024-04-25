function getProbeCalibration
%% function
cv = [ -5.1359    0.1399    0.0304   -0.0012    0.0032    0.0003];
fcn = @(cv,S,T) cv(1) + cv(2)*S + cv(3)*T + cv(4)*S.^2 + cv(5)*S*T + cv(6)*T.^2;

%%
% Load "data" 
load('labview_logs\probe\calibration\sweepT_-1_30degC_sweepS_0_60psu.mat')
home

% convert temperatures
data.T = conductivityTemperature(data.T);
% modify 0psu data such that near-freezing measurments are omitted
omit = data.t(:,1)/60<4;
data.S(omit,1) = inf;
data.T(omit,1) = inf;
% modify 10psu "
omit = data.t(:,2)/60<2.7;
data.S(omit,2) = inf;
data.T(omit,2) = inf;

% convert nans to inf for exclusion
data.S(isnan(data.S)) = inf;
data.T(isnan(data.T)) = inf;

valid = @(idx) ~isinf(data.S(:,idx));
% get 1st order polynomial fit coeffcients
fit1 = arrayfun(@(idx) coeffvalues(fit(data.T(valid(idx),idx),data.S(valid(idx),idx), ...
    'poly1')), ...
    1:length(data.s),'UniformOutput',false);

% get 3nd order polynomial fit coefficients
fit3 = arrayfun(@(idx) coeffvalues(fit(data.T(valid(idx),idx),data.S(valid(idx),idx), ...
    'poly3')), ...
    1:length(data.s),'UniformOutput',false);

fit1 = cell2mat(fit1');
fit3 = cell2mat(fit3');

% extrapolote regressions in salinity space
FIT3 = cell2mat(arrayfun(@(idx) coeffvalues(fit(data.s',fit3(:,idx),'poly3')),...
    1:width(fit3),'UniformOutput',false)');

s_ext = 0:10:100;
fit3_ex = cell2mat(arrayfun(@(idx) polyval(FIT3(idx,:),s_ext),1:height(FIT3),'UniformOutput',false)')';

% Create linear spaces
T = linspace(-5,30,1e3);
S1 = cell2mat(arrayfun(@(i) polyval(fit1(i,:),T),1:height(fit1),'UniformOutput',false)');
S2 = cell2mat(arrayfun(@(i) polyval(fit3(i,:),T),1:height(fit3),'UniformOutput',false)');
S3_ex = cell2mat(arrayfun(@(i) polyval(fit3_ex(i,:),T),1:height(fit3_ex),'UniformOutput',false)');


% Interpolate S
s0 = s_ext;%data.s;
S  = linspace(0,100,.5e3);
% meshes
[Tm0,Sm0]   = meshgrid(T,s0);
[Tm,Sm]     = meshgrid(T,S);
Si          = interp2(Tm0,Sm0,S3_ex,Tm,Sm);

% Make exclusion rules for liquidus
Tf = liquidus('T',S);
for i=1:numel(S)
    omit = T<Tf(i);
    Si(i,omit) = nan;
end

% Plotting
lbl = arrayfun(@num2str,data.s,'UniformOutput',false);
tiledlayout(2,1)
% Plot raw data
nexttile
plot(data.S,data.T)
addlabels('x','C_S (V)','y','T (degC)','title','Measurements')
L=legend(lbl,'Location','eastoutside','Box','off');
title(L,'NaCl Concentration (‰)')

% Plot lookup table
nexttile
imagesc(S,T,Si','AlphaData',~isnan(Si)')
set(gca,'Color','none')
addlabels('y','T (degC)','x','S (‰)','title','Lookup table')
addColorbar('cmap','balance','pivot',0,'title','C_S (V)')
axis xy

% export to struct
lookup = struct('CS',Si,'S',S,'T',T);
% save('labview_logs\probe\calibration\conductivityLookupTable',"lookup")

