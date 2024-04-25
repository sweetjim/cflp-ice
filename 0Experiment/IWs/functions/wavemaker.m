function [y,t,details,export]=wavemaker(N,varargin)
%% WAVEMAKER generates waveforms for use within the wavemaker LabView program
% Parameters:
%   N           -   [double]
%                   buoyancy frequency [rad/s]
%   Le          -   [double]
%                   e-folding length (exp(-x^2/Le^2), DEFAULT is 8.5cm)
%   'theta'     -   [double]
%                   wave angle [deg]
%   'fr'        -   [double]
%                   lateral Froude number
%   'amplitude' -   [double]
%                   wave amplitude (mm)
%   'export'    -   [char]
%                   exporting routine
%   'label'     -   [string]
%                   export label prefix e.g. 'folder/filename'
%   'silent'    -   [string]
%                   suppresses output messages and figure display
%% Inputs
theta       = []; % wrt horizontal
fr          = [];
omega       = [];
amplitude   = [];
export      = false;
suppress    = false;
label       = [];
Le          = 8.5e-2;%3e-2;
parseInput(varargin)


fileloc = 'IWs\wavemaker_profiles\';
%% Assumed geometry
k0 = 1/Le;

if ~isempty(amplitude)
    fr = (amplitude*1e-3)*k0/N;
else
    amplitude = fr*N/k0*1e3;
end

if numel(amplitude)>1
    fprintf('Detected: amplitude parameter space\n')
    for i=1:numel(amplitude)
       wavemaker(N,'omega',omega,'theta',theta,'amplitude',amplitude(i),'export','label',label);
    end
    return
end
%% Input
% theta = @(N,omega) rad2deg(asin(omega/N));
if isempty(omega)
    omega   = N*sin(deg2rad(theta));
elseif isempty(theta)
    theta   = rad2deg(asin(omega/N));
end
period  = 2*pi/omega;
kappa   = sqrt(1+theta^2);
cg      = N^2/kappa/omega*cos(theta)*sin(theta)*[sin(theta) -cos(theta)];


% decay factor
nu      = 1e-6;
alpha   = nu*(k0*N).^3./(2*omega^3*sqrt(N^2-omega^2)); % (Lighthill 1978; Waves in fluids)
x       = linspace(0,1e2,1e3);
decay   = exp(-alpha*x); 
decay   = x(find(decay<1e-2,1,'first')); % distance to 99% viscous attentuation
%% Time resolution
% max 600 pts
t_res   = [50 100 200 500]*1e-3; % 50ms precision
t_res   = t_res(find(round(period./t_res)<600,1,'first'));

t           = linspace(0,round(period,3),period/t_res);
y           = sin(omega*t).*amplitude;%/(24.5/4);
y([1 end])  = 0;


if ~suppress
    f = findall(0,'Type','Figure','Name','Wavemaker profile');
    if isempty(f)
        f = figure('Name','Wavemaker profile','IntegerHandle','off');
    end
    set(0,'CurrentFigure',f)
    plot(t,y,'.')
    addlabels('x','Time (s)','y','Amplitude (mm)')


    s{1}=sprintf('θ:\t\t\t%.3g (° wrt. hztl)\nω:\t\t\t\t%.3g (rad/s)',theta,omega);
    s{2}=sprintf('ω/N:\t\t%.3g',omega/N);
    s{3}=sprintf('Fr:\t%.3g\nAmplitude:\t\t%.3g (mm)',fr,amplitude);
    s{4}=sprintf('2π/ω:\t\t\t%.3g (s)\ndt:\t\t\t\t%.3g (ms)',period,t_res*1e3);
    delete(findall(gcf,'Type','annotation'))
    annotation(gcf,'textbox','EdgeColor','none','String',s)
    shg
end
%%
details             = struct;
details.theta       = theta;
details.omega       = omega;
details.N           = N;
details.amplitude   = amplitude;
details.fr          = fr;
details.period      = period;
details.t_res       = t_res;
details.decay       = decay;

if ~suppress
    fprintf('\n')
    fprintf('θ:\t\t\t\t%.3g (° wrt. hztl)\nω:\t\t\t\t%.2g (rad/s)\n',theta,omega)
    fprintf('ω/N:\t\t\t%.3g\n',omega/N)
    fprintf('Fr:\t\t\t\t%.2g\nAmplitude:\t\t%.2g (mm)\n',fr,amplitude)
    fprintf('2π/ω:\t\t\t%.3g (s)\ndt:\t\t\t\t%.2g (ms)\n',period,t_res*1e3)
end
%% Output
if ~isempty(label)
    if contains(label,'today')
        label = strrep(label,'today',datestr(datetime('today'),'DD_mm_YY'));
    end
    if contains(label,{'\','/'})&&~isfolder(label)
        fileloc = fullfile(fileloc,label);
        if ~isfolder(fileloc)
            mkdir(fileloc)
        end
    end
end

tag = strcat(strrep(sprintf('N_%.2frads-theta_%.2fdeg-amp_%.1fmm-Fr_%.3g-dt_%ims',...
        N,theta,amplitude,fr,t_res*1e3),'.','_'),'.csv');
tag     = [strcat(datestr(datetime('today'),'yy_mm_dd'),'#') tag];
filename = fullfile(fileloc,tag);
if export
    writetable(...
        table(y'),...
        filename,...
        'WriteVariableNames',false)
end
if nargout==0
    clear y t details
    return
end
%% functions
function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                case 'Le'
                    Le    = namevalue;
                case {'m','m0'}
                    val = namevalue;
                    Le = 1/val;
                case 'omega'
                    omega = namevalue;
                case 'theta'
                    theta = namevalue;
                case 'amplitude'
                    amplitude = namevalue;
                case {'froude','fr'}
                    fr  = namevalue;
                case 'export'
                    export = true;
                case 'label'
                    label = namevalue;
                case 'silent'
                    suppress = true;
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

