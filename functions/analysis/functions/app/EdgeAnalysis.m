function analysis=EdgeAnalysis(args,varargin)
%% Get Edge Analysis
% Uses the results from the GUI program 'tuner' and applies it to every
% image listed in the active repository (i.e. 'args.files') using the
% script 'applytuner'. At each time step ('step') the ice-edge, a
% vectorized quantity 'h', is added to the time series vector 'H'. Once
% complete, a surface type plot is generated displaying the results.
%
% -------------------------------------------------------------------------
% %  Parameters:
% -------------------------------------------------------------------------
%  args: [struct] (Required)
%   The output from GUI program 'tuner' (hitting 'save' generates this
%   variable).
%
% -------------------------------------------------------------------------
% %  Singular arguments (Optional):
% -------------------------------------------------------------------------
%
%  rerun: [char]
%   Activation variable for re-running the script. For use when new 'args'
%   variable is called.
%
%  surf, waterfall, contour: [char]
%   Plotting types; default is 'surf'.
%
%  height, gradient: [char]
%   C-axis data plotting types; default is 'height'.
%
%  sec, min, hour: [char]
%   Delta-t to extract from file timestamps (for time vector); default is
%   'min'.
%
% -------------------------------------------------------------------------
% % Name-value arguments (Optional):
% -------------------------------------------------------------------------
%  loop: [int, vector]
%   Elements to loop through and construct H from. For use when trialling
%   potentially unstable edge detection setups.
%
%  step: [int]
%   Time step to loop through and construct H from. Same as above.
%
%  ds: [double]
%   Delta-x or Delta-z in terms of cm/pixel.
%
%  buffer: [double]
%   Percentage quantity (a value between 0 and 50) relative to the maximum
%  z-value that is used to delete "spiked-data". If set to 5, i.e. 5% of
%  max(z), then any  "spiked-data" (erroneous edge detection) within 5% of
%  the z-boundaries ([min(z)(1+0.05), max(z)(1-0.05)]) will be cropped.
%  Note, this rescales 'z' and 'h'.
%  Default is 5.
% 
% -------------------------------------------------------------------------
% % Outputs:
% -------------------------------------------------------------------------
%   analysis: [struct]
%       fields:
%           z,t: [double, vec]
%           Depth and time vectors.
%           
%           h: [double, matrix]
%           Ice thickness data.
% 
%           dt,dt_str: [int / double, char]
%           Time step and units.
% 
%           root, savepath: [char]
%           Root folder (experiment) and the location of the saved
%           workspace (i.e. analysis, args, and others to come).
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% Folder setup

[~,name,~]  = fileparts(args.files(1).folder);
root        = 'results';
savepath    = fullfile(root,name);

%%  Initialization and parsers
buffer      = 5; % percent
cam         = [110 60];
colormode   = 'height';
mode        = 'Surf';
dt_str      = 'min';
step        = 1;
rerun       = 0;
Zspan       = 1.1;

warning off
parseInput(varargin);
loop            = 1:step:length(args.files);
[args.out,h,~]  = applytuner(args,'all',1);

Zlength     = length(h);
ds          = Zspan/Zlength;
z           = linspace(0,Zlength*ds,Zlength);

persistent H T dt Z EXPT

if isempty(EXPT)
    EXPT = name;
end

if ~rerun 
    if strcmp(EXPT,name)
        rerun = 0;
    else
        rerun = 1;
        EXPT = name;
    end
end

if ~isempty(H)
    cond = sum(size(H))==0;
else
    cond = true;
end
%% Generate time-series
if (isempty(H)&&isempty(T))||cond||rerun
    disp('Building time-series')
    H       = zeros(Zlength,length(loop));
%     size(H)
    Hpf     = H;
%     parfor i=loop
    for i=loop
        [~,h,~]     = applytuner(args,'all',i);
        h(h<0)      = 0;
        Hpf(:,i)    = h;
    end

    spikes = true;
    H = Hpf;
    H = H-mean(H(:,end));
    H(H<0) = 0;
else
    spikes = false;
end
%% Spikes near z-lims
if spikes
    z_bounds    = buffer*1e-2.*[1 -1]+[min(z) max(z)]; % buffer region at zlims
    BOUNDS      = ~(z>z_bounds(1)&z<z_bounds(2));
    BOUNDS      = repmat(BOUNDS,[length(loop) 1])';
    
    hlim        = mean(max(H(:,1:2)));
    
    if sum(sum(H>hlim))==0
        hlim=mean(max(H));
    end
    
    Herase      = H>hlim;
    h           = imdilate(Herase,strel('rectangle',[80 1]));
    [~,i]       = max(z(imbinarize(sum(h.*BOUNDS,2))));

    Z   = z(i:end);
    H   = H(i:end,:);
    
    %% Barrier removal
    dt  = getdt(args.files,step+1,dt_str);
    T   = 0:dt*step:dt*(size(H,2)-1);
    H   = H(:,2:end);
    T   = T(2:end);
end

%% Outputs


analysis = struct(...
    'z',Z,...
    't',T,...
    'cam',cam,...
    'loop',loop,...
    'colormode',colormode,...
    'mode',mode,...
    'h',H*ds,...
    'H',H,...
    'ds',ds,...
    'dt',dt,...
    'dt_str',dt_str,...
    'root',name,...
    'savepath',savepath);

args.root       = name;
args.savepath   = savepath;
assignin('base','analysis',analysis);

if isfile(strcat(savepath,'.mat'))
    save(savepath,'analysis','-append')
else
    save(savepath,'analysis')
end

%% Functions
    function parseInput(varargin)        
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                %% Name arguments
                case 'rerun'
                    rerun = true;
                case {'Default','Gradient','dht','dhtt'}
                    colormode = lower(items{m});
                case {'Waterfall','Surf','Contoured'}
                    mode    = lower(items{m});
                case {'seconds','minutes','hours','days'}
                    dt_str  = items{m};
                %% Name-value arguments
                case 'loop'
                    loop    = namevalue;
                case 'step'
                    step    = namevalue;
                case 'Zspan'
                    Zspan   = namevalue;
                case 'ds'
                    ds      = namevalue;
                case {'buffer','zbuffer'}
                    buffer  = namevalue;
                case 'cam'
                    cam     = namevalue;
            end
            m = m+1;
            
            if m>length(items)
                break;
            end
        end
        function out = namevalue
            out = items{m+1};
            m   = m+1;
        end
    end
end

