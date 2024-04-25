function getEdgeAnalysis(args,varargin)
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
buffer  = 5; % percent
[step,colormode,mode,rerun,dt_str] = parseInput(varargin);
loop  = 1:step:length(args.files);

args.out = args.in;
Zspan   = 1.1;
Zlength = size(args.out,1);
ds      = Zspan/Zlength;


z       = linspace(0,Zlength*ds,Zlength);
%% Loop (Edge detection time series)

fig_title = 'Edge Analysis';
openFig(fig_title)

persistent H T fig dt EXPT Z

if isempty(EXPT)                    % Assign EXPT
    EXPT = name;
elseif strcmp(EXPT,name)&&~rerun    % If same EXPT and no rerun called, assume no rerun
    rerun = 0;
elseif ~strcmp(EXPT,name)&&~rerun   % If new EXPT and no rerun called, set as rerun
    rerun = 1;
end

if rerun||(isempty(H)&&isempty(T))
    fig     = getfig;
    dt      = getdt(args.files,step*2,dt_str);
    
    T       = 0:dt*step:dt*(loop(end)-1);
    H       = zeros(size(args.out,1),length(loop));
    
    parfor i=loop
        [~,h,~]     = applytuner(args,'all',i);
        h(h<0)      = 0;
        H(:,i)      = h;
%         displayProgress('Progress',i,loop(1),loop(end))
    end

    spikes = true;
    H = H-mean(H(:,end));
    H(H<0) = 0;
else
    spikes = false;
end

%% Spikes near z-lims
if spikes
    z_bounds    = buffer*1e-2.*[1 -1]+[min(z) max(z)]; % buffer region at zlims
    BOUNDS      = ~(z>z_bounds(1)&z<z_bounds(2));
    BOUNDS      = repmat(BOUNDS,[length(T) 1])';
    
    hlim        = mean(max(H(:,1:2)));
    
    if sum(sum(H>hlim))==0
        hlim=mean(max(H));
    end
    
    Herase      = H>hlim;
    h           = imdilate(Herase,strel('rectangle',[80 1]));
    [p,i]       = max(z(imbinarize(sum(h.*BOUNDS,2))));
%     line(xlim,p.*[1 1])
    
    Z   = z(i:end);
    H0  = H;
    H   = H(i:end,:);
    
    %% Barrier removal
    H   = H(:,2:end);
    T   = T(2:end);
    tplot = T;
end
%% Plotting

if rerun
    if ~isvalid(fig)
        fig = getfig;
    else
        clf(fig)
    end
else
    switch dt_str
        case 'sec'
            tplot = T*60;
        case 'min'
            tplot = T;
        case 'hour'
            tplot = T/60;
    end
    if ~isvalid(fig)
        fig = getfig;
    else
        clf(fig)
    end
end


ax  = axes(fig);

switch colormode
    case 'height'
        C       = H*ds;
        clabel  = '$h$ (m)';
    case {'gradient','dht'}
        [C,~]   = gradient(H*ds,tplot,Z);
        clabel  = sprintf('%s (m/%s)','$\partial_t h$',dt_str);
    case 'dhtt'
        [ht,~]  = gradient(H*ds,tplot,Z);
        [C,~]   = gradient(ht);
        clabel  = sprintf('%s (m/%s)','$\partial_{tt} h$',dt_str);
end

switch mode
    case 'waterfall'
        waterfall(ax,Z,tplot,H'*ds,C')
        view(110,60)
    case 'surf'
        surf(ax,Z,tplot,H'*ds,C')
        shading interp
        view(110,60)
    case 'contour'
        %%
        surf(ax,Z,tplot,H'*ds,C')
        hold(ax,'on')
        caxis(max(caxis).*[-1 1])
        cmocean('balance','pivot',0,2^5)
        contour3(ax,Z,tplot,H'*ds,'k')
        hold(ax,'off')
        shading interp
        %         view(90,90)
        view(110,60)
end
c                       = colorbar;
c.Label.String          = clabel;
c.Label.Interpreter     = 'latex';
c.TickLabelInterpreter  = 'latex';
c.FontSize              = 15;

ylabel(sprintf('%s (%s)','$t$',dt_str),'Interpreter','latex')
xlabel('$z$ (m)','Interpreter','latex')
zlabel('$h$ (m)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',15)

axis(gca,'tight')
%% Outputs
[C,~] = gradient(H*ds,tplot,Z);
analysis = struct(...
    'z',Z,...
    't',T,...
    'h',H*ds,...
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
    function [step,colormode,mode,rerun,dt_str] = parseInput(varargin)
        
        colormode   = 'height';
        mode        = 'surf';
        rerun       = false;
        dt_str      = 'min';
        step        = 1;
        
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                %% Name arguments
                case {'default','gradient','dht','dhtt'}
                    colormode = items{m};
                case {'waterfall','surf','contour'}
                    mode    = items{m};
                case 'rerun'
                    rerun   = true;
                case {'sec','min','hour'}
                    dt_str  = items{m};
                %% Name-value arguments
                case 'loop'
                    loop    = namevalue;
                case 'step'
                    step    = namevalue;
                case 'ds'
                    ds      = namevalue;
                case {'buffer','zbuffer'}
                    buffer  = namevalue;
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

