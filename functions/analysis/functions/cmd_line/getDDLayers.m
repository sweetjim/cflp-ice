function getDDLayers(analysis,varargin)
%% Get Double Diffusive Layers
% Uses the results from 'getEdgeAnalysis' and detects the locations and
% thicknesses of the double-diffusive convective layers.
%
% -------------------------------------------------------------------------
% %  Parameters:
% -------------------------------------------------------------------------
%  analysis: [struct] (Required)
%   The output from program 'getEdgeAnalysis'.
% 
% -------------------------------------------------------------------------
% %  Singular arguments (Optional):
% -------------------------------------------------------------------------
%  rerun: [char]
%   Activation variable for re-running the script. For use when new
%   'analysis' variable is called or when 'min / max' are changed. 
%
%  min, max: [char]
%   Detection type; either the maxima of the ice-face scollups
%   (corresponding to double-diffusive layer boundaries) or minima of
%   ice-face scollups (centrepoint of layer boundaries).
%   Default is 'max'.
%
%  gradient, peak: [char]
%   Peak detection method; either use a gradient-based method (preferred)
%   which takes the second numerical gradient of the spatial dimension (z),
%   or use an explicit peak finding algorithm ('islocalmin / islocalmax')
%   along the spatial dimension for each instance of time. 
%   Default is 'gradient'.
%
%  sec, min, hour: [char]
%   Delta-t to extract from file timestamps (for time vector); default is
%   'min'.
%  
%  reversed: [char]
%   Plotting argument that "inverts" the distance map. E.g.:
%   Default:    C-map = 'eta'
%   Reversed:   C-map = 'eta(max(eta)-1)'
%
%  noedge: [char]
%   Plotting argument that sets layers that extend to the boundaries to
%   zero thickness.
%   Default is false.
% 
% -------------------------------------------------------------------------
% % Name-value arguments (Optional):
% -------------------------------------------------------------------------
%  smooth: [int, vector]
%   2D smoothing factors in [z,t] on the ice-thickness data 'h'. Expressed
%   as either a scalar or vector (integers only).
%   Default is [7,3].
% 
%  minz, zmin, min_z: [double]
%   Minimum distance between peaks in meters for peak-acquisition
%   algorithms.
%   Default is 0.01.
%
%  resize: [int, vec]
%   Resizing coefficients on 'h','z', and 't' (applied after 2D smoothing)
%   in the form of either [t] or [z,t].
%   Default is 2 (stretch only 't' related data).
% 
%  thres: [double]
%   Threshold factor between 0 and 1; either sets for:
%   'gradient', the normalized gradient rejection factor, 
%   i.e., ||'hzz'||>thres,
%       where ||'hzz'|| isin [-1,1] but for modes 'min' or 'max' the target
%       range is always between [0,1].
%   'peak' the zerotime-normalized binary-labelled-peak rejection factor,
%   i.e., 
%                                   / < thres, rejected
%       ||length(eta_k)||_zerotime |
%                                   \ >= thres, kept
%       where k is an integer between 1 and N (number of detected peaks),
%       length(eta_k) is the time-length of the kTH peak, and
%       ||.||_zerotime is the time-length normalization factor relative to
%       the time at which 'h=0' or 'h(t(end))~=0' (the zerotime).
%   Default is (0.9,0.3) for ('gradient','peak').
% 
%  t_offset: [int]
%   Time offset from zero as an index (or image number).
%   Default is 10.
%
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

%% Initialization and parsers           
plottype = [];
noedge   = [];
fig_title = 'Double-Diffusive Layers Analysis';
persistent fig OUTgrad OUTpeak type_mode type_method BINARYMAP
[mode,min_z,prominance_w,smooth_x,smooth_y,method,resize,rerun,thres,t_offset] = parseInput(varargin);

if isempty(type_mode)
    type_mode = mode;
end
if isempty(type_method)
    type_method = method;
end
if ~strcmp(type_mode,mode)&&strcmp(type_method,method)
    rerun = true;
end

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

dt      = analysis.dt;
z       = analysis.z;
ds      = max(z)/length(z);

z_range = z>=0;
min_th  = find(z>min_z,1,'first');

z       = z(z_range);
t       = analysis.t(t_offset:end);                        % Assume that first image includes barrier
h       = smooth2a(analysis.h(z_range,t_offset:end),smooth_x,smooth_y);


if sum(resize)>0
    if numel(resize)==1
        h = imresize(h,[size(h,1) size(h,2)*resize]);
        t = imresize(t,[1 length(t)*resize]);
    else
        h = imresize(h,size(h).*resize);
        z = imresize(z,[1 length(z)*resize(1)]);
        t = imresize(t,[1 length(t)*resize(2)]);
    end
end

H       = h;
H(H==0) = NaN;

switch mode
    case 'max'
        clabel = '$\eta$ (m)';
    case 'min'
        clabel  = '$\eta($max$(\eta)-1)$ (m)';
end

%% Peak acquisition                     
LMax        = zeros(size(h));
LMin        = zeros(size(h));

melt        = load(analysis.savepath,'melt');
if isempty(fieldnames(melt))
    getMeltPrediction(analysis)
    load(analysis.savepath,'melt');
else
    melt = melt.melt;
end

normalize                   = melt.zerotime(z_range);
normalize(isnan(normalize)) = max(normalize);
normalize                   = normalize/max(normalize);

switch method                           
    case 'gradient'
        if isempty(OUTgrad)||rerun
            %% Gradient and labelling           
            [~,hz]  = gradient(H,t,z);
            [~,hzz] = gradient(hz,t,z);
            
            switch mode
                case 'max'
                    H_pa    = -hzz;
                case 'min'
                    H_pa    = hzz;
            end
            
            new_lim = 5*nanstd(H_pa,[],'all');
            H_pa(H_pa>new_lim)  = NaN;
            H_pa(H_pa<-new_lim) = NaN;
            
            H_pa    = smooth2a(H_pa,0,5);
            IN      = H_pa>0;
%             IN      = ~bwareaopen(~IN,300);
%             IN      = bwareaopen(IN,500);
            LBL     = bwlabeln(IN')';
            %% Peak/trough acquisition          
            IN   = LBL;
            for i=1:max(IN,[],'all')
                
                VAR = (IN==i);
                LT  = getLayerThickness(VAR,'window',5,'normalize','invert');
                if i==1;ETA = zeros(size(LT));end
                ETA = ETA + LT;
                
            end
            %% Length argument                  
            %             holes   = OUT>thres;
            %             OUT     = (OUT.*(bwareaopen(holes,500)));
            ETA     = removal(ETA>thres);
            ETA     = binarySmooth(ETA);
            dist_map = zeros(size(h'));
            
            for i=1:length(t)
                dist_map(i,:) = bwdist(ETA(:,i)>.9);
            end
            
            BINARYMAP = ETA;
            OUTPUT = dist_map*ds;%.*~isnan(H');           
            %% Output                           
            OUTgrad = OUTPUT;
        end
    case 'peak'
        if isempty(OUTpeak)||rerun
            %% Local minima / maxima acquisition
            H_pa    = H;
            
            for i=1:size(h,2)
                HT          = H_pa(:,i);
                localmax    = islocalmax(HT,'MinSeparation',min_th,'ProminenceWindow',prominance_w);
                localmin    = islocalmin(HT,'MinSeparation',min_th,'ProminenceWindow',prominance_w);
                LMax(:,i)   = localmax;
                LMin(:,i)   = localmin;
            end
            
            switch mode
                case 'max'
                    IN      = LMax;
                    INcomp  = LMin;
                case 'min'
                    IN      = LMin;
                    INcomp  = LMax;
            end
            %% Connectivity / labelling         
            % Compose a binary labelling map of the detected peaks (with
            % dilation)
            dilate = strel('rectangle',[5 3]);
            dilate = dilate.Neighborhood;
            ETA     = imdilate(IN,dilate);
            %% Length argument                  
            % If the time-length of a detected peak is less than 20% of the
            % normalized time-length (when h_peak=0) then it will be removed.
            
            % Normalized label deletion
            ETA     = removal(ETA);
            
            % Re-cast as binary image
            ETA = binarySmooth(ETA);
            dist_map = zeros(size(h'));
            for i=1:length(t)
                dist_map(i,:) = bwdist(ETA(:,i));
            end
            
            OUTPUT      = dist_map*ds.*~isnan(H');
            
            BINARYMAP   = ETA;
            OUTpeak     = OUTPUT;
        end
end
%% Layer thickness                      

switch method
    case 'peak'
        ETA     = OUTpeak;
    case 'gradient'
        ETA     = OUTgrad;
        ETA(isnan(H'))   = 0;
        IN      = ETA==0;
        
        for i=1:length(t)
            thick_max   = islocalmax(IN(i,:)*ds);
            width_t(i)  = mean(diff(z(thick_max)));
        end
end
%% End point fix
if noedge
    for i=1:length(t)
        lmin      = islocalmin(ETA(i,:));
        pts       = z([find(lmin,1,'first') find(lmin,1,'last')]);
        ETA(i,z<pts(1))=0;
        ETA(i,z>pts(2))=0;
    end
end
thickness = zeros(length(t),1);
for i=1:length(t)
    thick_max       = islocalmax(ETA(i,:)*ds);
    thickness(i)  = mean(diff(z(thick_max)));
end

ETA(isnan(H)')=NaN;

delta_eta_z = thickness;
grad_deta   = gradient(delta_eta_z,t);
spikes      = (grad_deta>std(grad_deta));

delta_eta_t = nanmean(ETA,1);
thickness   = islocalmax(delta_eta_t,'MinSeparation',find(z>mean(delta_eta_z),1,'first'));
mean_eta_t  = mean(diff(z(thickness)));


%% Plotting                             
switch plottype                         
    case 'default'
        PLOT = ETA;
    case 'reversed'
        PLOT = ETA.*(max(ETA,[],'all')-1);
        clabel = '$\eta($max$(\eta)-1)$ (m)';
end
clf
tile = tiledlayout(fig,2,2);

% -----------------------------------------------------------------------------
% ETA SURFACE
% -----------------------------------------------------------------------------

n(1) = nexttile(tile,[1 2]);%,[2 7]);

waterfall(n(1),z,t,h',PLOT);
switch plottype
    case 'reversed'
        colormap(n(1),flipud(parula))
    otherwise
        colormap(n(1),parula)
end
% switch method
%     case 'gradient'
%         hold on
%         w = waterfall(z,t(spikes),h(:,spikes)',h(:,spikes)'>0);
%         w.EdgeColor = 'w';
%         w.EdgeAlpha = 0.5;
%         hold off
% %         logfac = ceil(-log10(max(delta_H))+1);
% %         caxis([0 ceil(max(delta_H)*10^logfac)*10^-logfac])
% end
set(n(1),'ColorScale','log')
axis(n(1),'ij','tight')
view(10,50);
addColorbar('ax',n(1),'latex','label',clabel,'fs',15)
addlabels(n(1),'title','Layer thickness',...
    'x','$z$ (m)',...
    'y',sprintf('%s (%s)','$t$',analysis.dt_str),...
    'z','$h$ (m)',...
    'latex','fs',15)


% -----------------------------------------------------------------------------
% ETA SPATIAL AVERAGE
% -----------------------------------------------------------------------------

n(2) = nexttile(tile);%,9,[1 2]);
plot(n(2),t,delta_eta_z,t(spikes),delta_eta_z(spikes),'o')
addlabels(n(2),'title','Layer thickness (spatial-average)',...
    'y','$\langle\eta\rangle_z$ (m)',...
    'x',sprintf('%s (%s)','$t$',analysis.dt_str),...
    'latex','fs',15)
axis(n(2),'xy','tight')
ticks = yticks;
ticks(end+1) = round(mean(delta_eta_z),3,'decimals');
set(n(2),'YTick',sort(unique(ticks)))
line(n(2),xlim,mean(delta_eta_z).*[1 1],'linestyle','--','color','c')

child = get(n(2),'Children');
legend(n(2),child([1 2]),...
    {'$\langle\eta\rangle_z$',...
    '$\langle\partial_t \eta>\sigma_{\partial_t\eta}\rangle_z$'},...
    'Location','northwest','Interpreter','latex')

% -----------------------------------------------------------------------------
% ETA TIME AVERAGE
% -----------------------------------------------------------------------------


n(3) = nexttile(tile);%,19,[1 2]);
plot(n(3),z,delta_eta_t,z(thickness),delta_eta_t(thickness),'o')
% plot(n(3),z,delta_eta_t)
axis(n(3),'xy','tight')
ticks = xticks;
ticks(end+1) = round(mean_eta_t,3,'decimals');
set(n(3),'XTick',sort(unique(ticks)))

line(mean_eta_t.*[1 1],ylim,'linestyle','--','color','c')
child = get(n(3),'Children');
legend(n(3),child([1 2]),...
    {'$\langle\eta\rangle_t$','Layer edge'},...
    'Location','northwest','Interpreter','latex')

% ,...
%     '$\langle\eta\rangle_z$'},...

addlabels(n(3),'title','Layer thickness (time-average)',...
    'y','$\langle\eta\rangle_t$ (m)',...
    'x','$z$ (m)',...
    'latex','fs',15)
%% Outputs                              
ddlayer = struct(...
    'method',method,...
    'resize',resize,...
    'dt_str',analysis.dt_str,...
    't0',t_offset,...
    't',t,...
    'z',z,...
    'h',h,...
    'eta',ETA,...
    'binaryeta',BINARYMAP,...
    'layers',delta_eta_z,...
    'spikes',spikes,...
    'root',analysis.root,...
    'savepath',analysis.savepath);
save(ddlayer.savepath,'ddlayer','-append')
assignin('base','ddlayer',ddlayer)
%% Functions                            
    function OUT = removal(IN)
        
        switch method
            case 'gradient'
                OUT = bwdist(IN)<2;
                OUT = bwlabel(OUT')';
            case 'peak'
                OUT = bwlabel(IN')';
        end
        labels      = max(OUT,[],'all');
        
        for j=1:labels
            item        = OUT==j;
            highlight   = nanmean(IN.*item,2);
            range       = ([find(highlight>0,1,'first') find(highlight>0,1,'last')]);
            meanpos     = islocalmax(highlight,'MinSeparation',0);
            if numel(meanpos)>1
                meanpos = nanmean(range);
            end
            
            [~,c0] = find(item>0,1,'first');
            [~,cf] = find(item>0,1,'last');
            
            pos(j) = floor(meanpos);
            len(j) = diff([c0 cf])/length(t);
        end
        
        filt    = len./normalize(pos)';
        
        switch method
            case 'gradient'
                thres = 0;
        end
        cutoff  = filt>thres;
        remove  = filt.*(~cutoff)>0;
        
        for j=1:labels
            if remove(j)
                OUT(OUT==j) = NaN;
            end
        end
    end
    function OUT = binarySmooth(IN)
        bmap       = imdilate(IN,dilateMatrix('star',10,5));
        bmap_label  = bwlabel(bmap')';
        labels      = 1:max(bmap_label,[],'all');
        for j=1:max(bmap_label,[],'all')
            area(j) = bwarea(bmap_label==j);
        end
        area    = area/max(area);
        keep    = labels(area>std(area));
        bmap   = zeros(size(bmap));
        for j=labels
            if sum(j==keep)>0
                bkeep=bmap_label==j;
                bmap = bmap +  bkeep;
            end
        end
        OUT = bmap;
    end
%% Input parser                         
    function [mode,min_z,prominance_w,smooth_x,smooth_y,method,resize,rerun,thres,t_offset] = parseInput(varargin)
        
        mode            = 'max';
        prominance_w    = 1;
        smooth_x        = 7;
        smooth_y        = 3;
        min_z           = 1e-2;
        method          = 'gradient';
        resize          = 2;
        rerun           = false;
        thres           = 0.3;
        t_offset        = 10;
        plottype        = 'default';
        noedge          = false;
        
        %%
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                %% Name arguments
                case 'rerun'
                    rerun = true;
                case 'max'
                    mode  = items{m};
                case 'min'
                    mode  = items{m};
                case {'gradient','dhdz'}
                    method  = 'gradient';
                case {'peak'}
                    method  = 'peak';
                case {'reverse','reversed'}
                    plottype = 'reversed';
                case {'noedge','edgefix'}
                    noedge = true;
                %% Name-value arguments
                case 'smooth'
                    smooth = namevalue;
                    if numel(smooth)<1
                        smooth_x = smooth;
                        smooth_y = smooth;
                    else
                        smooth_x = smooth(1);
                        smooth_y = smooth(2);
                    end
                case {'zmin','minz','min_z'}
                    min_z           = namevalue;
                case {'prominance','prom','p'}
                    prominance_w    = namevalue;
                case {'resize'}
                    resize          = namevalue;
                case {'thres'}
                    thres           = namevalue;
                case {'toffset','timeoffset','offset'}
                    t_offset        = namevalue;
            end
            m = m+1;
            if m>length(items);break;end
        end
        
        switch method
            case 'gradient'
                if thres==0.3
                    thres = 0.9;
                end
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