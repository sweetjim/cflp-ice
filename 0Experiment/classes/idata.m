classdef idata < handle
    % IDATA is a "super" class object for general use with the
    % ice-related experiments at the Australian National University's
    % Geophysical Fluid Dynamics Laboratory.
    %%
    properties (Access = public)
        Expt        % Experiment title                                  (string)
        Details     % Details of experiment                             (string)
        Path        % Path to experiment from cd                        (string)
        Folder      % Details of experiment directory contents          (class)
        Analysis    % Analysis attributes                               (class)
        IW          % IW attributes                                     (class)
        Controls    % GUI input control parameters                      (struct)
        ROI0        % Persistent ROI w/ rotation argument               (vector)
        ROI         % Current ROI                                       (vector)
        Mask        % Mask polygon                                      (vector)
        Memory      % Details for adaptive control memory               (class)
        Dimensions  % Physical dimensions of the image and time series  (struct)
        Image       % Unfiltered image output                           (array)
        Output      % Recent algorithm output                           (array)
        Edge        % Recent edge detection output                      (vector)
        Offset      % Edge offset (true zero ice thickness)             (vector)
        Images      % Loaded images                                     (struct)
    end
    properties (Hidden)
        times                       % Datetime array of image files
        offsetMask  = struct(...    % Additional edge offset (2D)
            'poly',[],...
            'vector',[])
        analysisMask = [];
        hasLoggedData   = false;    
        root            = 'D:\PhD\Data'
        dirs            = struct(...
            'photos','ablation_images',...
            'results','ablation_results',...
            'julabo','')
        matlabroot      = cd %'D:\PhD\MATLAB_Drive';
        cloudroot       = 'G:\My Drive\Data transfer\';
        zoffset         = 0.05;
        history = struct('image',[],'edge',[],'dilation',[],'connect',[],'output',[],'index',1)
    end
    
    methods
        %% CONSTRUCTOR
        function id = idata(varargin)
            %% Constructor for icedata class
            % -------------------------------------------------------------------------
            %  Parameters :
            % -------------------------------------------------------------------------
            %
            % folder: [char] (Optional)
            %   call: ('folder',folderobject)
            %   Creates an empty ICEDATA class (all attributes are empty)
            %   DEFAULT is empty.
            %
            % load: [char] (Optional)
            %   call: ('load',folderpath)
            %   Searches the root directory for a folder with the given path (under
            %   './ablation'), attempts to load in any previously saved
            %   ICEDATA files, otherwise it creates a "blank" ICEDATA
            %   object (only attributes EXPT, PATH, FOLDER are non-empty).
            %
            % -------------------------------------------------------------------------
            % -------------------------------------------------------------------------
            %%
            
            % Constructor
            readstate   = false;
            folder2read = {};
            loadExpt    = false;
            parseInput(varargin)
            if readstate
                readFolder(id,folder2read)
%                 convert2TIFF(id)
            elseif loadExpt
                loadExperiment(id,folder2read);
            end
            %% Input parser
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case 'folder'
                            readstate   = true;
                            folder2read = namevalue;
                        case 'load'
                            loadExpt    = true;
                            folder2read = namevalue;
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
    end
    methods
        %% MAIN FUNCTIONS
        function varargout      = parseimage(id,ax,imageIndex)
            %% Parse image
            % This subroutine applies the edge detection algorithm to a
            % singular image from the time series.
            
            if isempty(id.times)
                id.times   = datetime({id.Folder.Files.date});
            end
            %% Viewing mode toggle
            if nargin==2
                if isa(ax,'double')
                    imageIndex  = ax;
                    ax        = gca;
                end
            elseif nargin==1
                ax = gca;
                imageIndex = 1;
            end
                        
            % Run constructor
            alg         = icealgorithm(id);
            % Load image
            alg.loadImage(imageIndex,true);
            
            % Run sequencer
            try
                id.Output  = alg.Sequencer(imageIndex);
            catch ME
                errorHandler;
            end
            
            
            id.Image   = alg.original;
            id.Edge    = alg.iceedge;
            if nargout==1
               varargout{1} = id.Edge;
               return
            end
            
            views = id.Controls.viewing;
            
            
            logicalmap = logical(sum([...
                views.edge,...
                views.dilation,...
                views.connectivity]));
            
            % Draw lengths
            [x,z,ds] = applyMetric(id);
            
            
            % Handle OVERLAY
            if views.overlay
                if (~views.connectivity)&&(logicalmap)
                    imagesc(ax,x,z,id.Image.*~id.Output);
                else
                    imagesc(ax,x,z,id.Image);
                end
            else
                imagesc(ax,x,z,id.Output);
            end
            
            
            % Handle CLIM and colors
            if ~logicalmap||views.overlay
                if views.overlay
                    caxis(ax,[0 Inf])
                else
                    caxis(ax,id.Controls.values.image.clim)
                    if isfield(id.Controls.values.image,'transient')
                        if id.Controls.values.image.transient>0
                            caxis(ax,[-1 1]*10)
                        end
                    end
                end
            else
                caxis(ax,[0 1]);
            end
            
            axis(ax,'ij')
            colormap(ax,flipud(gray))
            
%             addlabels('ax',axes,...
%                 'title',sprintf(...
%                 't=%s',...

                title(ax,datestr(diff([id.times(1) id.times(imageIndex)]),'HH:MM:SS')) %#ok<*DATST> 
            
            % Handle ICEEDGE DATA
            if (views.connectivity)&&(~isempty(id.Edge))
                hold(ax,'on')
                plot(ax,id.Edge/ds,z,'r','LineWidth',2)
                hold(ax,'off')
            end
            
            %% OUTPUTS
            if nargout>1
                varargout{1} = id.Image;
                varargout{2} = x;
                varargout{3} = z;
                varargout{4} = diff([id.times(1) id.times(imageIndex)]);
                varargout{5} = id.Plume.maskEdge(imageIndex);
                varargout{6} = id.Edge/ds;
            end
            return
            
            %% OLD CODE TO BE INCLUDED
            switch app.viewing %#ok<*UNRCH> 
                case 'Original'
                    method      = 'input';
                    app.overlay = 'image';
                    
                case {'Dilate','Combined','Detection','Filter'}
                    method      = lower(app.ViewingModeDropDown.Value);
                    app.overlay = 'image';
                    
                otherwise
                    method      = 'all';
                    app.overlay = 'pipeline';
            end
            
            % Undelay toggle
            if strcmp(app.UnderlayoriginalSwitch.Value,'On')
                app.overlay     = 'overlay';
            end
            % Time-series neighbour check and plotting
            if strcmp(app.ChecktimeneigborsSwitch.Value,'On')
                j = 1;
                for i=app.timerange
                    [out,h]     = applytuner(app.args,i,method);
                    if i==app.timerange(1)
                        output      = zeros(size(out));
                        h_output    = zeros(size(out,1),length(app.timerange));
                    end
                    if ~isempty(h)
                        h_output(:,j) = h;
                    end
                    output = out+output;
                    j = j+1;
                end
                
                switch app.overlay
                    case 'pipeline'
                        cla(ax)
                        out = applytuner(app.args,round(app.ImagenumberSlider.Value),method,'input');
                        imagesc(ax,out)
                        hold(ax,'on')
                        plot(ax,(h_output),1:size(h_output,1))
                        hold(ax,'off')
                        
                        %% Ordering of lines
                        h = get(ax,'Children');
                        for i=1:length(h)
                            isline(i) = isa(h(i),'matlab.graphics.chart.primitive.Line');
                        end
                        h = h(isline);
                        %% Color order
                        cmap = cmocean('balance');
                        c = 1-imresize(cmap,[length(h) 3],'nearest');
                        for i=1:length(h)
                            h(i).Color =  c(i,:);
                        end
                        
                        %% Legend
                        if length(h)>1&&strcmp(app.ViewingModeDropDown.Value,'Connectivity')
                            try
                                legend(ax,...
                                    h([1 median(1:length(h)) end]),{'Next','Current','Previous'},...
                                    'location','northeast');
                            catch
                                if length(h)==2&&app.ImagenumberSlider.Value==1
                                    legend(ax,...
                                        h([1 end]),{'Next','Current'},...
                                        'location','northeast');
                                elseif length(h)==2&&app.ImagenumberSlider.Value==length(app.args.files)
                                    legend(ax,...
                                        h([1 end]),{'Current','Previous'},...
                                        'location','northeast');
                                end
                            end
                        end
                        
                    otherwise
                        imagesc(ax,output)
                end
                
                axis(ax,'tight')
                
            else
                applytuner(app.args,round(app.ImagenumberSlider.Value),method,app.overlay,'ax',ax);
            end
            %% Error handler
            function errorHandler
                % Catch new inputs that saved files do not already have
                
                for ii=1:length(ME.stack)
                    errorLine   = ME.stack(ii).line;
                    errorFile   = ME.stack(ii).file;
                    file        = splitlines(fileread(errorFile));
                    errorCode   = strrep(file{errorLine},' ','');
                    
                    if contains(errorFile,'PhD')
                        break
                    end
                end
                matlab.desktop.editor.openAndGoToLine(errorFile,errorLine);
                error('%s\nLine: %i\nCode: %s\n\n',ME.message,errorLine,errorCode)
                
                %% MISC
                
                if contains(ME.message,'Unrecognized field name')
                    
                    % Attempt to deconstruct code
                    try
                        code    = strrep(extractAfter(errorCode,'='),';','');
                        fields  = split(code,'.');
                        command = fields{1};
                        
                        % Read-back any previous definitions
                        beforeLine = file(1:errorLine);
                        for k=fliplr(1:(length(beforeLine)-1))
                            codeLine = strrep(file{k},' ','');
                            if contains(codeLine,strcat(command,'='))
                                % Assuming code is of the form:
                                % command = id.command1.command2;
                                parentfields = ...
                                    split(...
                                    strrep(...
                                    extractAfter(codeLine,'=')...
                                    ,';',''),...
                                    '.');
                                break
                            end
                        end
                        
                        % Look at commands
                        if any(contains(parentfields,{'Controls','controls'}))
                            tmpControls = icedata;
                            tmpControls.loadControls([],[]);
                            subcommand = parentfields{end};
                            defaultValue = eval(...
                                sprintf('tmpControls.Controls.values.%s.%s;',subcommand,fields{end}));
                            if numel(defaultValue)==1
                                statement = sprintf('id.Controls.values.%s.%s=%s;',...
                                    subcommand,fields{end},num2str(defaultValue));
                            elseif numel(defaultValue)==2
                                statement = sprintf('id.Controls.values.%s.%s=[%s];',...
                                    subcommand,fields{end},num2str(defaultValue));
                            end
                            eval(statement)
                            fprintf("Added / changed a variable in 'icedata'\n")
                            fprintf('\t>> %s\n',statement)
                            
                        end
                        
                    catch
                        []; %#ok<*VUNUS> % Debug
                    end
                else
                    errorCode
                    ME
                end
                
            end
        end
        function [H,varargout]  = parseAnalysis(id,varargin)
            %% Parse analysis
            % This subroutine applies the edge detection analysis filters
            % algorithm to the data in the time series.
            % -------------------------------------------------------------------------
            %  Parameters : SINGULAR
            % -------------------------------------------------------------------------
            % time_fmt: [char]
            %  Time format ('Seconds'/'Minutes'/'Hours'/'Days')
            %  DEFAULT IS MINUTES 
            %   
            % caxis: [char]
            %   Output field ('thickness','gradient','acceleration')
            %   Time-derivative fields are with respect to time. 
            %   DEFAULT IS 'thickness'.
            % 
            % usemask: [char]
            %   Attempts to apply the Analysis-Mask if it exists.
            %   DEFAULT IS FALSE.
            %
            % -------------------------------------------------------------------------
            %  Parameters : NAME-VALUE
            % -------------------------------------------------------------------------
            %   'sample': [int]
            %   Sampling frequency at the specified time format (natural
            %   numbers)
            %   DEFAULT IS 1 
            %
            %  'maxtime': [double]
            %   Maximum time. 
            %   DEFAULT IS INF
            % -------------------------------------------------------------------------
            %  Outputs:
            % -------------------------------------------------------------------------
            %   H: Ice thickness data (H = H(Z,T))
            %   
            %   VARARGOUT:
            %   Z:          z-axis depth vector
            %   T:          x-axis time vector
            %   CLABEL_FMT: c-axis formatting structure
            %   omit:       c-axis omit binary map
            % -------------------------------------------------------------------------
            % -------------------------------------------------------------------------
            %%
            useparent   = false;
            parent      = strrep([extractBefore(id.Path,'short'),'long'],'images','results');
            smooth_zt   = [0 0];
            limiters    = [1 1];
            method      = 'NaN';
            dilation    = 0;
            t_fmt       = 'Minutes';
            caxis_fmt   = 'Thickness';
            usemask     = false;
            sampling    = 1;
            maxT        = inf;
            parseInput(varargin)
            
            offset = id.Offset;
            if useparent&&~contains(id.Path,'long')
                try
                    datap    = load([parent,'.mat']).data;
                    data     = datap.Analysis;
                    [~,~,ds] = id.applyMetric;
                    offset   = datap.Offset/ds;
                    [];
                catch
                    error('Invalid parent folder')
                end
            else
                data        = id.Analysis;
            end
            if isempty(data)
               error('No ice-edge analysis has been acquired for:\n%s\n',id.Path)
            end
            hasControls = any(contains(fieldnames(data),'controls'));
            
            % Apply metric tensor
            [X,Z,ds] = id.applyMetric(data.H);
%             if ~isempty(id.Dimensions.ds)
%                 Z       = (data.Z-1)/ds;
%                 
%             end
            data.H  = data.H/ds;
            [T,CLABEL_FMT] = timeAxes(t_fmt,data);
            
            if hasControls
               %% Gather controls
               items = fieldnames(data.controls);
               for i=1:numel(items)
                   labels = eval(sprintf('fieldnames(data.controls.%s)',items{i}));
                   for j=1:numel(labels)
                       value = eval(sprintf('data.controls.%s.%s',items{i},labels{j}));
                       switch items{i}
                           case 'smooth'
                               smooth_zt(j) = value;
                           case 'cutoff'
                               switch labels{j}
                                   case 'h'
                                       limiters(1)  = value/100;
                                   case 'method'
                                       method       = value;
                                   case 'dtdh'
                                       limiters(2)  = value/100;
                                   case 'dilation'
                                       dilation     = value;
                               end         
                       end
                   end                   
               end
            end
            smooth_zt = round(smooth_zt);
            
            %% APPLY FILTERS
            H           = data.H;
            nanned = imdilate(isnan(H),dilateMatrix("rect",20,2));
            H(nanned)=nan;
            H = fillmissing(H,'movmean',40,2);
            if any(limiters<1)
                % Maxima
                max_h       = max(H,[],'all');
                [dt,~]      = gradient(H,T,data.Z);
                max_dtdh    = max(abs(dt),[],'all');

                % Limiters
                hmax        = limiters(1)*max_h;
                dtdhmax     = limiters(2)*max_dtdh;

                cutoff  = H<hmax;
                limit   = dt<dtdhmax;
                omitmap = imdilate(~(cutoff.*limit),dilateMatrix('line',...
                    dilation,0));

                % Method unfiltered
                switch method
                    case 'NaN'
                        H(omitmap==1)=NaN;
                    otherwise
                        H(omitmap==1)=NaN;
                        
                        H = fillmissing(H,'movmean',30,2);
%                             H = fillmissing(H,'linear',2);
                        H(H<0)=0;
                end
            end

            % Apply smoothing
            H = smooth2a(H,smooth_zt(1),smooth_zt(2));
            H = fillmissing(H,'linear');
            % Remove 'islands'
            map = H>2e-2;
            regions = bwlabel(map);
            omit = (regions~=1); %#ok<NASGU> 
            omit = nan;
            
            % Mask
            if usemask
                try
                    H = id.maskAnalysis(H);
                catch
                    warning('Mask not detected.')
                end
            end

            % Flip Z s.t. z=0 is surface.
            if mean(diff(Z))>0
                Z = fliplr(Z);
            end


            if isempty(maxT)
                maxT = inf;
            end
            if numel(maxT)==1
                maxT = [0 maxT];
            end
            if ~isinf(maxT)
                map = within(T,maxT);
                T   = T(map);
                H   = H(:,map);
            end

            % Sampling
            if sampling~=1
                map     = roundto(T,sampling);
                H       = H(:,map);
                T       = T(map);
            end
            
            %% OUTPUTTING
            
            % Varargout CDATA
            switch lower(caxis_fmt)
                case 'thickness'
                    CDATA       = H;
                    CLIM        = [0 max(H,[],'all')];
                    CLABEL      = '$h$';
                case 'gradient'
                    [dt,~]      = gradient(H,T,Z);
                    CDATA       = dt;
                    CLIM        = std(dt,[],'all','omitnan').*[-1 1];
                    CLABEL   = '$\partial_t h$';
                case 'acceleration'
                    [dt,~]      = gradient(H,T,Z);
                    [dtt,~]     = gradient(dt,T,Z);
                    CDATA       = dtt;
                    CLIM        = std(dtt,[],'all','omitnan').*[-1 1];
                    CLABEL      = '$\partial_{tt} h$';
            end

            
            % Varargout
            switch nargout                    
                case 2
                    % Assume [H,CLABEL_FMT]
                    varargout{1} = CLABEL_FMT;
                case 3
                    % Assume [H,Z,T]
                    varargout{1} = Z;
                    varargout{2} = T;
                case 4
                    % Assume [H,Z,T,CSTRUCT]
                    varargout{1} = Z;
                    varargout{2} = T;
                    CSTRUCT = struct(...
                        'cdata',CDATA,...
                        'clim',CLIM,...
                        'label',CLABEL,...
                        'fmt',CLABEL_FMT,...
                        'map',map,...
                        'offset',offset);
                    varargout{3} = CSTRUCT;
                case 5
                    % Assume [H,Z,T,CSTRUCT,omit]
                    varargout{1} = Z;
                    varargout{2} = T;
                    CSTRUCT = struct(...
                        'cdata',CDATA,...
                        'clim',CLIM,...
                        'label',CLABEL,...
                        'fmt',CLABEL_FMT,...
                        'offset',offset);
                    varargout{3} = CSTRUCT;
                    varargout{4} = omit;
            end
            
            %% Functions
            function [t,CLABEL_FMT]=timeAxes(time_fmt,data)
                %% CONVERTS DATETIME ARRAY TO DD||HH||MM||SS DURATION ARRAY FROM EPOCH
                epoch = data.T(1);
                switch time_fmt
                    case 'Seconds'
                        t = seconds(data.T-epoch);
                        CLABEL_FMT = 's';
                    case 'Minutes'
                        t = minutes(data.T-epoch);
                        CLABEL_FMT = 'min';
                    case 'Hours'
                        t = hours(data.T-epoch);
                        CLABEL_FMT = 'hr';
                    case 'Days'
                        t = days(data.T-epoch);
                        CLABEL_FMT = 'd';
                end
            end
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch lower(items{m})
                        %% Name arguments
                        case {'seconds','minutes','hours','days'}
                            t_fmt  = items{m};                  
                        case {'thickness','gradient','acceleration'}
                            caxis_fmt = items{m};
                        case 'sample'
                            sampling = namevalue;
                        case 'mask'
                            usemask = true;
                        case {'maxt','maxtime'}
                            maxT = namevalue;
                        case {'parent','useparent'}
                            useparent = true;
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
        function varargout      = applyMetric(id,ref)
            if nargin<2
                ref = id.Image;
            end
            if isempty(ref)
                ref = id.Image;
            end
           % Draw lengths
           if isempty(id.Dimensions)
               id.getDimensions;
           end
            ds      = id.Dimensions.ds;
            [Z,X]   = size(ref);
            if isempty(ds)
                ds  = 1;
                x   = 1:X;
                z   = 1:Z;
            else
                ds  = id.Dimensions.ds;
                x   = linspace(0,X/ds,X);
                z   = linspace(0,Z/ds,Z);
            end
            
            shift = id.Offset;
            if isempty(shift)
                shift = 0;
            end
            
            z = fliplr(z);
            x = x-shift/ds;
            varargout{1} = x;
            varargout{2} = z;
            varargout{3} = ds;
        end
        function output         = parseOpticalFlow(id,varargin)
            timeseries = false;
            imageindex = 1;
            openApp    = false;
            loop       = [];
            useSubsample = false;
            usePadding = false;
            useSweep   = false;
            useSamples = false;
            parseInput(varargin)

            if openApp
                app = findall(0,'Type','Figure','Name','Optical Flow: Velocity Field');
                pos = [];
                if ~isempty(app)
                    pos = app.Position;
                    close(app)
                end
                app_velocityfield(id)
                if ~isempty(pos)
                    app = findall(0,'Type','Figure','Name','Optical Flow: Velocity Field');
                    app.Position = pos;
                end
                return
            end

            if isempty(id.OpticalFlow)
                warning('No OpticalFlow attributes exist.')
                return
            end
            
            ops = id.OpticalFlow;
            alg = icealgorithm(id);
            fprintf('Processing\n')

            if useSubsample
                useSubsample = 'usesubsample';
            end

            output = alg.getVelocity(...
                imageindex,...
                useSubsample,...
                'padding',usePadding, ...
                'loaded');

            
            if useSweep
                sweepParameters
                return
            end

            if ~timeseries
                plotting(nargout)
                return
            end
            if isempty(ops.tsStep)
                warning on
                warning('No timeseries step has been provided. Cancelling operation.')
                return
            end
            ops.ts = true;
            output = alg.VelocityTimeSeries(...
                'parallel',...
                'step',ops.tsStep,...
                'padding',usePadding,...
                'loaded');
            ops.ts = false;

            %% Reset
            

            %% Nested functions
            function plotting(nout)
                if nout>0
                    return
                end

                X = output.X;
                Z = output.Z;
                tiledlayout(2,2)
                nexttile;
                imagesc(X,Z,output.im1)
                addlabels('x','x (m)','y','z (m)','title','Input')
                addColorbar('cmap','gray')

                nexttile;
                imagesc(X,Z,output.U)
                hold on
                vis_flow(output.u,output.v,'X',X,'Z',Z)
                hold off
                addlabels('x','x (m)','y','z (m)','title','Velocity magnitude')
                addColorbar('cmap','thermal')

                nexttile;
                imagesc(X,Z,output.vor)
                addlabels('x','x (m)','y','z (m)','title','Vorticity')
                addColorbar('cmap','balance','pivot',0)

                nexttile;
                imagesc(X,Z,output.Q)
                addlabels('x','x (m)','y','z (m)','title','Second invariant')
                addColorbar('cmap','balance','pivot',0)
                clear output
            end
            function sweepParameters
                fprintf('Beginning sweep of parameters\n')
                alpha = [25 50 100];
                its   = [50 100 500];
                samples = 25; %#ok<NASGU>

                u = zeros([size(output.u),numel(alpha),numel(its)]);
                v= u;
                u = permute(u,[3 4 1 2]);
                v = permute(v,[3 4 1 2]);
                K = numel(its);

                if ~useSamples
                    if ~timeseries
                        q = parallel.pool.DataQueue;
                        parfevalOnAll(@clear, 0,'all');
                        parfevalOnAll(@warning, 0,'off','all');
                        afterEach(q,@displayProgressparfor)
                        w = waitbar(0,'Processing', ...
                            'Name','Parallel computation', ...
                            'Message','Preparing pool');
                        displayProgressparfor(w,numel(alpha))
                        parfor ii = 1:numel(alpha)
                            for jj=1:K
                                output = getVelocity(alg,imageindex,...
                                    'dt',ops.dt,...
                                    'method',ops.method,...
                                    'alpha',alpha(ii),...
                                    'iterations',its(jj),...
                                    'lambda',ops.lambda,...
                                    'filter',ops.filter,...
                                    'illumination',ops.illumination,...
                                    'resize',ops.resize,...
                                    'field',ops.field,...
                                    'rolling',ops.window,...
                                    'mask',ops.mask,...
                                    'alongslope',ops.alongslope,...
                                    'padding',usePadding); %#ok<PFBNS>

                                u(ii,jj,:,:) = output.u;
                                v(ii,jj,:,:) = output.v;
                                send(q,[])
                            end
                        end
                    else
                        if isempty(ops.tsStep)
                            warning('No timeseries step has been provided. Cancelling operation.')
                            return
                        end
                        for ii = 1:numel(alpha)
                            for jj=1:K
                                output = VelocityTimeSeries(alg,...
                                    'parallel',...
                                    'step',ops.tsStep,...
                                    'dt',ops.dt,...
                                    'method',ops.method,...
                                    'alpha',alpha(ii),...
                                    'iterations',its(jj),...
                                    'lambda',ops.lambda,...
                                    'filter',ops.filter,...
                                    'illumination',ops.illumination,...
                                    'resize',ops.resize,...
                                    'field',ops.field,...
                                    'rolling',ops.window,...
                                    'mask',ops.mask,...
                                    'alongslope',ops.alongslope,...
                                    'padding',usePadding);
                                u(ii,jj,:,:) = mean(output.u,3);
                                v(ii,jj,:,:) = mean(output.v,3);
                                fprintf('%i/%i\n',ii+jj-1,numel(alpha)*K)
                            end
                        end
                    end
                end
                u = permute(u,[3 4 1 2]);
                v = permute(v,[3 4 1 2]);
                assignin('base','u',u)
                assignin('base','v',v)
                assignin('base','out',output)
                fprintf('Plotting...\n')
                id.retrieveDetails;
                details = id.Details; %#ok<NASGU>
                %% Plotting
                tile=tiledlayout(numel(alpha),numel(its),'TileSpacing','compact','Padding','compact');
                for i=1:numel(alpha)
                    for j=1:numel(its)
                        nexttile
                        imagesc(output.X,output.Z,v(:,:,i,j)./(2*std(v(:,:,i,j),[],"all")))
                        caxis([-1 1])
                        cmocean('balance','pivot',0)
                        %                         xlim(XLIM)
                        %                         ylim(YLIM)
                    end
                end

                title_str = sprintf('%s=%i%s','$T_a$',25,' ($^\circ$C)');
                addlabels('ax',tile, ...
                    'x','Diffusion $(\alpha)\rightarrow$', ...
                    'y','$\leftarrow$ Iterations', ...
                    'title',title_str,...
                    'latex','fs',12)

                addlabels('ax',tile,'array','x',alpha,'y',its,'latex')
                c=colorbar;
                c.Title.String = '$U/2\sigma_U$';
                c.Title.Interpreter = 'latex';
                c.TickLabelInterpreter = 'latex';
            end
            %% Input parser
             function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch lower(items{m})
                        %% Name arguments
                        case {'ts','timeseries'}
                            timeseries = true;
                        case 'index'
                            imageindex = namevalue;
                        case {'openapp','open'}
                            openApp = true;
                        case {'padding','usepadding'}
                            usePadding = namevalue;
                        case 'subsample'
                            useSubsample = true;
                        case 'sweep'
                            useSweep = true;
                        case 'sweep-ts'
                            useSweep    = true;
                            useSamples  = true;
                        case 'loop'
                            loop = namevalue;
                            if numel(loop)==1
                                loop = 1:loop:numel(id.Folder.Files);
                                continue
                            end
                            loop(loop<1)=1;
                            loop(loop>numel(id.Folder.Files))=numel(id.Folder.Files);
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
        function parseAlgorithm(id,varargin)
            %% Parse algorithm
            % This subroutine applies the edge detection algorithm to a
            % all images from the time series.
            % -------------------------------------------------------------------------
            %  Parameters :
            % -------------------------------------------------------------------------
            %
            %  pool: [char] (Optional)
            %   call: ('pool'/'parallel'/'parpool')
            %   Runs the subroutine with a parallel pool.
            %   DEFAULT is 'off'.
            %
            %  step: [int/char] (Optional)
            %   call: ('step',1) OR ('dt','dt-dt_fmt')
            %   Sets the iteration steps either through image index count
            %   ('step') or time format ('dt').
            %   DEFAULT is 'step' = 1;
            % -------------------------------------------------------------------------
            % -------------------------------------------------------------------------
            %%
            wb          = [];
            poolstate   = false;
            step        = 1;
            parseInput(varargin)
            
            % Run constructor
            loop        = 1:floor(step):length(id.Folder.Files);
            
            
            % Construct analysis attributes
            % Draw lengths
            [Z,~]   = size(id.Image);
            ds      = 1;
            ds_fmt  = 'px';
            z       = 1:Z;
                        
            h = zeros(length(id.Edge),length(loop));
            try
                t = id.times(loop);
            end
            %%
            tic
            if poolstate
                q = parallel.pool.DataQueue;
                parfevalOnAll(@warning, 0,'off','all');
                if isempty(wb)
                    w = [];%waitbar(0,'Processing','Name','Serial computation');
                else
                    w = wb;
                    try
                        w.Title = 'Parallel computation';
                    end
                end
                LOOP = loop;
                if mean(diff(loop))>1
                   LOOP = 1:numel(loop);
                end
                parfor i=LOOP
                    alg = icealgorithm(id);
                    alg.loadImage(loop(i),true);
                    % Run sequencer
                    try
                        alg.Sequencer(loop(i));
                        h(:,i) = alg.iceedge;
                    catch ME
                        ME
                    end
                    send(q,[])
                end
                
            else
                
                alg = icealgorithm(id);
                j = 1;
                flag = NaN(1,length(loop));
                if isempty(wb)
                    w = [];%waitbar(0,'Processing','Name','Serial computation');
                else
                    w = wb;
                    try
                        w.Title = 'Serial computation';
                    end
                end
                displayProgressparfor(w,numel(loop))
                for i=loop
                    displayProgressparfor(w)
                    alg.loadImage(i,true);
                    % Run sequencer
                    try
                        alg.Sequencer(loop(i));
                        h(:,i) = alg.iceedge;
                    catch ME %#ok<NASGU> 
                        flag(j) = i;
                        j = j+1;
                    end
                end
                flag(isnan(flag)) = []; %#ok<NASGU> 
                %                 warning('Error at index: %.1f\nError%s\n',flag,ME.message)
            end
            toc
            %%
            id.Analysis = struct(...
                'H',h/ds,...             % Ice thickness array (px)
                'Z',z,...                % Depth vector (px)
                'T',t,...                % Absolute time vector
                'dt_fmt','',...
                'ds_fmt',ds_fmt...
                );           
            %% Input parser
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case {'pool','parallel','parpool'}
                            poolstate = true;
                        case {'step','dt'}
                            switch items{m}
                                case 'step'
                                    step = namevalue;
                                case 'dt'
                                    dt_fmt = split(items{m},'-'); %#ok<NASGU> 
                                    
                            end
                        case 'wb'
                            wb = namevalue;
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
        function loadOpticalFlow(id)
            files       = dir([strrep(fileparts(id.Path),'images','results') '\*.mat']);
            file2load   = files(contains({files.name},'OFtimeseries'));
            file        = load(fullfile(file2load.folder,file2load.name));
            file.u      = mean(file.u,3);
            file.v      = mean(file.v,3);

            calibrateResults;
            id.OpticalFlow.results = file;

            %% Nested functions
            function calibrateResults
                of = file;
                if ~isfield(id.OpticalFlow,'calibration')
                    warning('No calibration file detected.')
                    return
                end
                X = of.X;
                Z = of.Z;
                Vim = of.v;
                Uim = of.u;

                cb  = id.OpticalFlow.calibration;
                t   = cellfun(@str2double,{cb.roi.label})';
                pos = {cb.roi.pos};
                z_mean  = cellfun(@(y) y(:,2),pos,'UniformOutput',false);
                z_mean  = cell2mat(cellfun(@(x) 0.5*(x(1:end-1)+x(2:end)),z_mean,'UniformOutput',false)');
                x_mean  = cellfun(@(y) y(:,1),pos,'UniformOutput',false);
                x_mean  = cell2mat(cellfun(@(x) 0.5*(x(1:end-1)+x(2:end)),x_mean,'UniformOutput',false)');
                x_mean  = x_mean-min(X);

                

                Xp = arrayfun(@(x) find(X>=x,1,'first'),x_mean,'UniformOutput',false);
                Zp = arrayfun(@(x) find(Z<=x,1,'first'),z_mean,'UniformOutput',false);
                Xomit = ~cellfun(@isempty,Xp);
                Yomit = ~cellfun(@isempty,Zp);
                omit = logical(Yomit.*Xomit);
                Xp = cell2mat(Xp(omit));
                Zp = cell2mat(Zp(omit));
                v = cell2mat(cellfun(@(x) x(:,1),cb.vel,'UniformOutput',false)');
                u = cell2mat(cellfun(@(x) x(:,2),cb.vel,'UniformOutput',false)');
                U = hypot(u,v);
                U = U(omit);
                
                vim = zeros(size(Xp));
                uim = vim;
                for i=1:numel(Xp)
                    ZC = Zp(i)+[-2:2]; %#ok<*NBRAK1> 
                    XC = Xp(i)+[-2:2];
                    ZC = ZC(within(ZC,1:size(Vim,1)));
                    XC = XC(within(XC,1:size(Vim,2)));
                    vim(i) = mean(Vim(ZC,XC),'all');
                    uim(i) = mean(Vim(ZC,XC),'all');
                end
                UIM = hypot(vim,uim);
                
                Ucorr = U./UIM;
                [count,range]=histcounts(Ucorr,1:1e3:1e6);
                [~,idx] = max(count);
                scaling = mean(Ucorr(within(Ucorr,range(idx:idx+1))));
                file.scaling = scaling;
            end
        end
        function saveData(id,localpath,variable) %#ok<INUSD> 
            %% Check CD and ROOT
            if nargin<2;localpath = id.Path;end
            if nargin<3;variable = [];end%#ok<NASGU> 


            if isempty(localpath)||all(~localpath)
                warning on
                warning('Unspecified file path')
                return
            end
            savepath        = strcat(strrep(localpath,id.dirs.photos,id.dirs.results),'.mat');
            [path,expt,~]   = fileparts(savepath);
            
            if ~isfolder(path)
                mkdir(path);
            end

           
            data = id;

            % NEED TO HAVE SAVED "DATA" AS -STRUCT TO WORK
%             if ~isempty(variable)
%                 if ~any(cellfun(@(x) strcmp(x,variable),fieldnames(id)))
%                     warning on
%                     warning('The variable "%s" does not exist in "icedata". Variable was not saved.',variable)
%                     return
%                 end
%                 fprintf('Appending .mat file\n')
%                 m = matfile(fullfile(path,expt),'Writable',true);
%                 fprintf('\tLocally\n')
%                 save(fullfile(path,expt),sprintf('m.%s',variable),'-append')
%                 fprintf('\tCloud\n')
%                 save(fullfile(cloudpath,expt),sprintf('m.%s',variable),'-append')
%                 return
%             else
%             end

            % Save locally
            save(fullfile(path,expt),'data')
            fprintf('Saved to local drive\n')
            
            % Upload to Google Drive cloud (omitted)
            return
            cloudpath = strrep(path,[id.root '\'],id.cloudroot);

            if ~isfolder(cloudpath)
                mkdir(cloudpath);
            end
            
            save(fullfile(cloudpath,expt),'data')
            fprintf('Saved to cloud\n')
            assignin('base','data',data)
        end
        function retrieveDetails(id)
           %% RETRIEVES LOGGED DATA FROM GOOGLE DRIVE PATH
           root     = 'labview_logs'; 
           paths    = {'calibration','profile'};
           paths    = fullfile(cd,root,paths); 
           date_fmt = "yy_mm_dd";
           if isempty(id.times)
                id.times   = datetime({id.Folder.Files.date});
           end
           t0       = id.times(1);         % Initiation date
           date     = datestr(datetime(t0),date_fmt);
           id.Details.date = date;
           warning on
           profileReader                    % Attempts to load densimeter data (initial rho(z))
           calibrationReader                % Attempts to load stratifier calibration data
%            layerPredictions
            function calibrationReader
                %% READS IN THE STRATIFIER CALIBRATION PROFILE
                files    = dir(fullfile(paths{1},'\*.mat'));
                
                variables = {'t','s','vpp'};     %#ok<NASGU> % Variables listed in filename 
                
                % Filename date tag should match initiation date
                try 
                    %% Get filename variable attributes
                    file2read        = files(contains({files.name},date));
                    if isempty(file2read)
                        warning('No calibration profile exists for:\n%s\n',id.Path)
                        id.Details.calibration = [];
                        return
                    end
                    
                end
            end
            function profileReader
                %% READS IN THE LOGGED DENSITY PROFILE
                files    = dir(fullfile(paths{2},'\*.mat'));
                
                variables = {'t','Tf','Tj','s','vpp'};     % Variables listed in filename
                
                % Filename date tag should match initiation date
                try
                    %% Get filename variable attributes
                    file2read        = files(contains({files.name},date));
                    if isempty(file2read)
                        warning('No density profile exists for:\n%s\n',id.Path)
                        % Look for details in filename
                        try 
                            id.Path = strrep(id.Path,'psi','psu');
                            filedetails = split(id.Path,'\');
                            filedetails = extractAfter(filedetails{end-1},date);
                            temp        = str2double(extractBetween(filedetails,'_','dC'));
                            sal         = split(extractBetween(filedetails,'dC-','psu'),'-');
                            id.Details.temperature = temp(1);
                            id.Details.salinity    = str2double(sal');
                        end
                        id.Details.profile     = [];
                        return
                    end
                    
                    [~,filedetails,~]= fileparts(file2read.name);
                    details          = lower(split(extractAfter(filedetails,'#'),';'));
                    
                    readSingular     = @(x) str2double(x(~isletter(x)));
                    readRange        = @(x) str2double(split(x,'-'))';
                    
                    for i=1:numel(variables)
                        active  = variables{i};
                        idx     = contains(details,active);
                        if any(idx)
                            terms = split(details(idx),'_');
                            if isempty(terms)
                                break
                            end
                            
                            value = terms{2};
                            if numel(terms)>2
                                value = terms{3};
                            end
                            
                            switch active
                                case {'t'}
                                    id.Details.temperature = readSingular(value);
                                case 'Tf'
                                    [];
                                case {'Tj'}
                                    id.Details.forcing = readSingular(value);
                                case 's'
                                    id.Details.salinity = readRange(value);
                                case 'vpp'
                                    id.Details.vpp   = str2double(value);
                            end
                            
                        end
                    end
                    
                    %% Load in profile
                    filedata = load(fullfile(paths{2},file2read.name)).fluid;
                    id.Details.profile = filedata;
                    
                    % Use dimensions
                    [~,z]                       = id.applyMetric;
                    id.Details.profile.rhoz    = imresize(fliplr(id.Details.profile.rho'),size(z));
                    id.Details.profile.z       = z;
                    fit                         = polyfitn(z,id.Details.profile.rhoz,1);
                    id.Details.profile.rhofit  = polyval(fit.Coefficients,z);
                    id.Details.profile.drdz    = -fit.Coefficients(1);
                    id.Details.profile.drdzUn  = fit.ParameterStd(1);

                    tmpZ = id.Details.profile.z; %#ok<NASGU> 
                    tmpRho = id.Details.profile.rhofit; %#ok<NASGU> 
                    % Find 'ideal' density
                    rho_ideal = density( ...
                        linspace(id.Details.salinity(1),id.Details.salinity(end),numel(id.Details.profile.rhofit)), ...
                        id.Details.temperature);

                    % Find 'ideal' salinity
                    S_ideal = density_find('S',id.Details.temperature,...
                        [min(id.Details.profile.rhofit) max(id.Details.profile.rhofit)]);
                    
                    id.Details.profile.rho_ideal   = rho_ideal;
                    id.Details.profile.s_infer     = S_ideal;
                    id.Details.N = sqrt(abs(9.81/density(0,0)*id.Details.profile.drdz));
                catch ME
                    errorHandler(ME)
                end
            end
            function layerPredictions
                %% CALCULATES THE LAYER THICKNESS FROM PREDICTIVE SCALINGS
                %% Two prediction scalings
                %   Huppert & Turner, strong salinity gradient, depth
                %   averaged:
                %       h_HT = <0.66*diff([rho_int
                %       rho_inf])/d_zrho_inf>_z
                %
                %   Magorrian & Wells, weak salinity gradient, rest-plume,
                %   evaluated at initiation site z0
                %       h_MW = 2.83*diff([rho_int
                %       rho_inf])/d_zrho_inf|_z=z0
                %%
                try %#ok<*TRYNC> 
                    s       = id.Details.salinity;
                    prof    = id.Details.profile;

                    % Scalings
                    diffrho = @(rho0) abs((rho0-prof.rhofit)/prof.drdz);
                    h   = @(rho0) 0.66 * diffrho(rho0);
                    lrho= @(rho0) 2.83 * diffrho(rho0);


                    h1  = h(density(0,0));
                    plot(h1,prof.z)
                    h2  = h((density(linspace(s(1),s(2),numel(prof.rhofit)),0)));

                    h3  = lrho(density(0,0));
                    [];
                end
            end
        end
    end
    methods
        %% QUICK PLOTTING
        function seeHovmoller(id,ax)
            if nargin<2
                ax = gca;
            end
            [H,Z,t]    = id.parseAnalysis('Seconds', ...
                'sample',5*60, ...
                'maxtime',10*60*60, ...
                'mask', ...
                'parent');
            H(isnan(H))=0;
            H       = smooth2a(H,4,5);
            omit    = gradient(H)>1e-3;
            H(omit)=nan;
            contourf(ax,hours(seconds(t)),Z,1-H./H(:,1),'LevelStep',1/20)
            caxis([0 1])
            axis(ax,'ij')
            addColorbar('cmap','ice','reverse','colorbar')
        end
        function seePsi(id,layers)
            alg     = icealgorithm(id); 
            %mask    = alg.MaskIce(1);
            try
                mask = (arrayfun(@(x) alg.MaskIce(x),1:50:500,'UniformOutput',false));
                mask = ~logical(sum(cat(3,mask{:}),3)<numel(mask));
            catch
                f = figure('Visible','off');
                iceedge       = id.Edge;
                iceedge(1)    = 1;
                iceedge(end)  = 1;
                [x,z]         = id.applyMetric;
                mask      = poly2mask(iceedge,1:numel(z),numel(iceedge),numel(x));
                mask      = ~imdilate(mask,ones(3));
                delete(f)
            end
            shg,clf
            a(1) = axes; a(2) = axes;
            x = id.Velocity.x;
            z = fliplr(id.Velocity.z);
            [~,Cb] = contourf(a(1),x,z,id.Velocity.phix_t*1e6);%,'ShowText','on');
            Cb.LineColor = 'none';
            caxis(a(1),[-1 1]*5)%,round(caxis(a(1))))
            imagesc(a(2),x,z,id.Image,'AlphaData',~mask)
            a(2).Color = 'none';
            c(1)=addColorbar('ax',a(2),'cmap','gray','invisible');
            axis(a,'ij')
            a(2).XTick = [];
            a(2).YTick = [];
            c(2) = addColorbar('ax',a(1),'cmap','balance','pivot',0,...
                'title','CCW $$\leftarrow \Psi (\times10^{-6}) \rightarrow$$ CW',...
                'latex','fs',15,'levels',1+diff(caxis(a(1))));
%             set(c,'Location','southoutside')
%             c(2).TickLabels{1}     = sprintf('%s\n(%s)',c(2).TickLabels{1},'CCW');
%             c(2).TickLabels{end}   = sprintf('%s\n(%s)',c(2).TickLabels{end},'CW');
            id.retrieveDetails;
            temp    = id.Details.temperature;
            s       = id.Details.salinity;
            maxT    = round(minutes(max(id.Velocity.t)));
            title_str = sprintf('$T_a = %i%s C$, $S_a = [%i,%i]$ psu, t = [0,%i] min',temp,'^\circ',s(1),s(2),maxT);
            addlabels('ax',a(1),'latex','fs',15,...
                'x','Position (m)',...
                'y','Depth (m)',...
                'title',title_str)
            %% Attempt to show layers
            try
                arrayfun(@(x) line(gca,xlim,z(x).*[1 1],'LineStyle','--','Color','k'),layers)
            end
            %% Rendering
            if nargin<3
                return
            end
            latexformat('ratio',[1.5 2])
            lbl = extractBetween(id.Path,'regimes\','\ablate');
            filename = sprintf('%s%s_strmfcn.png','renders\ice\', ...
                lbl{1});
            warning off
            out2latex(filename,'r',300)
            disp('Image saved')
            warning on

        end
        function seeInput(id,idx)
           id.parseimage(gca,idx)
           c=addColorbar('latex','title','$I$ (gv)','fs',15,'cmap','gray','invert');
           addlabels('latex','x','$x$ (m)','y','$z$ (m)')
           c.Location = 'northoutside';
           c.Color = 'k';
        end
        function seeTurbulence(id)
%             [~,z,ds] = id.applyMetric;
%             T       = hours(id.times-id.times(1));
            plume   = id.Plume.turb;
            plume   = fillmissing(plume,'linear');
            [H,Z,T] = id.parseAnalysis('Hours','gradient','mask');
            H = smooth2a(H,50,2);
            imagesc(T,Z,log10(plume))
            CLIM = caxis;
            hold on
            [~,dHdZ] = gradient(H,T,Z);
            contour(T,Z,H,'LineColor',[1 1 1]*.5,'LevelStep',1e-2);
            contour(T,Z,dHdZ,'k','LevelList',0)
            legend({'h (1cm)','$\partial_z h=0$'},'Interpreter','latex')
            hold off
            caxis(CLIM)

%             shift = id.Plume.ROI0(1) - id.Plume.ROI(1);
%             in = id.Plume.loadImage(idx)-id.Plume.getMeanState('i',idx,'rolling',5);
%             imagesc((0:size(in,2)-1)/ds+shift/ds,z,imgaussfilt(in,1))
%             addlabels('latex','x','$x$ (m)')
%             warning off
%             caxis([-10 10])
%             c=addColorbar('latex','title','$\bar{I}-I$ (gv)','fs',15,'cmap','balance','pivot',0);
%             warning on
%             c.Location = 'northoutside';
%             c.Color = 'k';
        end
        function seePlume(id,idx)
           [x,z,ds] = id.applyMetric;
           shift = id.Plume.ROI0(1) - id.Plume.ROI(1);
           in = id.Plume.getResidualMeanState('i',idx);
           imagesc((0:size(in,2)-1)/ds+shift/ds,z,in)
           hold(gca,'on')
%            plot(gca,(id.Edge)/ds,z,'r','LineWidth',2)
           hold(gca,'off')
           addlabels('latex','x','$x$ (m)')
           warning off
           c=addColorbar('latex','title','$\overline{(|\bar{I}-I|)}$ (gv)','fs',15,'cmap','gray','invert');
           warning on
           c.Location = 'northoutside';
           c.Color = 'k';
           caxis([0 20])
%            drawrectangle('Position',id.Plume.ROI0/ds,'InteractionsAllowed','none','Color','w');
        end
        function seeBoth(id,idx)
           clf
           set(gcf,'Units','inches')
           set(gcf,'Position',[0    0.2292   11.2812   11.3854])
           tile = tiledlayout(1,3,'TileSpacing','compact');
           n(1) = nexttile;
           id.seeInput(idx)
           n(2) = nexttile;
           id.seeTurbulence(idx)
           n(3) = nexttile;
           id.seePlume(idx)
           set(n,'FontSize',15)
           t = datestr(id.times(idx)-id.times(1),'hh:mm:ss');
           time = sprintf('%s = %s','$t$',t);
           annotation('textbox',[0.2323,0.8502,0.1169,0.0271],...
               'String',...
               time,'HorizontalAlignment','center','FitBoxToText','off',...
               'Interpreter','latex','FontSize',15,'BackgroundColor','w')
%            addlabels('ax',tile,'title',,'latex','fs',20)
        end
        function seeAblation(id)
            %%
            [H,Z,T] = id.parseAnalysis('Hours');
            T       = T-T(1);
            dHdT    = gradient(imgaussfilt(H,[1 3]),T,Z);

            tiledlayout(10,10)
            n(1)=nexttile([7 7]);
            imagesc(T,Z,dHdT)
            caxis(2*std2(dHdT).*[-1 .1])
            colorbar
            colormap(1-cmocean('balance','pivot',0))
            n(2)=nexttile([3 7]);
            plot(T,mean(dHdT,1))
            n(3)=nexttile([7 3]);
            plot(mean(dHdT,2),Z)
            axis(n,'tight')
            axis(n(3),'ij')
        end
        function seeIceThickness(id,tmax)
            if nargin<2
                tmax = 5;
            end
            t0 = id.characteristicTime; %s
            t0 = t0/60^2;
            [h,z,t] = id.parseAnalysis('Hours', ...
                'maxtime',2*t0,...
                'mask','parent', ...
                'sample',t0/10);
            h = smooth2a(h,30,0);
            idx=find(t>t0,1);
            h = h./h(:,1);
            linedata=plot(h,z);
            axis ij
            addColorbar('lineplot', ...
                'colorbar', ...
                'cmap','haline', ...
                'linemap',t/t0)
            hold on
            plot(h(:,idx),z,'k','LineWidth',1)

            ylim([0 1])
%             waterfall(z,t,h',repmat(t,[numel(z) 1])')
%             shading flat
%             view(180,0)
            addlabels('y','$z$ (m)','x','$h/h_0$','latex')
            warning on
        end
        function seeAnalysis(id)
           if isempty(id.Analysis)
              warning('No analysis has been performed')
              return
           end
%            figure
           %%
           
           [H,Z,T]  = id.parseAnalysis('Hours');
           T        = T-T(1);
           truezero = 0;
           if ~isempty(id.offsetMask.vector)
               truezero = repmat(id.offsetMask.vector,[1 size(H,2)]);
           end
           time_intvl   = .2;
           relative     = double(H>truezero).*H;
           Thr          = [0 diff(abs(gradient(mod(T,time_intvl))))];
           Thr(1)       = max(Thr);
           Thrmax       = max(Thr);
           Thrstd       = std(Thr);
           range        = Thr>Thrmax-0.25*Thrstd;
           range(1)     = 1;
           h            = relative(:,range);
           nanMap       = h<3e-2;
           h(isnan(h))  = 0;
           t            = T(range);
           [dhdt,~]     = gradient(h,t,Z);
           
           relative     = h-h(:,1);
           relative(nanMap) = NaN;
           tiledlayout(1,5)
           nexttile([1 2])
           waterfall(Z,t,relative',repmat(t,[numel(Z) 1])')
           shading flat,view(180,0)
           addColorbar('cmap','rain','limits',[0 max(t)],'invert','title','Time (hours)')
           caxis([0 6])
           axis tight
           addlabels('y','Time (hours)','x','Depth (m)','z','Ice thickness')
           
           nexttile([1 2])
           surf(t,Z,h,dhdt)
           caxis(std2(dhdt).*[-1 1])
           zlim(std2(dhdt)*2.*[-1 1])
           shading flat
           view(0,90)
           addlabels('x','Time (hours)','y','Depth (m)')
           axis tight
           colorbar;
           set(gca,'Colormap',1-(cmocean('balance','pivot',0)))
           
           nexttile
%            dhdt(nanMap) = NaN;
           plot(mean(dhdt,2,'omitnan'),Z)
           axis tight
           
        end
        function seeOpticalFlow(id)
            of = id.OpticalFlow;
            if ~isfield(of,'results')
                id.loadOpticalFlow
            end
            r = of.results;
            U = hypot(r.u,r.v);
            Un = std(U,[],'all');
            tiledlayout(1,2)
            nexttile
            imagesc(r.X,r.Z,r.u/Un)
            addColorbar('cmap','balance', ...
                'pivot',0, ...
                'limits',[-3 3], ...
                'location','southoutside', ...
                'title','$||u||$','latex')
            nexttile
            imagesc(r.X,r.Z,r.v/Un)
            addColorbar('cmap','balance', ...
                'pivot',0, ...
                'limits',[-3 3], ...
                'location','southoutside', ...
                'title','$||v||$','latex')
        end
        function seeVelocity(id,idx,dt,varargin)
            %% Inputs   
            type    = 'none';
            resize  = [];
            repeat  = '';
            smoothing = 50;
            mask    = false;
            parseInput(varargin) 
            %%
            alg = icealgorithm(id);
            if ~isempty(id.Velocity)
                alg.velocity = id.Velocity;
            end
            id.Velocity = alg.showVelocity(idx,dt,type,repeat,...
                'resize',resize,'smoothing',smoothing,'mask',mask);
            %% Functions
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case {'velocity','vorticity','Q','none','vectors'}
                            type = items{m};
                        case {'repeat'}
                            repeat = 'repeat';
                        case {'resize','rescale'}
                            resize = namevalue;
                        case {'smooth','smoothing'}
                            smoothing = namevalue;
                        case 'mask'
                            mask = true;
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
        function [im,X,Z,t,idx,im0] = getAnomaly(id,idx,range)
            if nargin<2
                idx = 1;
            end
            if nargin<3
                range = 5;
            end
            T = id.times-id.times(1);
            switch class(idx)
                case 'duration'
                    idx = find(T>=idx,1,'first');
                case 'char'
                    idx = str2double(join(split(idx,'tau')));
                    etime = id.characteristicTime;
                    idx = find(T>=(seconds(etime*idx)),1,'first');
            end
            if isempty(idx)
                idx = numel(T);
            end
            t = T(idx);
            im = id.Plume.getMeanState('index',idx,'rolling',range,'anomaly');
            im = rescaleMax(im,true);
            [X,Z] = id.applyMetric(im);

            a   = icealgorithm(id);
            im0 = a.loadImage(idx,true);
            im0 = rescaleMax(im0,true);
        end
        function seeIceAndAnomaly(id,mode,varargin)
            if nargin<2
                mode = 'ice';
            end
            [im,X,Z,T,idx,im0] = id.getAnomaly(varargin{:});        
            omit = within(Z,[id.zoffset 1]);
            Z   = Z(omit);
            Z   = Z-min(Z);
            im  = im(omit,:);
            
            id.retrieveDetails
            DT = id.Details.temperature;
            ax1 = gca;
            rawimage = im0(omit,:);
            switch mode
                case 'ice'
                    imagesc(ax1,X,Z,rawimage)
                case 'anomaly'
                    imagesc(ax1,X,Z,im)
                    caxis(ax1,[-1 1]*.2)
                otherwise
                    mask = id.Plume.maskEdge(idx);
                    mask = imerode(mask,dilateMatrix('rect',5,20));
                    imagesc(ax1,X,Z,rawimage)
                    ax2 = axes('Position',ax1.Position);
                    imagesc(ax2,X,Z,im,'AlphaData',~logical(mask))
                    ax2.Color='none';
                    ax2.XTick = [];
                    ax2.YTick = [];
                    caxis(ax2,[-1 1]*.2)
            end
            xlim([0 .16])
            colormap gray
            addlabels(ax1,'x','$x$ (m)','z','$z$ (m)','latex','fs',15, ...
                'title',sprintf('%s=%i%s\nt=%s','$T_a$',DT,symbol('degClatex'),string(T,'hh:mm')))
            
        end
    end
    methods
        %% ANALYSIS
        function [etime,M,phi,t,ss,ss_un,h0,Mun] = characteristicTime(id)
            %%
            try
                [H,Z,t]    = id.parseAnalysis('Seconds','sample',5*60,'maxtime',10*60*60,'gradient','mask','parent');
            catch
                [H,Z,t]    = id.parseAnalysis('Seconds','sample',5*60,'maxtime',10*60*60,'gradient','mask');
            end
            Zomit       = ~within(Z,[0 1]);%(Z>offset&Z<1);
            H(Zomit,:)  = [];
            Z(Zomit)    = [];
            Z           = Z-min(Z);
            dz      = mean(diff(Z));


            h       = sum(fillmissing(H,'movmean',20,2),1,'omitnan')*dz;
            h0      = h(1);
            h       = h-h0;
            h0      = -h0;
            phi     = h./h0;

            etime_pt = 1/exp(1);
            t   = imresize(t,[1 1e3]);
            t   = t-min(t);
            phi = imresize(phi,[1 1e3]);
            etime_idx = find(phi>etime_pt,1,'first');
            etime_idx(arrayfun(@isempty,etime_idx)) = size(t,1);
            etime = t(etime_idx);
            if isempty(etime)
                etime = nan;
                fit = polyfitn(t(t<60^2),phi(t<60^2)*h0,1);
            else
                fit = polyfitn(t(t<etime),phi(t<etime)*h0,1);
            end

            % uncertainty
            etime_1_5 = find(phi>etime_pt,1,'first');
            if isempty(etime_1_5)
                etime_1_5=numel(t);
            end
            Mun_slope_std = (arrayfun(@(x)  polyfitn(t(1:x),phi(1:x)*h0,1).ParameterStd(1),2:etime_1_5,'ErrorHandler',@(~,~) nan));
            Mun_slope = (arrayfun(@(x)  polyfitn(t(1:x),phi(1:x)*h0,1).Coefficients(1),2:etime_1_5,'ErrorHandler',@(~,~) nan));

            Mun = std(Mun_slope,'omitnan');
            % gradient approach
            Mun_g = [mean(gradient(phi(t<etime),t(t<etime))),std(gradient(phi(t<etime),t(t<etime)))];
            

            M   = fit.Coefficients(1);
%             Mun = fit.ParameterStd(1);
            ss  = polyval(fit.Coefficients,t)/h0;
            ss_un = fit.ParameterStd(1);
            if ~nargout
                plot(t,phi,t,ss,'--k')
            end

        end
        function H = maskAnalysis(id,H)
             if ~isempty(id.analysisMask)
                omit     = id.analysisMask;
                H(omit)  = nan;
                return
             end
             warning off backtrace
             warning('No Analysis mask detected.')
             warning on backtrace
        end
        function multivariateRegression(id)
            %% MULTIVARIATE REGRESSION attempts to find the solution to
            %  m = c1*I^n*(S(z)-0)+c2*I^n*(Tinf-Tice)
            %% LOAD IN DATA
            [H,Z,T] = id.parseAnalysis('Hours');
            hfhrs   = roundto(T,.25);
            m       = gradient(H(:,hfhrs),T(hfhrs),Z);
            
            if ~isempty(id.analysisMask)
                omit     = id.analysisMask(:,hfhrs);
                m(omit)  = nan;
            end
            
            I       = movmean(id.Plume.turb,mean(diff(hfhrs)),2,'omitnan');
            scaling = id.Plume.scaling;            
            I       = I(:,hfhrs)*scaling;
            %% MASK 'I' LOOSELY
            lastidx = find(mean(~isnan(m)),1,'last');
            I(:,lastidx:end)=nan;
            if ~isempty(id.analysisMask)
                I(omit) = nan;
            end
            %% Regression
            id.retrieveDetails
            % Let m = c1*I^n*(S(z)-0)+c2*I^n*(Tinf-Tice)
            Tice    = -15;
            Tinf    = id.Details.temperature;
            Sz_exp  = linspace(id.Details.salinity(1),id.Details.salinity(2),numel(Z));
            rho     = id.Details.profile.rhofit;
            sreal   = [density_find('S',Tinf,min(rho)) density_find('S',Tinf,max(rho))];
            Sz_rea  = linspace(sreal(1),sreal(2),numel(Z));
            
            % constants
            c1      = linspace(0,1,20);
            c2      = -10:11;
            n       = linspace(0,5,20);
            sn      = linspace(0,1,20);
           
            
            merr = zeros([numel(c1) numel(sn) numel(n)]);
            for i=1:numel(c1)
                for j=1:numel(sn)
                    for k=1:numel(n)
                        tmp = c1(i).*I.^n(k).*Sz_rea.^sn(j)';%+c2(j).*I.^n(k)*(Tinf-Tice);
                        merr(i,j,k) = sqrt(sum(sum((m+tmp).^2,'omitnan')));
                    end
                end
                displayProgress('Progress',i,1,numel(c1))
            end
            
            [val,idx] = min(merr,[],3);
            %
            merr2 = merr(:,:,idx);
            imagesc(c1,sn,val)
            addlabels('x','c1','y','sn','title',sprintf('n=%.2f',mean(n,[1 2])))
            colorbar
            
            % Pseudo-inverse method?
            
            
            
%             c = pinv(I)*m(:,1);

           

            %% Plotting
            clf
            waterfall(Z,T(hfhrs),(m)',I'),view(-160,30)
            ylim([0 6])
            addColorbar('title','$I$ (m$^2$ s$^{-1}$)','latex','fs',20,'cmap','thermal')
            addlabels('x','Depth (m)','y','Time (hours)','z','$\dot{m}$ (m s$^{-1}$)','latex','fs',20)
            set(gca,'ColorScale','log')
            %%
%             clf,
%             imagesc(T(hfhrs),Z,m),colorbar,cmocean('balance','pivot',0),hold on,clim=caxis;
%             [~,cf]=contour(T(hfhrs),Z,log10(I),'LevelStep',.3,'ShowText','on','LabelSpacing',5e3);
%             caxis(clim)
%             xlim([0 6])
        end
        function getCorrelation(id)
            plume_analysis_test(id)
        end
        function logCorrelation(id,in)
            if isempty(id.Correlation)
                id.Correlation = in;
                return
            end
            id.Correlation(numel(id.Correlation)+1) = in;
        end
        function decomposeOFLayers(id,saving)
            %% DECOMPOSEOFLAYERS
            % Script used to convert velocity field time-series data into
            % layerwise components.
            % REQUIRES:
            %   1. OpticalFlow outputs:
            %       1.1. Timeseries output 
            %            - see PARSEOPTICALFLOW('timeseries')
            %       1.2. Layer table data
            %            - see PARSEOPTICALFLOW('open') -> Layer tab

            %%
            if nargin<2
                saving = false;
            end
            of = id.OpticalFlow;

            if ~isfield(of,'layers')
                error('No layer table detected')
            end
            try
                fprintf('Attempting to load OpticalFlow timeseries data\n')
                ofdata = dir([fileparts(strrep(id.Path,'images','results')) '\*timeseries.mat']);
                ofdata = load(fullfile(ofdata.folder,ofdata.name));
            catch
                error('No OpticalFlow timeseries data detected')
            end
            alg=icealgorithm(id);
            fprintf('Decomposing layers\n')
            alg.layerwise(ofdata,saving)
        end
        function [u0,v0,Z,X] = loadOFlayers(id)
            fname = [strrep(id.Path,'images','results') '_OFlayers.mat'];
            if ~isfile(fname)
                error('No layerwise OpticalFlow data detected')
            end
            o     = load(fname).out;
            out   = o.layers;
            zr    = cell2mat(arrayfun(@(x) x.Zrange,out,'UniformOutput',false)');
            zr    = cumsum(zr,2);

            X   = out(1).X;
            Z   = out(1).Z;
            u0  = o.original.u0;
            v0  = o.original.v0;
             % recompose
            uk = flipud(cell2mat(arrayfun(@(x) flipud(x.u0),out,'UniformOutput',false)'));
            vk = flipud(cell2mat(arrayfun(@(x) flipud(x.u0),out,'UniformOutput',false)'));

        end
        function decomposeMRLayers(id,xint)
            %% DECOMPOSEMRLAYERS
            % Script used to decompose melting rate data into layerwise
            % components
            % REQUIRES:
            %   1. Edge detection outputs:
            %       1.1 Analysis data
            %% Preamble
            a       = icealgorithm(id);
            % melting rate from dynamic ROI
            [m,z]   = a.dynamicROIMeltingRate;


            % melting rate from edge detection
            [M,Z]   = id.meltingrate;

            % find alignment matrix
            m2      = mean(m(z<=max(Z),:),2);
            z       = z(z<=max(Z));
            [b,c]   = xcorr(M,m2);
            [~,idx] = max(b);
            m       = circshift(m2,c(idx));


            of      = id.OpticalFlow;
            % guards
            if isempty(of)
                error('No OpticalFlow data detected')
                return
            end
            if ~isfield(of,'layers')
                error('No layer table detected')
                return
            end
            d = of.layers;
            d = d(contains(d.Type,'final'),:);

            %% Layerwise
            % split into layerwise coordinates
            pts = 1e3;
            etaZ = linspace(0,1,pts);
            zp0 = sort([0 d.z']);
            zp  = fliplr(arrayfun(@(x) find(z<=x,1,'first'),zp0));
            ZP  = fliplr(arrayfun(@(x) find(Z<=x,1,'first'),zp0));
            mr  = zeros(pts,numel(zp));
            MR  = mr;
            for k=1:numel(zp)-1
                mr(:,k) = -imresize(mean(m(zp(k):zp(k+1),:),2),[pts 1]);
                MR(:,k) = -imresize(mean(M(ZP(k):ZP(k+1),:),2),[pts 1]);
            end

            % load velocity fields
            fname = [strrep(id.Path,'images','results') '_OFlayers.mat'];
            if ~isfile(fname)
                error('No layerwise OpticalFlow data detected')
            end
            o     = load(fname).out;
            out   = o.layers;
            zr    = cell2mat(arrayfun(@(x) x.Zrange,out,'UniformOutput',false)');
            zr    = cumsum(zr,2);

            X   = out(1).X;
            Z   = out(1).Z;
            u0  = o.original.u0;
            v0  = o.original.v0;

            % recompose
            uk = flipud(cell2mat(arrayfun(@(x) flipud(x.u0),out,'UniformOutput',false)'));
            vk = flipud(cell2mat(arrayfun(@(x) flipud(x.u0),out,'UniformOutput',false)'));

            % streamfunctions
            if nargin<2
                xint = 2.5e-2; %cm
            end
            xint = X>xint;

            psik = flipud(cell2mat(arrayfun(@(x) flipud(flowfun(x.u0,x.v0*0)),out,'UniformOutput',false)'));
            psi0 = cumsum(u0(:,xint));%,v0(:,xint)*0);
            psi1 = flowfun(u0(:,xint)*0,v0(:,xint));
            %%
            TL=tiledlayout(1,4);
            n(1)=nexttile;
            imagesc(X,Z,u0/std(u0,[],'all'))
            n(2)=nexttile;
            imagesc(X,Z,v0/std(v0,[],'all'))
            n(3)=nexttile;
            imagesc(X(xint),Z,psi0/std(psi0,[],'all'))
            n(4)=nexttile;
            imagesc(X(xint),Z,psi1/std(psi1,[],'all'))
            addColorbar(TL,'cmap','balance','pivot',0, ...
                'limits',[-1 1],'levels',11,'latex','location','southoutside')
            addlabels(n([1 2]),'latex',...
                'c',{'$u/\sigma_u$','$w/\sigma_w$'})
            addlabels(n([3 4]),'title',...
                {sprintf('Residual %s = %.2f%s','$\Psi_u$',1e2*mean(psi0(end,:)/mean(abs(psi0),[1 2])),'\%'),...
                sprintf('Residual %s = %.2f%s','$\Psi_w$',1e2*mean(psi1(end,:)/mean(abs(psi1),[1 2])),'\%')},'latex',...
                'c','$\Psi/\sigma_\Psi$')
            [];
            %% Plotting (hidden)
% 
%             %%
%             tiledlayout(3,size(zr,1),'Padding','tight','TileSpacing','tight');
%             for k=1:size(zr,1)
%                 n(k)=nexttile(k,[2 1]);
%                 imagesc(out(k).X, ...
%                     out(k).etaZ, ...
%                     out(k).v/(2*std(out(k).v,[],'all')))
% %                 streamslice(out(k).u,out(k).v,10)
%                 caxis([-1 1])
%                 cmocean('balance','pivot',0)
%                 
%             end
%             for  k=1:size(zr,1)
%                 nn(k)=nexttile;
%                 plot(mr(:,k),etaZ,'k',...
%                     MR(:,k),etaZ,'r')
%                 r = cumsum(out(k).Zrange)/max(Z);
%                 title(sprintf('%.1f:%.1f%s',r(1),r(2),'%'))
% %                 xlim([0 inf])
%             end
%             L = legend({'Dynamic ROI','Analysis'},'Location','best');
%             title(L,'Melting rate source')
%             set([n nn],'XTick',[],'YTick',[])
        end
        function etaZ = getLayerPosition(id,time)
            noFile = true;
            etaZ = [];
            if ~isempty(id.OpticalFlow)
                if isfield(id.OpticalFlow,'layers')
                    noFile = false;
                end
            end
            if noFile
                return
            end

            if nargin<2
                time = 'last';
            end
            d = id.OpticalFlow.layers;
            
            layers  = max(d.Layer);
            init    = d(contains(d.Type,'initial'),:);
            final   = d(contains(d.Type,'final'),:);
            min_t   = min(init.t);
            max_t   = max(final.t);

            switch time
                case 'first'
                    t = min_t;
                case 'last'
                    t = max_t;
                otherwise
                    t = time;
            end

            % get current position for all times
            etaZ     = cell2mat(arrayfun(@(x) interp1( ...
                [init(x,:).t final(x,:).t], ...
                [init(x,:).z final(x,:).z],t), ...
                1:layers,'UniformOutput',false));
        end
        function [ax,M,U,z] = UvsM(id,plottype)
            [M,z,~,Mun]=id.meltingrate;
            etaZ = id.getLayerPosition('last');
            [~,Mav] = id.characteristicTime;
 
            alg       = icealgorithm(id);
            [m0,z0]   = alg.dynamicROIMeltingRate;
            m   = movmean(interp1(z0,-mean(m0,2),z)',200);
            mun = movmean(interp1(z0,-std(m0,[],2),z)',200);
            m(m<0)=0;
%             z = z+.02;
            id.retrieveDetails
            U       = nan;
            omitU   = true;
%             if id.Details.temperature>15
                try
                    omitU = false;
                    [U,Z,Uun]   =id.manualVelocity;
%                     keepZero    = U==0;
%                     keepZero    = logical(interp1(Z,double(keepZero),z));
                    U           = interp1nonunique(Z,U,z);
                    Uun         = interp1nonunique(Z,Uun,z);
                    if isrow(Uun), Uun = Uun';end

                end
%             end

            U = U*1e3;      % mm/s
            M = M*1e3*60;   % mm/min
            m = m*1e3*60;
            mun = mun*1e3*60+.05;
            Mav = Mav*1e3*60;
            Mun = Mun*1e3*60;
            
            omit = z>0;%>.2;
            z = z(omit);
            M = M(omit);
            Mun = Mun(omit);
            m = m(omit);
            mun = mun(omit);

            upwell = [];
            id.retrieveDetails
            if ~omitU
%                 Uun = abs(std(U-movmean(U,50),'omitnan'));
                U   = movmean(U,10)';
%                 U(keepZero) = 0;
                U   = U(omit);
                Uun = Uun(omit);
                Ubar = [U-Uun/2 U+Uun];
                Ubar(Ubar<=0)=0;

                try
                Uof = id.OpticalFlow.slice.v;%hypot(id.OpticalFlow.slice.v,id.OpticalFlow.slice.u);
                Zof = id.OpticalFlow.slice.Z;
                omit = within(Zof,[min(z) max(z)]);
                Uof = interp1(Zof(omit),Uof(omit),z);
                upwell = Uof>0;
                upwell = movmean(upwell,50)==1;
                end
            end
            z = z-min(z);
            switch id.Details.temperature
                case 10

                otherwise
                    M = wmean([M m],[.75 1],2);
            end
            Mbar = [M-Mun M+Mun];
            Mbar(Mbar<=0)=0;
            mbar = [m-mun m+mun];
            mbar(mbar<=0)=0;

            if nargin<2
                plottype = 'both';
            end

            ax = gca;
            switch plottype
                case 'z'
                    ax=plotZ(ax);
                case 'scatter'
                    plotScatter(ax)
                case 'both'
                    tiledlayout(2,1,'Padding','loose','TileSpacing','compact')
                    set(gcf,'Color','w')
                    n(1)=nexttile;
                    n(1).Visible = 'off';
                    n(2)=nexttile;
                    plotScatter(n(2))
                    plotZ(n(1))
                case 'U'
                    plot(U,z,'k')
                    p1=shadedErrorBarPoly(z,Ubar(:,1),Ubar(:,2),'y','k');
                    if id.Details.temperature<20
                        p1(2)=shadedErrorBarPoly(z,U.*upwell',U*0,'y','b');
                    end
                    addlabels('x','$U$ (mm s$^{-1}$)','latex','fs',12.5,...
                    'title',sprintf('%s=%i%s','$T_a$',id.Details.temperature,'$^\circ$C'))
                    axis ij
                    xlim([0 max(U)*(1.5)])
                    ylim([0 1.1])
                case 'M'
                    hold on
                    plot(M,z,'r')
                    p2(1)=shadedErrorBarPoly(z,Mbar(:,1),Mbar(:,2),'y','r');
                    plot(mean(M,'omitnan')*[1 1],Z([1 end]),'--r')
                    addlabels('x','$\dot{m}$ (mm min$^{-1}$)','y','$z$ (m)', ...
                    'latex','fs',12.5)
                    axis ij
                    xlim([0 max(M)*(1.5)])
                    ylim([0 1.1])
                otherwise
                    return
            end
            if ~nargout
                clear ax
            end
            %%
            function plotScatter(ax)
                fit     = polyfitn(sort(U),sort(M),1);
                xtmp    = linspace(0,max(U));
                scatter(ax,U,M,1,z,'.')
                addColorbar(ax, ...
                    'cmap','thermal', ...
                    'title','$z$ (m)', ...
                    'latex', ...
                    'levels',9, ...
                    'limits',[0 0.9])
                hold(ax,'on')
                plot(ax,xtmp,xtmp*fit.Coefficients(1),'--k')%polyval(fit.Coefficients,xtmp),'--k')
                hold(ax,'off')
                xlim([0 inf])
                ylim(xlim)
                addlabels('x','$U$ (mm/s)','y','$\dot{m}$ (mm/min)', ...
                    'title',sprintf('%s=%i%s','$\Delta T_a$',id.Details.temperature,'($^\circ$C)'),...
                    'latex','fs',12.5)
            end
            function ax=plotZ(ax)
                set(ax,'Visible','off')
                a = xxaxis(M,z,U,z,ax);
                p1 = [];
                hold(a,'on')
                p2(1)=shadedErrorBarPoly(z,Mbar(:,1),Mbar(:,2),'y','r',a(1));
%                 plot(a(1),mean(M,'omitnan')*[1 1],Z([1 end]),':r')
%                 plot(a(2),mean(U,'omitnan')*[1 1],z([1 end]),':k')
                if ~omitU
                    p1=shadedErrorBarPoly(z,Ubar(:,1),Ubar(:,2),'y','k',a(2));
                    if id.Details.temperature<20
                        p1(2)=shadedErrorBarPoly(z,U.*upwell',U*0,'y','b',a(2));
                    end
                end
                hold(a,'off')
                set([p1 p2],'EdgeAlpha',0)
                
                axis(a,'ij')
                addlabels(a(2),'x','$U$ (mm s$^{-1}$)','latex','fs',12.5,...
                    'title',sprintf('%s=%i%s','$T_a$',id.Details.temperature,'$^\circ$C'))
                addlabels(a(1),'x','$\dot{m}$ (mm min$^{-1}$)','y','$z$ (m)', ...
                    'latex','fs',12.5)
                a(2).YTick = [];
                ylim(a,[0 1.1])
                pause(1-8)
                a(1).Tag = sprintf('m-%i',id.Details.temperature);
                a(2).Tag = sprintf('u-%i',id.Details.temperature);
                ax = a;

                % rescale xlim s.t. mean(U) coincides with mean(M)
                %%
                if ~omitU
                    xlim(a(2),[0 mean(U,'omitnan')*4])
                end
                xlim(a(1),[0 mean(M,'omitnan')*4])
                line(a(2),mean(U,'omitnan').*[1 1],ylim(a(2)),'LineStyle',':','Color','k')
            end
        end
        function [U,Z,Uun] = manualVelocity(id)

            [~,z0]      = id.applyMetric;
            mv          = id.OpticalFlow.calibration;
            t           = id.times(cellfun(@str2double,{mv.roi.label}))-id.times(1);
            [t,map]     = sort(t);
            mv.vel      = cellfun(@(x) hypot(x(:,1),x(:,2)),mv.vel,'UniformOutput',false);
%             valid       = cellfun(@uniqueIDX,mv.z,'UniformOutput',false);
            
            U           = cell2mat( ...
                arrayfun(@(idx) interp1(mv.z{idx},mv.vel{idx},z0),1:numel(mv.z),'UniformOutput',false) ...
                );
            
            U = U(:,map);
            UUn = mean(U,2,'omitnan')+[-1 1].*std(U,[],2,'omitnan');
%             U           = cat(1,mv.vel{:});
%             Z           = cat(1,mv.z{:});
%             [~,idx]     = unique(Z);
%             Z(idx)      = Z(idx)+1e-5;
%             [Z,sorting] = sort(Z);
%             
%             U   = U(sorting,:);
            U   = movmean(mean(U,2,'omitnan'),10);
%             U   = interp1nonunique(Z,U,z0);
            try
                Uun = id.OpticalFlow.uncertainty.U/2;
                Zun = id.OpticalFlow.uncertainty.Z;
                Uun = interp1(Zun,Uun,z0);
            catch
                Uun = std(U,'omitnan')*ones(size(U))*1e3;
            end
            Z   = z0;
            function idx = uniqueIDX(in)
                [~,idx] = unique(in);
            end
            function out = errFcn(~,varargin)
                out = nan(size(z0));
            end
        end
        function [M,Z,T,Mun] = meltingrate(id,trange,samplingrate)
            if nargin<2
                trange = [0 60];    % first hour
            end
            if nargin<3
                samplingrate = 6;   % every 6 minutes
            end
            [H,Z,T,C]   = id.parseAnalysis( ...
                'Minutes',...
                'gradient',...
                'sample',samplingrate,...
                'maxT',trange,...
                'parent');
            T       = T*60;

%             try
%                 M = layerFollowing(id,C.cdata,T);
%             end

            M       = -C.cdata/60;
%             Mun     = std(M,[],2);
            MAV     = mean(M,2); % m/s
            Mun        = mean(M,2,'omitnan');
            q           = [1/4 3/4]; % quantile interval
            for k=1:numel(q)-1
                Mun(:,k)   = std( ...
                    quantile(M,linspace(q(k),q(k+1), ...
                    round(numel(T)*diff([q(k) q(k+1)]))),2), ...
                    [],2); %#ok<AGROW>
            end

            MAV(MAV<0) = 0; % discard regrowth
            M = MAV;
            function M = layerFollowing(id,M,T)
                %% Preamble
                of = id.OpticalFlow;
                if isempty(M)
                    return
                end
                if ~isfield(of,'layers')
                    return
                end
                d = of.layers;
                %% Compare
                a         = icealgorithm(id);
                % melting rate from dynamic ROI
                [m0,z0]   = a.dynamicROIMeltingRate;

                %% Layer tracking
                layers  = max(d.Layer);
                init    = d(contains(d.Type,'initial'),:);
                final   = d(contains(d.Type,'final'),:);
                min_t   = min(init.t);
                max_t   = max(final.t);
                t = T(within(T,min_t:max_t));
                % get current position for all times
                pos     = arrayfun(@(x) interp1( ...
                    [init(x,:).t final(x,:).t], ...
                    [init(x,:).z final(x,:).z],t), ...
                    1:layers,'UniformOutput',false);
                % convert to pixel coordinates
                pos = cellfun(@(y) arrayfun(@(x) find(Z<x,1,'first'),y),pos,'UniformOutput',false);
                posf= cellfun(@(x) x(end),pos);
                %%
                Mav = zeros(layers-1,numel(t));
                Zav = zeros(layers-1,1);
                pts = 1e3;
                m   = zeros(pts,layers-1,numel(t));
                mf  = zeros(pts,layers-1);
                etaZ = linspace(0,1,pts);
                for i=1:layers-1
                    posi  = pos{i};
                    posi2 = pos{i+1};
                    mf(:,i) = imresize(mean(M(posf(i+1):posf(i),:),2),[pts 1]);
                    for j=1:numel(t)
                        range   = sort([posi(j) posi2(j)]);
                        range   = range(1):range(2);
                        mtmp    = M(range,j);
                        mtmp(mtmp>0)=0;
                        Mav(i,j) = mean(mtmp);
                        m(:,i,j) = imresize(mtmp,[pts 1]);
                    end
                    Zav(i) = Z(round(mean([posi posi2])));
                end
            [];
            end

        end
        function tbl = OFlayers(id)
            DRHO        = @(S,T) density(S,T*0)-density(S,T);
            eta         = @(S,T,drhodz) 0.66*(DRHO(S,T))./drhodz;
            
            d           = id.OpticalFlow.layers;
            init        = d(contains(d.Type,'initial'),:);
            final       = d(contains(d.Type,'final'),:);

            id.retrieveDetails;
            D           = id.Details;
            eta_pred    = eta(mean(D.salinity),D.temperature,-D.profile.drdz);
            H           = 0.66; Hun = 0.06;
            Drho        = DRHO(mean(D.salinity),D.temperature);
            DrhoUn      = mean(arrayfun(@(x) DRHO(mean(D.profile.s_infer),D.temperature+2*x),-1:1)-Drho);
            DrDz        = D.profile.drdz;
            DrDzUn      = D.profile.drdzUn;

            pred_un     = eta_pred*sqrt((Hun/H)^2+(DrhoUn/Drho)^2+(DrDzUn/DrDz)^2);
            eta_obs     = abs(mean([mean(diff(init.z)) mean(diff(final.z))]));
            obs_un      = mean([std(diff(init.z)) std(diff(final.z))]);
            tbl = table(D.temperature, ...
                max(d.Layer)-1,...
                eta_obs*1e3, ...
                obs_un*1e3, ...
                eta_pred*1e3, ...
                pred_un*1e3,...
                eta_obs/eta_pred,...
                eta_obs/eta_pred*sqrt((obs_un/eta_obs)^2+(pred_un/eta_pred)^2),...
                'VariableNames', ...
                {'DTa(degC)','eta_count','eta_obs(mm)','obs_un(mm)','eta_pred(mm)','pred_un(mm)','eta_obs/eta_pred','obs_un/pred_un'});


        end
        function theoryEstimates(id,eta)
            id.retrieveDetails
            d   = id.Details;
            p   = d.profile;
            Ta  = d.temperature;
            Z   = p.z;
            Sa      = linspace(d.salinity(1),d.salinity(2),numel(Z));
            Sa_inf  = linspace(p.s_infer(1),p.s_infer(2),numel(Z));
            DTa = Ta-0;
            STt = property('StT');
            LL  = property('LL');
            g   = 9.81;
            alpha = property('water','alpha','T',DTa);
            N   = d.N;
%             b   = g*alpha*DTa/eta;
            
            etaS        = calculateLayers(Ta,Sa,Z,[],true);
            etaS_inf    = calculateLayers(Ta,Sa_inf,Z,[],true);
            plot(etaS,Z,'k',...
                etaS_inf,Z,'r')
%             delT = (mean(id.meltingrate)/(-STt/LL*sqrt(eta*g*alpha)^.5))^1.5; % estimated temperature loss over layer depth
[];

        end
    end
    methods
        %% IMAGE SAVING and MOVIE
        function saveImages(id,stride,interval_min,saving)
            %% SAVEIMAGES saves a sequence of images to the class
            % image = image(x,z,t)
            % where t = [interval_1:interval_1+stride; 
            %            interval_2:interval_2+stride
            %               ...
            %            interval_n:interval_n+stride]
            %%
            if nargin<4
                saving = false;
            end
            if stride<1
                error('"Stride" must be greater than 1.')
            end
            T       = minutes(id.times-id.times(1));
            T(T>30) = [];
            loop    = roundto(T,interval_min)-stride;
            loop(1) = 1;
            
            alg     = icealgorithm(id);
            alg.loadImage(1);
            %%
            clc
            images  = nan([stride*numel(loop) size(alg.output)]);
            Tlist   = nan(1,stride*numel(loop));

            fprintf('Images to be saved:\t%i\t(stride = %i, interval = %i min)',stride*numel(loop),stride,interval_min)
            
            k = 1;
            for i=1:numel(loop)
                for j=0:stride-1
                    images(k,:,:) = alg.loadImage(j+loop(i));
                    Tlist(k)      = T(j+loop(i)); 
                    k = k+1;
                end
                displayProgress('Loading images',i,1,numel(loop),'delete')
            end
            id.Images = struct;
            id.Images.im = permute(images,[2 3 1]);
            id.Images.T = Tlist;

            if saving
            id.saveData;
            end
        end
        function filename = createMovie(id,filename,makegif)
            if nargin<2
                error('Define a filename')
            end
            if nargin<3
                makegif = false;
            end
            id.retrieveDetails;
            lblhov  = sprintf('$%s=%i%s$C','\Delta T_a',id.Details.temperature,'^\circ');%'Hovmoller';
            lblim   = '';
            warning off
            offset      = 0.05;
            [H,Z,Tt,C]  = id.parseAnalysis('Minutes','sample',3,'maxtime',5*60,'gradient','mask');
            Zomit       = ~(Z>offset&Z<1);%plot(Zomit,Z)
            H(Zomit,:)  = [];
            Z(Zomit)    = [];
            Z   = Z-min(Z);
            dz  = mean(diff(Z));
            T0  = hours(id.times-id.times(1));
            TMP = smooth2a((H-H(:,1))./H(:,1),10,5);
            dt  = gradient(TMP,Tt/60,Z);
            TMP = -TMP;
            TMP(isnan(TMP))=1;

            fig = figure(1);
            set(gcf,'Color','w')
            tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
            n(1)=nexttile;
            n(2)=nexttile;

            fs = 15;
            [im,x,z]=id.parseimage(1);
            imagesc(n(2),x,z-offset,im),colormap gray,caxis([0 255])
            addlabels('ax',n(2), ...
                'x','$x$ (m)','y','$z$ (m)','title',lblim,'latex', ...
                'fs',fs)
            contourf(n(1),Tt/60,Z,TMP,'LevelStep',.05,'Color','k','EdgeAlpha',.5);
            addColorbar('ax',n(1),'cmap','ice','reverse','levels',20, ...
                'title','$1-h(z,t)/h_0$','latex','fs',15,'colorbar')
            caxis(n(1),[0 1])
            axis(n,'ij')
            addlabels('ax',n(1),'x','$t$ (hours)','y','$z$ (m)','latex','title',lblhov,'fs',fs)

            filename0 = filename;
            filename = [fullfile('..\Data\Movies',filename) '.avi'];
            V = VideoWriter(filename);
            V.Quality = 100;
            open(V)
            for i=C.map'
                figure(2)
                [im,~,~,time]=id.parseimage(i);
                figure(1)
                n(2).Children.CData = im;
                ylim(n(2),[0 max(Z)])
                caxis(n(2),[0 255])
                delete(findall(n(1),'Tag','line'))
                line(hours(time).*[1 1],ylim(n(1)),'Parent',n(1),'Tag','line','Color','r')
                writeVideo(V,getframe(gcf));
                pause(.1)
            end
            close(V)
            if makegif
                fprintf('Generating gif:\n')
                avi2gif('filename',filename, ...
                   'destination',['renders\ice\' filename0 '.gif'], ...
                   'delay0',.5, ...
                   'delayend',.5, ...
                   'fps',.05, ...
                   'forever')
                fprintf('Gif saved to \n\t%s\n',['renders\ice\' filename0 '.gif'])
            end
        end
    end
    methods (Hidden)
        %% FOLDER/DIRECTORY FUNCTIONS
        function varargout  = checkDirectory(id)
            path = split(id.root,'\');
            root = split(pwd,'\'); %#ok<*PROP>
            match = matches(root,path);
            
            if all(match)
                varargout{1} = id.dirs.photos;
                varargout{2} = id.dirs.results;
                return
            end
            
            base = id.root;
            
            id.dirs.photos     = fullfile(base,id.dirs.photos);
            id.dirs.results    = fullfile(base,id.dirs.results);
            
            varargout{1} = id.dirs.photos;
            varargout{2} = id.dirs.results;
        end
        function folder_obj = folderSearch(~,exptFolderName)
            % Find experiment folder under root
            folders = split(genpath('.'),';');
            list    = contains(folders,exptFolderName);
            list2   = contains(folders,'ablation');
            idx     = find(list.*list2);
            path    = folders{idx}; %#ok<FNDSB>
            
            folders = Folder('.'); %#ok<CPROPLC>
            parent  = folders;
            state   = true;
            clc
            while state
                folderSearchNest
            end
            folder_obj = parent;
            
            function folderSearchNest
                for i=1:length(parent.Children)
                    children = parent.Children{i};
                    if contains(path,children.Path)
                        parent = children;
                        folderSearchNest;
                        if length(parent.Children)==0 %#ok<ISMT>
                            state = false;
                            return;
                        end
                    end
                end
            end
        end
        function DT         = meanTime(id,fmt)
            %% MEANTIME returns the average time between images
            % whilst accounting for unsually large differences (e.g.,
            % images are 1s but there's a 5min pause somewhere)
            
            %%
            if nargin<2
                fmt = 'ms';
            end
                        
            date    = sort(id.times);
            date    = date-date(1);
            cent    = mean(diff(date));
            devi    = std(diff(date));
            outlier = max(diff(date));
            DT      = milliseconds(cent);
            
            
            
            if milliseconds(cent)<1e3
                DT = round(milliseconds(cent),1);
                DT = dtConvert(milliseconds(DT),fmt);
                return
            end
            
            %%
            if outlier>2*devi
                [~,locs]      = findpeaks(diff(seconds(date)));
                locs(end+1)   = 1;
                locs(end+1)   = numel(date);
                locs          = sort(locs);
                ranges        = numel(locs);
                
                dt_av = repmat(date(1),[1 numel(ranges)]);
                for j=1:ranges-1
                    range      = locs(j)+1:locs(j+1)-1;
                    dt_av(j)   = mean(diff(date(range)));
                end
                DT = round(milliseconds(mean(dt_av)),1);
                DT = dtConvert(milliseconds(DT),fmt);
%                 if nargout>0
%                     varargout{1} = mean(dt_av);
%                 end
            end
            
            function dtout = dtConvert(dt,fmt)
                switch fmt
                    case {'ms','milliseconds'}
                        dtout = milliseconds(dt);
                    case {'s','seconds'}
                        dtout = seconds(dt);
                    case {'m','min','minutes'}
                        dtout = minutes(dt);
                    case {'h','hr','hours'}
                        dtout = hours(dt);
                end
            end
        end
        function renameDirectory(id,newname,varargin)
           %
           appendstr    = {};
           treehandle   = [];
           parseInput(varargin)
           
           % Make newfolder under same path
           path     = fileparts(id.Path);
           oldname  = id.Expt;
           renameFolder(path,oldname,newname)
           
           % Rename results (local and cloud versions)
           path_results  = strrep(path,id.dirs.photos,id.dirs.results);
           renameFolder(path_results,oldname,newname)
           path_cloud    = strrep(path_results,[id.root '\'],id.cloudroot);
           renameFolder(path_cloud,oldname,newname)
           
           return
           %% Nested experiments
           % Nested experiments must have their id.Path and id.Folder
           % updated?
           
           list =  split(genpath(fullfile(path_results,newname)),';');
           
           
           
           tmp  = []
           tmp2 = []
           tmp3 = []
           tmp4 = []
           
           
           % Change attributes
           id.Path = newfolder;
           id.Expt = newname;
           
           % Save again in new folders
           id.saveData;
           
           return
           
           

           %% Input parser
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case isa(items{m},'matlab.ui.container.Tree')
                            treehandle = items{m};
                        case 'append'
                            appendstr = namevalue;
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
        function loadExperiment(id,exptFolderName)
            %% Loads the experiment data
            % Assumes the experiment directory is located somewhere in the
            % working directory.
            
            folder_obj = id.folderSearch(exptFolderName);
            readFolder(id,folder_obj)
            
            % Attempt to load saved filed
            loadpath = strrep(path,'images','results');
            file = what(loadpath);
            if ~isempty(file)
                load(file) %#ok<LOAD>
            else
                loadControls(id,[],[])
                roiRequest(id,'',1)
            end
            getDimensions(id)
            checkForUpdatedImages(id)
            
        end
        function loadControls(id,values,viewing)
            if isempty(values)&&isempty(viewing)
                values = struct(...               % Default values
                    'image',...
                    struct(...
                    'vignette',0,...
                    'blurx',0,...
                    'blury',0,...
                    'clim',[0 255]),...
                    'edge',...
                    struct(...
                    'algorithm','Sobel',...
                    'threshold',1e-3,...
                    'sigma',1e-3,...
                    'gradient',0),...
                    'dilation',...
                    struct(...
                    'shape','Line',...
                    'arg1',0,...
                    'arg2',0,...
                    'medfilt',1,...
                    'reject',0),...
                    'connect',...
                    struct(...
                    'connect',0,...
                    'smooth',1e-3),...
                    'ROI',...
                    struct(...
                    'set',false,...
                    'reset',false,...
                    'rotation',0,...
                    'flipUD',false,...
                    'flipLR',false));
                
                viewing = struct(...                % Default viewing
                    'image',true,...
                    'edge',false,...
                    'dilation',false,...
                    'connectivity',false,...
                    'overlay',false,...
                    'time',false,...
                    'timeRange',0);
            end
            id.Controls = struct('values',values,'viewing',viewing);
        end
        function readFolder(id,folder)
            switch class(folder)
                case 'Folder'
                    % Reads the FOLDER class object
                    id.Expt    = folder.Label;
                    id.Path    = folder.Path;
                    id.Folder  = folder;
                    if ~isempty(folder.Children)
                        id.times   = datetime({folder.Files.date});
                    end
                case 'matlab.ui.container.TreeNode'
                    readFolder(id,Folder(folder.UserData)) %#ok<CPROPLC> 
                    return
            end
            if ~isempty(folder.Files)
                id.times = datetime({folder.Files.date});
                [~,order] = sort(id.times);
                id.times = id.times(order);
                id.Folder.Files = id.Folder.Files(order);
            end
            getDimensions(id)
        end
        %% DATA FUCNTIONS
        function getDimensions(id)
            if ~isempty(id.Folder.Files)
                % Check if file is valid
                file        = id.Folder.Files(1);
                [~,~,ext]   = fileparts(fullfile(id.Folder.Path,file.name));
                if ~contains(ext,{'.tiff','.bmp'})
                    return
                end
                
                dates   = datetime({id.Folder.Files.date});
                dt      = seconds(mean(diff(dates)));
                
                default_dims = struct('ds',1,'ds_fmt','m','dt',dt,'dt_fmt','seconds');
                
                if isempty(id.Dimensions)
                    id.Dimensions = default_dims;
                else
                    if isempty(id.Dimensions.ds)
                        id.Dimensions = default_dims;
                    else
                        id.Dimensions.dt = dt;
                    end
                            
                end
                
            end
        end
        function details(id)
            %% BUILDS A STRUCTURE WITH DETAILS OF THE EXPERIMENT
            details             = struct;
            details.date        = [];
            details.temperature = [];
            details.salinity    = [];
            details.calibration = [];
            details.profile     = [];
            details.vpp         = [];
            details.julabo      = [];
            id.Details         = details;
        end
        function getJulabo(id,path)
            if nargin<2
                path = id.dirs.julabo;
                if isempty(path)
                    path = id.Path;
                end
            end
            [temp,time,~,id.dirs.julabo]=scanTemp(path);
            id.Details.julabo = struct('temp',temp,'time',time);
            id.saveData;
        end
        function getZero(id,ax)
            im  = [];
            f   = [];
            makeFigure;
            [z0,x0]=size(im);
            
            hasPoly = ~isempty(id.offsetMask.poly);
            pos0 = [x0/2 0; x0/2 z0];
            if hasPoly
                pos0 = id.offsetMask.poly;
            end
            line = drawpolyline(ax,'Position',pos0);
            addlistener(line,'ROIMoved',@lineROI);
            addlistener(line,'MovingROI',@lineROI);
            POS = line.Position;
            x = POS(:,1);
            y = POS(:,2);
            function lineROI(~,~)
                pos = line.Position;
                line.Position(1,2)      = 1;
                line.Position(end,2)    = z0;
                x = pos(:,1);
                y = pos(:,2);
            end
            function makeFigure
                f   = uifigure('Name','Drawing Mode: True zero');
                gl  = uigridlayout(f,'RowHeight',{'10x','1x'},'ColumnWidth',{'1x',80,80,'1x'});
                ax  = uiaxes(gl);
                ax.Layout.Column = [1 4];
                
                acceptButton = uibutton(gl,'Text','Accept','ButtonPushedFcn',@accept);
                acceptButton.Layout.Row = 2;
                acceptButton.Layout.Column = 2;
                
                useEdgeButton = uibutton(gl,'Text','Use Edge','ButtonPushedFcn',@useEdge);
                useEdgeButton.Layout.Row = 2;
                useEdgeButton.Layout.Column = 3;
                
                alg = icealgorithm(id);
                im  = alg.loadImage('end');
%                 im = alg.getOriginal;
                alg.Sequencer;
                imagesc(ax,alg.original)
                hold(ax,'on')
                shift = 0;
                if ~isempty(id.Offset)
                    shift = id.Offset;
                end
                plot(ax,id.Edge+shift,1:numel(id.Edge),'k')
                hold(ax,'off')
                axis(ax,'tight','xy')
                colormap(ax,flipud(gray))
                
                function accept(~,~)
                    [~,Z,ds] = id.applyMetric;
                    id.offsetMask.poly     = [x y];
                    id.offsetMask.vector   = interp1(y/ds,x/ds,Z)';
                    delete(f)
                end
                function useEdge(~,~)
                    id.offsetMask.vector = id.Edge+id.Offset;
                    delete(f)
                end
            end

        end
        %% IMAGE FUNCTIONS
        function drawMask(id)
            H = id.parseAnalysis('Hours');
            clf,imagesc(H),colorbar
            caxis([0 inf])
            poly = drawpolygon;
            assignin('base','poly',poly)
        end
        function saveMask(id,polydata)
            H = id.parseAnalysis('Hours');
            id.analysisMask = poly2mask(polydata.Position(:,1),polydata.Position(:,2),size(H,1),size(H,2));
            id.saveData;
        end
        function convert2TIFF(id)
            if isempty(id.Folder.Files)
                return
            end
            files   = id.Folder.Files;
            path    = id.Folder.Path;
            list    = files(contains({files.name},'bmp'));
            if isempty(list)
               return 
            end
            fprintf('')
            exts    = cell(size(list));
            for i=1:numel(list)
                file        = fullfile(path,list(i).name);
                [~,~,ext]   = fileparts(file);
                exts{i}     = ext; 
                im          = imread(file);
                imwrite(im,strrep(file,'bmp','tiff'))
                delete(file)
            end
            exts = unique(exts);
            fprintf('%i/%i files have been converted from %s to TIFF',numel(list),numel(files),strjoin(exts,','))
        end
        function roiRequest(id,mode,imageIndex)
            switch mode
                case 'setROI'
                    image = icealgorithm(id);
                    image.loadImage(imageIndex);
                    if isempty(id.ROI0)
                        id.ROI0 =  [1 1 ...
                            size(image.output,1) ...
                            size(image.output,2) ...
                            id.Controls.values.ROI.rotation];
                    end
                    tmproi = [];
                    if ~isempty(id.ROI)
                        tmproi = id.ROI(1:4);
                        roi1 = [1 1 size(image.output)];
                        roi2 = roi1([1 2 4 3]);
                        if all(tmproi==roi1)||all(tmproi==roi2)
                            tmproi = tmproi/2;
                        end
                    end
                    roi         = getROI(image.output,tmproi,'normal');
                    if all(roi==[1 1 size(image.output)])
                        return
                    end
                    id.ROI = roi;
                    id.ROI(end+1)  = id.Controls.values.ROI.rotation;
                    id.Mask    = [];
                case 'resetROI'
                    id.ROI     = [];
                    image       = icealgorithm(id);
                    image.loadImage(imageIndex);
                    id.ROI0    = [1 1 ...
                        size(image.output,2) ...
                        size(image.output,1) ...
                        id.Controls.values.ROI.rotation];
                    id.ROI = id.ROI0;
                    id.Mask    = [];
                    % currentRot  = id.Controls.values.ROI.rotation;
                    % resetRot    = id.ROI0(5);
                    % newROI0     = ROIrot90(id.ROI0,id.Output,diff([resetRot currentRot]));
                    % id.ROI     = newROI0;
                otherwise
                    image = icealgorithm(id);
                    image.loadImage(imageIndex);
                    if isempty(id.ROI0)
                        id.ROI0 =  [1 1 ...
                            size(image.output,1) ...
                            size(image.output,2) ...
                            id.Controls.values.ROI.rotation];
                    end
                    id.Output = image.output;
            end
        end
        function checkForUpdatedImages(id)
            %% Checks for path designation
            path     = id.Folder.Path;
            hasDrive = numel(split(path,':'))>1;
            if ~hasDrive
               id.Folder.Path = fullfile(id.root,path);
               id.Path = id.Folder.Path;
            end
            %% Updates folder file list
            files = dir(id.Folder.Path);
            files = files(3:end);
            
            if numel(files)==numel(id.Folder.Files)
                return
            end
            
            newfiles  = dir(id.Folder.Path);
            newfiles  = newfiles(3:end);
            [~,order] = sort(datetime({newfiles.date}));
            
            id.Folder.Files = newfiles(order);
            id.times   = datetime({id.Folder.Files.date});
            getDimensions(id) % refresh time resolution
        end
    end
end

