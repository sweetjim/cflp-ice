classdef iwdata<handle
    %IWDATA is a dynamic class for designing and recording information
    %pertaining to internal wave experiments.
    %% Global properties
    properties
        date        datetime
        label       string = ''
        N           (1,1) double
        ds          (1,1) double = 4.27e-4 % m/px
        omega       double
        m0          double = 1/(8.5e-2)
        ss_fcd      struct
        lat         struct
        files       
        times       datetime
        dimensions  struct = struct('x',[],'z',[],'x0',[],'z0',[])
        programs    struct = struct('y',[],'t',[],'details',[])
        experiment  labExperiment
    end
    properties (Hidden)
        ss_dic          struct
        ss_of           struct
        im0             (:,:) double
        im0index        (1,1) double = 1
        refimages       (:,1) double % index array of reference images (blank backgrounds)
        previousindex   logical = false
        automask        logical = false
        last_controls   struct
        last_command    cell
        hm              struct = struct('u',[],'v',[],'t',[]) % Hovmoller
        fcdComputations double = 1
        savepath        char {mustBeFolder(savepath)}=pwd
        maxIndex        (1,1) double
        imagepath       char
        mergeSettings   struct
        ncfilepath      char
        camValue        double
    end
    properties (Hidden,Access = private)
        fnameformat     string = 'cam<device_id|int>_<label|str>C<image_count|int>T<elapsed_time|hh:mm:ss>'
        fileLocation    char {mustBeFolder(fileLocation)}=fileparts(which('iwdata'))%fileparts(matlab.desktop.editor.getActiveFilename)
        hasAutomatedMask logical = false
        hasSavedFile    logical = false
        hasNCFile       logical = false
        maskfill        double
        oldpaths 
        oldlabel
        ed              % experiment designer
    end
    properties (Access = private, Hidden)
        VideoData       struct
        history         struct
        internalCommand logical=false
        %#ok<*TRYNC> 
    end
    %% Constructor / configuration
    methods
        function iw = iwdata(label,loadiwdata)
            iw.initialize;
            if nargin<1
                label    = 'untitled';
            end
            if nargin<2
                loadiwdata = false;
            end
            if loadiwdata
                if isempty(label)
                    load(uigetfile(iw.savepath)); %#ok<LOAD>
                end
                return
            end
            iw.changeDate
            iw.changeLabel(label)
%             if ~isfolder(iw.savepath)
%                 mkdir(iw.savepath)
%             end
        end
        function getN(iw,icedata_iw) % read stratification
            data    = icedata_iw;
            retrieveDetails(data)
            z       = data.Analysis.Z/data.Dimensions.ds;
            rho     = (imresize(data.Details.profile.rho,[numel(z) 1]));
            Ncalc   =@(rho) mean(sqrt(9.81/density(0,0)*abs(gradient(rho,z))));
            Nbulk   =@(rho) sqrt(9.81/density(0,0)*abs((rho(1)-rho(end))/(z(1)-z(end))));
            iw.N   = Ncalc(rho);
        end
        function getSetup(iw) % interactive designer / tank+wave emulator
            iw.runSetup;
        end
        function [x,z,t]=getDimensions(iw,image)
            if nargin<2
                image = iw.im0;
            end
            [y,x] = size(image);
            dim = [[0;0] iw.ds*[x-1 y-1]'];
            x   = linspace(dim(1),dim(3),size(image,2));
            z   = linspace(dim(1),dim(4),size(image,1));
            dt  = mean(diff(iw.times));
            t   = 0:dt:dt*(size(image,3)-1);
        end
    end
     %% Plotting
    methods
        function [imagedata,im,im0,same_state,renderused] = getImageData(iw,index,method,mode,rolling,userender)
            %% Input arguments
            if nargin<3
                method = 'Image';
            end
            if nargin<4
                mode = 'v';
            end
            if nargin<5
                rolling = 0;
            end
            if nargin<6
                userender = false;
            end

            %% Image loading / computations
            controls        = iw.ss_dic.controls;
            controls.index  = index;
            
            % Check for rendered data
            if ~iw.hasNCFile
                userender = false;
            end

            if any(strcmp(mode,{'Image','Difference'}))
                userender = false;
            end
            if userender&&strcmp(method,'fcd')
                userender = iw.nc_read(index,'rend',true);
            end
            renderused  = false;
            if userender
                iw.nc_read(index,mode);
                renderused = true;
            end
            
            im  = [];
            im0 = [];

            %same_state = false;
            same_state = isequaln(iw.last_controls,controls);
            if ~renderused&&~same_state
                [im,im0]    = iw.loadimage(index);
            end

            switch mode
                case 'Image'
                    imagedata = im;
                case 'Difference'
                    imagedata =  im-im0;
                otherwise
                    switch method
                        case 'fcd'
                            % Calculate field
                           
                            if ~same_state&&~userender
                                
                                if strcmp(mode,'Q')
                                    mode = 'h';
                                end
                                mode = lower(mode);
                                if rolling>0
                                    iw.ss_movmean_fcd(index,rolling);
                                else
                                    iw.ss_solver_fcd(im0,im);
                                end
                            end
                            imagedata = eval(sprintf('iw.ss_fcd.%s;',mode));
                        case 'of'
                            iw.ss_solver_of(im0,im)
                            imagedata = eval(sprintf('iw.ss_of.%s;',mode));
                    end
            end
% 
%             try
%                 imagedata(~iw.ss_fcd.mask)=nan;
%                 im0(~iw.ss_fcd.mask)=nan;
%             end
        end
        function imagedata = imageParser(iw,index,varargin)
            %% Input
            ax          = [];
            method      = 'fcd';
            mode        = 'Image';
            overlay     = 'None';
            userender   = false;
            renderdata  = false;
            rolling     = 0;
            flag        = false;
            parseInput(varargin)

            if any(strcmp(mode,{'Image','Difference'}))
                renderdata = false;
            end

            [imagedata,im,im0,same_state,renderused] = getImageData(iw,index,method,mode,rolling,userender);
            %% NetCDF Rendering
            if renderdata&&~iw.nc_check;iw.nc_create;end

            %% Date sorting
            if isempty(iw.times)
                iw.updateTimes
            end
            
            %% Outputs & rendering
            if renderdata&&~renderused
%                 iw.nc_render(index,im,'im')
                if strcmp(method,'fcd')
                    iw.nc_render(index,iw.ss_fcd.u,'u')
                    iw.nc_render(index,iw.ss_fcd.v,'v')
                end
            end
            if nargout>1
                return
            end
            %% Axis handling
            if isempty(ax)
                ax  = gca;
            end

            imdata = findall(ax,'Type','Image');
            if isempty(imdata)
                imdata=imagesc(ax,im);
            end
            imdata.AlphaData = true;
            controls.index  = index;
            
            %% Display 
            cbar = findall(ancestor(ax,'Figure'),'Type','Colorbar');
            hasCbar = ~isempty(cbar);
            switch mode
                case 'Image'
                    imdata.CData = im;
                    caxis(ax,[0 Inf])
                    colormap(ax,'gray')
                case 'Difference'
                    imdata.CData = im-im0;
                    caxis(ax,[-20 20])
                    if ~hasCbar
                        addColorbar('ax',ax,'cmap','balance','pivot',0,'invert')
                    else
                       colormap(ax,cmapCheck('balance'))
                       caxis(ax,[-1 1]*std(im-im0,[],'all','omitnan'))
                    end
                otherwise
                    if ~same_state && ~userender
                        switch method
                            case 'dic'
                                iw.ss_solver_dic(im0,im,...
                                    'blockspacing',controls.bs,...
                                    'blockoverlap',controls.bo,...
                                    'maxdisp',controls.md,...
                                    'minquality',controls.mq/1e3,...
                                    'normcorr',controls.nc,...
                                    'fast',controls.fast);
                                iw.ss_dic.im0 = im0;
                                iw.ss_dic.im  = im;
                            case 'fcd'
                                if iw.ss_fcd.fail
                                    imdata.CData = double(im)-double(im0);
                                    caxis(ax,[-20 20])
                                    if ~hasCbar
                                        addColorbar('ax',ax,'cmap','balance','pivot',0,'invert')
                                    end
                                    return
                                end
                        end
                        flag = true;
                    end
            end

            switch method
                case 'fcd'
                    ss = iw.ss_fcd;
                case 'dic'
                    ss = iw.ss_dic;
                case 'of'
                    ss = iw.ss_of;
                case 'lat'
                    ss = iw.lat;
                    if all(~strcmp(mode,{'Image','Difference'}))
                        mode = 'Difference';
                    end
                otherwise
                    ss = [];
            end
            plotting(ss,mode)
            
            % dimensions
            dim =  size(imdata.CData)*1/(iw.ss_dic.controls.resize*1e-2)*iw.ds;
            imdata.XData = [0 dim(2)];
            imdata.YData = [0 dim(1)];
            % alpha mask
            if ~any(strcmp(mode,{'Image','Difference'}))
                if ~userender
                    try
                        imdata.AlphaData = imresize(iw.ss_fcd.mask,iw.ss_dic.controls.resize*1e-2);
                    end
                end
            end
            %% Overlay
            hold(ax,'on')
            if ~strcmp(overlay,'None')
                U = smooth2a(gather(ss.u),ceil(size(ss.u,1)/200),ceil(size(ss.u,2)/200));
                V = smooth2a(gather(ss.v),ceil(size(ss.u,1)/200),ceil(size(ss.u,2)/200));
                zerocutoff = -1.3;
                U = U.*double(log10(abs(U))>zerocutoff);
                V = V.*double(log10(abs(V))>zerocutoff);
            end

            switch overlay
                case {'Velocity vectors','vectors'}
                    vis_flow(U,V,30, 1, 2, 'm',1,ax);
                case 'Streamlines'
                    streamslice(ax,U,V,2)
            end
            hold(ax,'off')
            %% Labels
            try
                title(ax,string(diff([iw.times(1) iw.times(index)])))
                [~,fname] = fileparts(iw.files.Files(index));
                subtitle(ax,fname,'Interpreter','none')
            end
            %% Flagging
            if flag
                iw.last_controls        = controls;
                iw.last_controls.index  = index;
            end
            %% Saving
            assignin('base','iw',iw)
            assignin('base','ss',ss)
            assignin('base','im',im)
            assignin('base','im0',im0)
            assignin('base','ax',ax)

            if ~nargout
                clear imagedata
            end
            %% Functions
            function plotting(in,mode) 
                scale   = 'linear';
                cmap    = 'balance';
                if isempty(in)
                    return
                end
                switch mode
                    case {'Image','Difference'}
                        return
                    case {'u','v'}
                        eval(sprintf('out=in.%s;',mode));
                        if strcmp(method,'dic')
                            out = smooth2a(out,5,5); %#ok<NODEF>
                        end
                        cmode = 'pos-neg';
                        if round(sum(colormap(ax),[1 2]))==384
                            colormap(ax,cmapCheck('balance'))
                        end
                    case {'H','h'}
                        eval(sprintf('out=in.%s;',lower(mode)));
                        cmode = 'pose-neg';
                    case {'Q','q','U'}
                        eval(sprintf('out=in.%s;',lower(mode)));
                        cmap    = 'gray';
                        cmode    = 'pos';
                    case {'vorticity','zeta'}
                        U = smooth2a(gather(ss.u),ceil(size(ss.u,1)/200),ceil(size(ss.u,2)/200));
                        V = smooth2a(gather(ss.v),ceil(size(ss.u,1)/200),ceil(size(ss.u,2)/200));
                        zerocutoff = -1.3;
                        U = U.*double(log10(abs(U))>zerocutoff);
                        V = V.*double(log10(abs(V))>zerocutoff);
                        out     = vorticity(U,V);
                        cmode = 'pos-neg';
                end
                alphamap = ~isnan(out);%.*in.q>0.5;
                imdata.CData        = out;
                imdata.AlphaData    = alphamap;
                set(ax,'ColorScale',scale)

                sigma2 = std(out,[],'all','omitnan')*2;
                if isempty(ax)
                    ax  = gca;
                end

                cmaxref = [-1 1]*.1;
                if ~isempty(ax.UserData)
                    cmaxref = ax.UserData;
                end
                switch cmode
                    case 'pos-neg'
                        cmax = sigma2*2.*[-1 1];
                    case 'pos'
                        cmax = sigma2*2.*[0 1];
                end

                if contains(mode,{'Q','q'})
                    caxis(ax,[-inf inf])
                end

                cbar = findall(ancestor(ax,'Figure'),'Type','Colorbar');
                set(ax,'Colormap',cmapCheck(cmap))
                if isempty(cbar)
                    try
                        addColorbar('ax',ax,'cmap',cmap,'pivot',0,'invert')
                    catch
                        caxis(ax,std(out,[],'all','omitnan')*2.*[0 1])
                        addColorbar('ax',ax,'cmap',cmap)
                    end
                end

                caxis(ax,cmax)
                return
                switch cmode
                    case 'pos'
                        caxis(ax,cmax)
                    otherwise
                        caxis(ax,[-1 1]*.5)
                end
                if all(abs(cmax)>abs(cmaxref))
                    %ax.UserData = cmax;
                    caxis(ax,[-1 1]*.5)
%                     caxis(ax,cmax)
                end
                if isempty(ax.UserData)
                    ax.UserData = cmax;
                end
            end
            %% Input parser
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case 'ax'
                            ax = namevalue;
                        case {'Image','u','v','U','vorticity','Q','zeta','q','Difference','H','h'}
                            mode = items{m};
                        case 'overlay'
                            overlay = namevalue;
                        case 'method'
                            method = namevalue;
                        case 'rolling'
                            rolling = namevalue;
                        case {'userender','usenc','nc'}
                            userender = namevalue;
                        case 'render'
                            renderdata = namevalue;
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
    %% Video writing
    methods (Access = public)
        function WriteVideo(iw,executeVideo,filename,index,varargin)
            fname = join([string(datetime,'yyyy_MM_dd_hh_mm') '_video.avi']);
            if nargin<2
                executeVideo = false;
            end
            if nargin<3
                filename = fname;
            end
            if nargin<4
                index = 1:iw.maxIndex;
            end
            if isempty(filename)
                filename = fullfile(pwd,fname);
            end
            if isfolder(filename)
                filename = fullfile(filename,fname);
            end
            if iscell(varargin)&&numel(varargin)==1
                varargin = varargin{:};
            end

            index(index<1)=1;
            index(index>iw.maxIndex)=iw.maxIndex;

            fig = findall(0,'Type','Figure','Tag','VideoWriter');
            if isempty(fig)
                fig = figure('Name','Video Writer', ...
                    'Tag','VideoWriter', ...
                    'HandleVisibility','callback', ...
                    'Color','w', ...
                    'IntegerHandle','off', ...
                    'ToolBar','none');
            end
            clf(fig)
            fig.Name = sprintf('Video Writer: %s',filename);
            ax = axes(fig);
            iw.imageParser(index(1),'ax',ax,varargin{:});
            axis(ax,'image')
            imageformatting(ax)
            cbar = colorbar(ax);

            if ~executeVideo
                return
            end
            delete(cbar)
            if numel(index)<2
                error('Cannot write a video of unitary or invalid indicies')
            end
            
            V           = VideoWriter(filename);
            V.FrameRate = 60;
            V.Quality   = 75;
            open(V)
            for i=index
                warning off backtrace
                iw.imageParser(i,'ax',ax,varargin{:});
                imageformatting(ax)
                writeVideo(V,getframe(ax))
                displayProgress('Writing video',i,index(1),index(end))
            end
            warning on
            close(V)
            %winopen(fileparts(filename))
            function imageformatting(ax)
                ax.XTick = [];
                ax.YTick = [];
                clim(ax,[-1 1]*.5)
                colormap(ax,cmapCheck('gray'))
            end
        end
    end
    %% Light attenuation technique (lat)
    methods
        function plot_lat(iw,ax)
            LAT = iw.lat;
            z = iw.metric(LAT.im0);
            plot(ax,z,mean(LAT.im0)-mean(LAT.im),'r',z,z*0,'--k')
            ylim(ax,[-1 1]*5)
            
        end
    end
    %% Synthetic Schlieren (ss)
    methods % Synthetic Schlieren functions
        function b              = ss_buoyancy_field(iw,u,v)
            %% Calculates the buoyancy field from displacement vector field
            W       = 0.2;          % [m] width of control volume
            B       = 0.22;         % [m] distance of reference map from control volume
            gamma   = W^2/2+W*B;    % [m2] adjustment for physical dimensions
            rho0    = 1e3;          % [kg/m3] reference density
            ds      = iw.ds;        % [px2m] image spatial resolution
            beta    = 0.184;

            if nargin==1
                u = iw.ss_fcd.u;
                v = iw.ss_fcd.v;
            end

            % calculate density anomaly gradient fields
            b = fftinvgrad(-u*ds,-v*ds)*rho0/(beta*gamma);
            % da_x    = cumsum(u*ds*rho0/(beta*gamma))*ds;
            % da_y    = cumsum(v*ds*rho0/(beta*gamma))*ds;
            % 
            % b = 9.81/1e3*(rho_ref+cat(3,da_x,da_y));
        end
        function [u,v]          = ss_movmean_fcd(iw,index,window)
            %% Rolling average
            if window==0
                return
            end

            % check for rolling step
            controls = iw.ss_dic.controls;
            try 
                step = controls.fcdrollingstep;
            catch
                step = 1;
            end

            % Get valid indicies
            indicies                        = index-window:step:index+window;
            indicies0                       = indicies;
            indicies(indicies<1)            = 1;
            indicies(indicies>iw.maxIndex)  = iw.maxIndex;
            
            % check memory for previous data
            [keep,omit] = checkHistory(iw,indicies);
            if any(keep)
                indicies = indicies(omit);
            end
            indicies = unique(indicies);

            % Initialise arrays
            [im,IM0]    = iw.loadimage(indicies);
            U = zeros(size(im));
            V = U;
            H = U;

            % Loop over indicies
            for i=1:numel(indicies)
                [u,v] = ss_solver_fcd(iw,IM0,im(:,:,i));
                U(:,:,i)=u;
                V(:,:,i)=v;
                H(:,:,i)=iw.ss_fcd.h;
            end

            % efficiency logging
            iw.fcdComputations = numel(indicies)*iw.fcdComputations;


            % Update array with other renders in memory
            if any(keep)
                iw.history.U(:,:,omit) = U;
                iw.history.V(:,:,omit) = V;
                iw.history.H(:,:,omit) = H;

                % indicies in memory
                mem = iw.history.indicies(~omit);
                % indicies to add
                add = indicies0(~ismember(indicies0,mem));
                

                indicies = sort([mem add]);
                U = iw.history.U;
                V = iw.history.V;
                H = iw.history.H;
            end

            % Subtract mean field
            u = mean(U,3)-U(:,:,indicies==index);
            v = mean(V,3)-V(:,:,indicies==index);
            h = mean(H,3)-H(:,:,indicies==index);

            % Output to fcd structure
            iw.ss_fcd.u = u;
            iw.ss_fcd.v = v;
            iw.ss_fcd.h = h;
            iw.ss_fcd.U = hypot(u,v);

            % Update history
            iw.history = struct( ...
                'indicies',indicies,...
                'index',index, ...
                'U',U, ...
                'V',V, ...
                'H',H,...
                'resize',controls.resize);

            %% Subroutines
            function [keep,omit] = checkHistory(iw,indicies)
                keep = false;
                omit = false;
                if isempty(fieldnames(iw.history))
                    return
                end
                if iw.ss_dic.controls.resize~=iw.history.resize
                    return
                end
                indh = iw.history.indicies;
                if numel(indh)~=numel(indicies)
                    return
                end
                try
                    keep    = within(indh,index+[-1 1]*window);
                    base    = iw.history.indicies==iw.history.index;
                    omit    = ~xor(keep,base);
                end
            end
        end
        function [u,v,h]          = ss_solver_fcd(iw,Iref,Idef,notRecursion)
            %% Controls
            c = iw.ss_dic.controls;
            useAutomatedMask = all([iw.hasAutomatedMask,iw.automask]);
            usehamming  = false;
            thresh      = 0.5;
            kmult       = 1;
            try
                usehamming = c.hamming;
            end
            try
                thresh = c.thresh;
            end
            try
                kmult= c.kmult;
            end
                        
            %% Attempt mask subroutine
            if nargin<4
                notRecursion = true;
                % Efficiency logging
                iw.fcdComputations = 0;
            end
            fcd = iw.ss_fcd;


%             m = fcd.mask;
%             p = iw.ss_mask_fill.*~m;
%             Iref(~m)=0;
%             Idef(~m)=0;
%             Iref = Iref+p;
%             Idef = Idef+p;

            if useAutomatedMask
                MASK = imresize(iw.ss_fcd.mask,size(Iref));
                softmask = imgaussfilt(double(MASK),10);
                Idef = Idef.*softmask;
                Iref = Iref.*softmask;
            end

            if ~isempty(fcd.mask)&&notRecursion&&~iw.hasAutomatedMask
                [u,v] = maskRoutine(iw,fcd,Iref,Idef);
                return
            else
                iw.fcdComputations = iw.fcdComputations+1;
            end
            %% Main routine
            
            % get two independent carrier peaks from reference image
            iw.ss_fcd.fail = false;
            kr = [];
            ku = [];
            try
                kr = iw.ss_fcd.K(1,:);
                ku = iw.ss_fcd.K(2,:);
            end
            try
                [kr, ku] = findorthcarrierpks(Iref, 4*pi/min(size(Iref)), Inf,thresh,usehamming);
                iw.ss_fcd.K = [kr;ku];
            catch
                if isempty(ku)
                    %                 shortWarning('No carrier signal detected.')
                    iw.ss_fcd.fail = true;
                    return
                end
            end
            
            % extract carrier signals from reference image and store them for later use
            krad    = sqrt(sum((kr-ku).^2))/2*kmult;%2;
            fIref   = fft2(Iref);
            cr      = getcarrier(fIref, kr, krad);
            cu      = getcarrier(fIref, ku, krad);

            % gpu acceleration
%             cr = gpuArray(cr);
%             cu = gpuArray(cu);

            
            fIdef   = fft2(Idef);
            [u,v]   = fcd_dispfield(fIdef,cr,cu);
            
            % Remask
            if useAutomatedMask
                Iref = Iref.*MASK;
                Idef = Idef.*MASK;
                u = u.*MASK;
                v = v.*MASK;
            end

            % get displacement field and height profile
            h       = fftinvgrad(-u,-v);%,'bcFix','mirror');            
            %% Structure memory
            iw.ss_fcd.im0  = Iref;
            iw.ss_fcd.im   = Idef;
            iw.ss_fcd.cr   = cr;
            iw.ss_fcd.cu   = cu;
            iw.ss_fcd.u    = u;
            iw.ss_fcd.v    = v;
            iw.ss_fcd.U    = sqrt(u.^2+v.^2);
            iw.ss_fcd.h    = h;

            %% Subroutines
            function [u,v]=maskRoutine(iw,fcd,Iref,Idef)
                % Initialize arrays
                u0      = zeros(size(Iref));
                v0      = u0;
                h0      = u0;

                % Get mask and distinct regions of mask
                mask    = imresize(fcd.mask,iw.ss_dic.controls.resize*1e-2);
                props   = regionprops(mask);
                
                % Loop over each distinct region
                for i=1:numel(props)
                    % Get bounding boxes and ensure limits are valid
                    roi         = round(props(i).BoundingBox);
                    roi(roi<1)  = 1;
                    if roi(3)>size(Iref,2);roi(3)=size(Iref,2);end
                    if roi(4)>size(Iref,1);roi(4)=size(Iref,1);end

                    % Apply bounding box to image
                    iref = imcrop(Iref,roi);
                    idef = imcrop(Idef,roi);

                    % Solve FCD on region
                    iw.ss_solver_fcd(iref,idef,false);
                    xrange = roi(1):roi(3)+roi(1);
                    zrange = roi(2):roi(4)+roi(2);
                    if iw.ss_fcd.fail
                        continue
                    end

                    % Integrate results to main image
                    u0(zrange,xrange) = iw.ss_fcd.u;
                    v0(zrange,xrange) = iw.ss_fcd.v;
                    h0(zrange,xrange) = iw.ss_fcd.h;
                end

                % Save results to FCD reference structure
                iw.ss_fcd.im0 = Iref;
                iw.ss_fcd.im  = Idef;
                iw.ss_fcd.u = u0;
                iw.ss_fcd.v = v0;
                iw.ss_fcd.U = hypot(u0,v0);
                iw.ss_fcd.h = h0;

                % Function outputs
                u = u0;
                v = v0;
            end
        end
        function [u,v]          = ss_solver_of(iw,Iref,Idef)
            OPS = iw.ss_of;
            output = velocity_field(...
                Iref,Idef,...
                OPS.method,...
                OPS.alpha,...
                OPS.iterations,...
                OPS.lambda,...
                OPS.illumination,...
                OPS.filter,...
                OPS.diffmethod);
            iw.ss_of = catstruct(iw.ss_of,output);
        end
        function status         = ss_plot_fcd(iw,ax)
            %%
            fcd = iw.ss_fcd;
            [rows,cols] = size(iw.im0);
            kxvec = fftshift(kvec(cols));
            kyvec = fftshift(kvec(rows));
            wr = hann(rows,'periodic');
            wc = hann(cols,'periodic');
            win2d = wr(:)*wc(:)';
            im0   = double(iw.im0);
            fftIm = fftshift(abs(fft2((im0-mean(im0(:))).*win2d)));
            imdata=findall(ax,'Type','Image');
            imagesc(ax,kxvec, kyvec, fftIm,[0,max(fftIm(:))/50/4]);
            if isempty(imdata)
                set(ax, ...
                    'XTickLabelMode','auto','XTick',unique(round(kxvec)), ...
                    'YTickLabelMode','auto','YTick',unique(round(kyvec)),...
                    'ColorScale','log','Colormap',flipud(gray))
                addlabels(ax,'x','$k_x$','y','$k_y$','latex')
            end
            hold(ax,'on')
            try
                fcd.cr.plot('Parent',ax,'color','r','marker','none')
                fcd.cu.plot('Parent',ax,'color','r','marker','none')
            catch
                %                 warning('An error occured')
            end
            hold(ax,'off')
            axis(ax,'tight')
            status = [];
        end
        function [u,v,Xq,Yq,q]  = ss_solver_dic(iw,Iref,Idef,varargin)
            useparfor = false;
            fast = 1;   % use cubic+symmetric (slow) or linear+none (fast) in solver
            bs = 4;     % resolution
            bo = 4;     % blocksize = blockspacing + 2*borderoverlap
            md = 4;     % sets size of search window = blocksize + 2*maxdisp
            mq = 0.01;  % minimum quality of the correlation to be a valid vector, NaN otherwise
            nc = 1;     % use normalized cross correlation? relatively slow but robust
            parseInput(varargin)
            disp('Processing')
            if mq==0
                mq=nan;
            end

            switcher = contains(class(Iref),'gpu')||fast;
            warpmethod = 'linear';
            if ~switcher
                warpmethod = 'cubic';
            end
            %% Solver
            tic
            useparfor = 0;
            if useparfor
                [blocksz,newsize] = parallelize(Iref,Idef);
                Iref    = imresize(Iref,newsize);
                Idef    = imresize(Idef,newsize);
                divisions = newsize./blocksz;
                overlap = 6+floor(bs/2)+md.*[1 1];

                Irefbk  = blockify(Iref);
                Idefbk  = blockify(Idef);

                U       = zeros([blocksz prod(divisions)]);
                V = U; %XQ = U; YQ = U; Q = U;

                % permute for parfor indexing
                Irefbk = permute(Irefbk,[3 1 2]);
                Idefbk = permute(Idefbk,[3 1 2]);

                tic
                parfor ip=1:prod(divisions)
                    [u,v,Xq,Yq,q]=ss_solver_dic_block(Irefbk(ip,:,:),Idefbk(ip,:,:),bs,bo,md,mq,nc,warpmethod,switcher);
                    U(ip,:,:) = u;
                    V(ip,:,:) = v;
                end
                toc
                []
                return
                fhandle = @(block_struct) blockfcn(block_struct,Idef,bs,bo,md,mq,nc,warpmethod,switcher);
                tic
                out     = blockproc(Iref,blocksz,fhandle,...
                    'BorderSize',6+floor(bs/2)+md.*[1 1],...
                    'UseParallel',1,...
                    'DisplayWaitbar',0);
                clc
                toc
                imagesc(out)
                [];
            else
                [u,v,Xq,Yq,q]=solver_routine(Iref,Idef,bs,bo,mq,nc);
            end
            toc

            %% iwect memory
            iw.ss_dic.q = q;
            iw.ss_dic.u = u;
            iw.ss_dic.v = v;
            iw.ss_dic.U = sqrt(u.^2+v.^2);
            iw.ss_dic.Xq = Xq;
            iw.ss_dic.Zq = Yq;
            %% Functions
            function [u,v,Xq,Yq,q]=solver_routine(Iref,Idef,bs,bo,mq,nc)
                [u, v, br, bc, q] = dic_dispfield(Iref, Idef, bs, bo, md, [], mq, nc);


                filtering(switcher)

                % warp cycle(s) if required
                Iref_w = interpimwarp(Iref, u, v, bc, br, warpmethod);
                smalldisp = 2; % allow small displacement around current solution
                [du, dv, ~, ~, q] = dic_dispfield(Iref_w, Idef, bs, bo, md, smalldisp, mq, nc);
                u = u + du;
                v = v + dv;

                filtering(switcher)

                % final interpolation for display
                [Xq, Yq] = meshgrid(1:size(Iref,2), 1:size(Iref,1));
                % u = interp2(bc,br,u,Xq,Yq,'linear',0);
                % v = interp2(bc,br,v,Xq,Yq,'linear',0);
                [Iref_w, u, v] = interpimwarp(Iref, u, v, bc, br, warpmethod);

                if switcher
                    try
                        [du, dv] = of_dispfield(gpuArray(Iref_w),gpuArray(Idef), .2);
                        du = gather(du);
                        dv = gather(dv);
                    catch
                        [du, dv] = of_dispfield(Iref_w,Idef, .2);
                    end
                else
                    roi = ~isnan(u);
                    [du, dv] = of_dispfield(Iref_w,Idef, .1, roi);
                end
                u = u + du;
                v = v + dv;
                clc
                function filtering(switcher)
                    if switcher
                        try
                            u = medfilt2(u);
                            v = medfilt2(v);
                            return
                            %                         catch
                            %                             flag = true;
                        end
                    end
                    u = medfilt2(gather(u),'symmetric');
                    v = medfilt2(gather(v),'symmetric');
                end
            end

            function [blocksize,newsize] = parallelize(Iref,Idef,bs,mq)
                % Initiate parallel pool
                if isempty(gcp)
                    parpool('IdleTimeout',600)
                end
                tmp         = gcp;
                divs        = tmp.NumWorkers;
                clear tmp
                blocks      = factor(divs);
                c           = arrayfun(@(x)length(find(blocks == x)), unique(blocks), 'Uniform', false);
                dim         = (cell2mat(c).*unique(blocks));
                blocksize   = ceil(size(Iref)./dim);
                newsize     = blocksize.*dim;

            end
            function out = blockify(Iref)
                IrefB   = zeros([blocksz divisions]);
                % organise blocks in x-y-i-j array
                for i=1:divisions(1)
                    for j=1:divisions(2)
                        range_x = blocksz(1).*[i-1 i];x = range_x(1)+1:range_x(2);
                        range_y = blocksz(2).*[j-1 j];y = range_y(1)+1:range_y(2);
                        IrefB(:,:,i,j) = Iref(x,y);
                    end
                end

                % shift j into i array for 1D-parfor loop
                % s.t.  x-y-(i1,i2,..in)-(j1,j2,...jn) becomes
                %       x-y-(i1,j1,j2...jn,i2,j1,j2,...jn,...)
                %  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
                % |i1,j1 | i1,j2 |  ...  | i1,jn |
                % |i2,j1 | i2,j2 |  ...  | i2,jn |
                % | ...  |  ...  |  ...  | i1,jn |
                %  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
                %
                % becomes
                %  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
                % |i1,j1 | i1,j2 |  ...  | i1,jn | i2,j1 | i2,j2 |  ...

                IrefB2 = zeros([blocksz prod(divisions)]);
                IrefB2(:,:,1:divisions(2)) = squeeze(IrefB(:,:,1,:));
                for i=2:divisions(1)
                    to_cat      = squeeze(IrefB(:,:,i,:));
                    sequence    = divisions(2)*(i-1)+1:divisions(2)*(i);
                    IrefB2(:,:,sequence) = to_cat;
                end

                out = IrefB2;
            end
            function out = blockfcn(BS,idef,bs,bo,md,mq,nc,warpmethod,switcher)
                im0blk = double(BS.data);
                startInd = BS.location;
                endInd   = BS.location+BS.blockSize-1;
                imblk = idef(startInd(1):endInd(1), startInd(2):endInd(2));
                [uu,vv,Xxq,Yyq,qq]=ss_solver_dic_block(im0blk,imblk,bs,bo,md,mq,nc,warpmethod,switcher);
                out      = uu;
            end
            %% Input parser
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case 'blockspacing'
                            bs = namevalue;
                        case 'borderoverlap'
                            bo = namevalue;
                        case 'maxdisp'
                            md = namevalue;
                        case 'minquality'
                            mq = namevalue;
                        case 'normcorr'
                            nc = namevalue;
                        case 'parallel'
                            useparfor = true;
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
        function roi            = setROI(iw,type)
            if nargin<2
                type = 'set';
            end

            switch type
                case 'set'
                    if isempty(iw.im0)
                        iw.im0 = handleMerged(iw,1);
                    end
                    try
                        roi = getROI(iw.im0,iw.ss_dic.controls.roi);
                    catch
                        roi = getROI(iw.im0);
                    end
                otherwise
                    roi = [];
            end

            if ~isempty(iw.ss_fcd.mask)
                iw.hasAutomatedMask = false;
                iw.ss_fcd.mask      = false;
                shortWarning('FCD mask has been deleted')
            end

            iw.ss_dic.controls.roi = roi;
        end
        function ss_mask_fcd(iw,automatic,setAutoparameters)
            %% Opens a program to define masking areas 
            % or automatically determines a masking space
            if nargin<2
                automatic = false;
            end
            if nargin<3
                setAutoparameters = false;
            end

            if automatic
                index   = iw.im0index;
                samples = 30;
                if iw.hasAutomatedMask
                    q = questdlg('Overwrite automated mask?', ...
                        'Automated mask decision', ...
                        'Yes','No','No');
                    if strcmp(q,'No')
                        return
                    end
                end
                if setAutoparameters
                    q = inputdlg({'Reference image (index)','Block size (px)'}, ...
                        'Automated mask inputs',[1 45;1 45],{num2str(index),num2str(samples)});
                    if isempty(q)
                        return
                    end
                    try
                        index   = eval(q{1});
                        index(index<0)=0;
                        index(index>iw.maxIndex)=iw.maxIndex;
                    catch
                        error('Invalid reference index input')
                    end
                    try
                        samples = eval(q{2});
                    catch
                        error('Invalid block size input')
                    end
                end

                ss_mask_fcd_automated(iw,index,samples)
                return
            end

            fig = findall(0,'Type','Figure','Name','FCD Mask');
            fcdapp = app_fcdMask(iw);

            if iw.hasAutomatedMask
                iw.ss_fcd.mask = false;
                iw.hasAutomatedMask = false;
                iw.automask = false;
            end

            try
                % load user-defined masks
                if isfield(iw.ss_fcd,'maskrois')
                    % from images.roi data
                    fcdapp.rois = iw.ss_fcd.maskrois;
                elseif isfield(iw.ss_fcd,'mask')
                    % from logical mask data
                    mask  = iw.ss_fcd.mask;
                    if isempty(mask)
                        return
                    end
                    masks = regionprops(mask);
                    masks(1).type      = 'images.roi.rectangle';
                    masks(1).inverted  = 'Add';
                    for i=1:numel(masks)
                        masks(i).position = masks(i).BoundingBox;
                        masks(i).type      = 'images.roi.rectangle';
                        masks(i).inverted  = 'Add';
                    end
                    fcdapp.rois = masks;
                end
            end
        end
        function pattern = ss_mask_fill(iw)
            %% Fills masked areas with a repeated background pattern from the image            % get background pattern
            
            % begin search
            window = round(size(iw.im0)/10);
            val = arrayfun(@(idx) fcn(iw.im0,window,idx), ...
                min(window):min(window):min(window)*10,'UniformOutput',false)';
            val = val(~cellfun(@isempty,val));
            pattern = iw.im0(array(val{2}(1,:)),array(val{2}(2,:)));
            return
            pattern = repmat(pattern,[10 10]);

            [x,y] = size(iw.im0);
            pattern = pattern(array([1 x]),array([1 y]));

            function val = fcn(Iref,window,idx)
                xrange = idx:idx+window(1);
                yrange = idx:idx+window(2);
                try
                    Iref    = Iref(xrange,yrange);
                catch
                    val = [];
                    return
                end
                try
                    findorthcarrierpks(Iref, 4*pi/min(size(Iref)), Inf);
                catch
                    val = [];
                    return
                end
                val = [min(xrange) max(xrange)
                    min(yrange) max(yrange)];
            end

        end
        function hovmoller(iw)
            t       = minutes(iw.times-iw.times(1));
            loop    = roundto(t,1);
            %             loop    = loop(1:50);
            %             loop    = 1:numel(iw.files.Files);
            im0     = loadimage(iw,iw.im0index);

            u       = zeros([numel(loop) size(im0)]);
            v       = u;
            j       = 1;
            Ucalc   = @(u,v) sqrt(u.^2+v.^2);
            %%
            for i=1:numel(loop)
                [ui,vi]=ss_solver_fcd(iw,im0,iw.loadimage(loop(i)));
                %                 imageParser(iw,i,'u')
                %                 tmp = iw.ss_fcd;
                u(i,:,:) = ui;
                v(i,:,:) = vi;
                displayProgress('Progress',i,1,numel(loop))
                %                 j       = j+1;
            end
            u = permute(u,[2 3 1]);
            v = permute(v,[2 3 1]);
            %%
            iw.hm.u = u;
            iw.hm.v = v;
            iw.hm.U = Ucalc(u,v);
            iw.hm.t = loop;
            %%
            imagesc(t(loop),1:size(im0,2),squeeze(mean(iw.hm.U,1)))
            addlabels('x','Time (minutes)','y','Depth (px)')
        end
    end
    methods (Access = private)
        function ss_mask_fcd_automated(iw,index,samples)
            %SS_MASK_FCD_AUTOMATED attempts to find regions on the image
            %which have a common synthetic schlieren background by a
            %sub-sampling algorithm

            ref         = iw.loadimage(index,true);
            bim         = blockedImage(ref,'BlockSize',samples.*[1 1]);
            im          = apply(bim,@beginSearch);
            im          = gather(im);
            im          = rescaleMax(im);

            % look for most common carrier
            [c,v]       = imhist(im,500);
            % omit false readings
            c(1)        = []; 
            v(1)        = [];

            [~,idx]=max(c);
            w = floor(pulsewidth(c));
            valid = within(im,v(idx-w:idx+w));%(idx-1:idx+1));

            iw.ss_fcd.mask = valid;
            
            iw.hasAutomatedMask = true;
            function d = beginSearch(bim)
                d = bim.Data;
                try
                    
                    val = findorthcarrierpks(d,4*pi/min(size(d)), Inf,.5,true);
                    d(:,:)= hypot(val(1),val(2));
                catch
                    d(:,:)=0;
                end
            end
        end
    end
    %% Initialization / design
    %% Image processing
    methods % Image loading (imagedatastore)
        function [im,im0] = loadimage(iw,index,noresize,displayProgess)
            if nargin<3
                noresize = false;
            end
            if nargin<4
                displayProgess = false;
            end
            
            [im0,im] = handleMerged(iw,index);
            controls = iw.ss_dic.controls;

            im0 = transformImage(im0);
            im = arrayfun(@(x) transformImage(im(:,:,x)), ...
                1:size(im,3),'UniformOutput',false);
            im = cat(3,im{:});

            function im=transformImage(im)
                if ~isempty(controls.roi)
                    im   = imcrop(im,controls.roi);
                end

                if controls.rotation~=0
                    im  = rot90(im,controls.rotation);
                end

                if controls.resize~=100&&~noresize
                    im   = imresize(im,controls.resize/1e2);
                end
            end
            
        end
        function [im0,im] = handleMerged(iw,index)
            % Switch between "imMerge" class (two indepedendent datastores)
            % and "ImageDatastore" class

            switch class(index)
                case 'double'
                    % limit index
                    index(index<1)=1;
                    index(index>iw.maxIndex)=[];
                    if isempty(index)
                        index = iw.maxIndex;
                    end
                case {'char','string'}
                    index = 1:floor(iw.maxIndex/str2double(index)):iw.maxIndex;
            end

            % switch image reading mode (merged image/multistream image or
            % singular image)
            condition = isa(iw.files,'imMerge');
            refindex = iw.im0index;
            if iw.previousindex
                refindex = index-1;
            end
            refindex(refindex<=0)=1;
            
            if numel(index)==1
                iw.im0  = switchReadMode(iw,refindex);
            end
            im0     = iw.im0;
            im      = switchReadMode(iw,index);

            function im = switchReadMode(iw,index)
                if numel(index)>1
                    im = readDatastoreAtIndices(iw.files,index);
                    return
                end

                switch condition
                    case true
                        im = iw.files.readImageAtIndex(index);
                    case false
                        im = readimage(iw.files,index);
                end

                im = im2double(im);
            end
        end
        function [x,z,ds] = metric(iw,ref)
            ds = iw.ds;
            if isempty(iw.im0)
                X = [];
                Z = [];
                return
            end
            if nargin<2
                ref = iw.im0;
            end
            [Z,X]   = size(ref);
            if isempty(ds)
                ds  = 1;
                x   = 1:X;
                z   = 1:Z;
            else
                ds  = iw.ds;
                x   = linspace(0,X*ds,X);
                z   = linspace(0,Z*ds,Z);
            end
        end
        
    end
    
    %% Wave detection
    methods 
        function [W,Z,t,T,Fk,om,m] = wave_analyseColumn(iw,index,x_pos,isPlotting)
            %%WAVE_ANALYSECOLUMN extracts wave signals from the
            %%displacement field of the synthetic schlieren output
            %%(vertical) along a fluid column.
            % Inputs:
            %   index   -   [1xn, int] range of indicies for image->displacement field conversion
            %   x_pos   -   [1xn, double] range of horizontal coordinates (pixel default) to average over

            if nargin<4
                isPlotting = false;
            end

            if numel(index)<=1
                error('Incorrect argument. Index must be a vector array of integers of length greater than one.')
            end
            
            im      = iw.loadimage(index(1));
            [X,Z]   = iw.getDimensions(im);
            Xmax    = size(im,2);
            
            x_pos(x_pos<0)=0;
            x_pos(x_pos>Xmax)=Xmax;
            
            index(index<1)=1;
            index(index>iw.maxIndex)=iw.maxIndex;

            loop    = index;
            t       = iw.times(loop);
            t       = seconds(t-t(1));                  
            %T       = t./(2*pi/(iw.omega*iw.N));        % wave period
            W       = zeros(size(im,1),numel(loop));

            % reference image check
            sums = zeros(1,numel(loop));

            for i=1:numel(loop)
                try
                    [tmp,imd]=iw.getImageData(loop(i),'fcd');
                catch
                    continue
                end
                sums(i) = sum(imd,[1 2]);
                W(:,i)  = mean(tmp(:,x_pos),2,'omitnan');
                displayProgress('Analysing fluid column',i,1,numel(loop));
            end

            % blank images (assuming a white background) will be associated
            % with the largest sums
            valid       = within(sums,max(sums)+std(sums)/2*[-1 1]);
            
            W(valid)    = nan;
            t(valid)    = nan;

            W2              = W-mean(W,2,'omitnan');
            t2              = t;

            if ~(all(valid))
                W2(:,valid)   = [];
                t2(:,valid)   = [];
            end
            W2 = fillmissing(W2,'linear');

            [Fk,om,m]       = ffts2(W2,t2,Z);      

            if ~nargout||isPlotting
                tiledlayout(1,3)
                n(1)=nexttile;
                imagesc(X,Z,tmp)
                line(X(x_pos(1))*[1 1],ylim,'LineStyle','--','Color','k')
                line(X(x_pos(end))*[1 1],ylim,'LineStyle','--','Color','k')
                line(min(x_pos)*[1 1],ylim,'Color','k','LineStyle','--')
                line(max(x_pos)*[1 1],ylim,'Color','k','LineStyle','--')

                n(2)=nexttile;
                imagesc(t2,Z,W2)
                addlabels('x','$t$ (s)','y','$z$ (m)')
                cmocean('balance','pivot',0)
                axis(n(1),'ij')

                n(3)=nexttile;
                imagesc(om,m,fftshow(Fk))
                addlabels('x','$\omega$ (s$^{-1}$)','y','$A$ (m$^{-1}$)')
                axis(n(2),'xy')
                clim([-15 0])
                if ~isPlotting
                    clear W Z t T Fk om m
                end
            end
        end
    end
    %% Experiment Designer integration / experiment setup
    methods
        function integrateExperimentDesign(iw,ed)
            if nargin<2
                ed = experimentdesigner('',true);
            end
            iw.ed = ed;
            iw.N  = ed.inputs.N;
        end
        function [y,t,details]=waveGenerator(iw,omegaN,varargin)
            if isempty(iw.N)
                iw.N=input('Please provide a value for N: ');
            end
            [y,t,details,export]=wavemaker(iw.N, ...
                'omega',omegaN*iw.N, ...
                'Le',1/iw.m0,...
                varargin{:});

            if export
                prg = struct('y',y,'t',t,'details',details);
                if numel(iw.programs)==1
                    iw.programs(1)=prg;
                else
                    iw.programs(end+1)=prg;
                end
            end
            if ~nargout
                clear y t details
            end
        end
    end
    %% Backlighting
    methods 
        function ss_pattern_Movie(~)
            %%SS_PATTERN_MOVIE opens a program that enables the creation of
            %%a synthetic schlieren background in video format with
            %%the potential for encoding variable/alternating backgrounds
            app_patternMovieDesigner;
        end
        function [valid,loop] = ss_findreferenceImages(iw,stride,range,silent)
            %%SS_FINDREFERENCEIMAGES opens each image (with a stride of 1
            %%by default) in the image repository and assigns the index
            %%value to the array 'refimage' depending on the background
            %%image state
            maxRange = false;
            if nargin<2
                stride = 1;
            end
            if nargin<3
                range = [1 iw.maxIndex];
                maxRange = true;
            end
            if isempty(range)
                range = [1 iw.maxIndex];
                maxRange = true;
            end
            if nargin<4
                silent = false;
            end
            range(range<1)=1;
            range(range>iw.maxIndex)=iw.maxIndex;

            loop = range(1):stride:range(end);
            sums = zeros(1,numel(loop));
            j = 1;
            for i=loop
                sums(j) = sum(iw.loadimage(i),[1 2]);
                j = j+1;
                if ~silent
                    displayProgress('Reading images',i,loop(1),loop(end))
                end
            end

            % blank images (assuming a white background) will be associated
            % with the largest sums
            valid = within(sums,max(sums)+std(sums)/2*[-1 1]);


            if maxRange
                iw.refimages = loop(valid);
            end
            if ~nargout
                clear valid
            end

        end
        function copyRefImages(iw)
            %%COPYREFIMAGES copies all the reference images into a new
            %%folder titled "allstreamsref" in the savepath directory
            if isempty(iw.refimages)
                return
            end
            allstream = fullfile(iw.savepath,'allstreamsref');
            
            if ~isfolder(allstream)
                mkdir(allstream)
            end
            
            valid = iw.files.Files(iw.refimages);
            fprintf('Copying files...')
            cellfun(@(x) copyfile(x,allstream),valid)
            fprintf('Complete\n')
        end
    end
    methods (Hidden) % Synthetic schlieren pattern generator
        function pattern = ss_pattern(~,varargin)
            %SS_PATTERN generates a regular pattern
            %('dot','lattice','checker','row','column') for use as a synthetic schlieren
            %background. 
            lambda      = 6;
            mode        = 'dot'; % or 'checker','lattice','row','column'
            isSaving    = false;
            image       = [];
            rotation    = 0;%rad2deg((pi/4)*.2);
            parseInput(varargin)
            %
            % s = diam * 72/25.4;   % diameter, in points units (1 point = 1/72 inch = 25.4/72 mm)
            %a4 @ 600dpi = [4961 7016] px = [21 29.7] cm
                       
            if isempty(image)
                x0 = 1920;%4961;
                y0 = 1080;%7016;
            else
                [x0,y0] = size(image,[1 2]);
            end
            switch mode
                case {'checker','lattice'}
                    val = ceil([y0 x0]/(lambda*2));
                    Iref = double(checkerboard(lambda,val(1),val(2))>1/2);
                    if strcmp(mode,'lattice')
                        Iref = double(abs(del2(Iref))>0);
                    end
                    pattern = Iref;
                    if rotation
                        pattern = imrotate(pattern,rotation,"bilinear",'crop');
                    end
                    pattern = pattern(1:y0,1:x0);
                    
                case 'dot'
                    pxX = [];
                    [x,y] = meshgrid(0:x0-1,0:y0-1);

                    k0 = 2*pi/(lambda-1);%*pxpercma4/pxpercm);
                    
                    xp  = cosd(rotation)*x + sind(rotation)*y;
                    yp = -sind(rotation)*x + cosd(rotation)*y;
                    Iref = .5+(cos(k0*xp) + cos(k0*yp))/4;
                    pattern = Iref;
                case 'row'
                    y = double(mod(0:y0-1,lambda)==0);
                    pattern = repmat(y,[x0 1])';
                case 'column'
                    x = double(mod(0:x0-1,lambda)==0);
                    pattern = repmat(x,[y0 1]);
            end

            if isSaving
                savepath = filepartsn(which('iwdata.m'),1);
                fname = sprintf('iw_refpattern_%ipx_%iX%iY_%s.png', ...
                    lambda,x0,y0,mode);
                imwrite(pattern,fullfile(savepath,'reference_patterns',fname),'png')
            end

            if ~nargout
                clear pattern
            end
            return
            
            %% - physical printing - %
            % desired resolution 2pi/k = 3.4px  in sensor
            % 50mm lens, 3088x2064px, 3.5m dist, ~66 px/cm
            pxpercm     = 66;
            pxpercma4   = 4961/21;
            h=figure(2);
%             set(h,'units','normalized','outerposition',[0 0 1 1])
%             clf(h,'reset')
            ax=axes(h);
            imagesc(ax,Iref')%,shg,colormap gray
            axis(ax,'image')
            set(ax,'Unit','normalized','Position',[0 0 1 1]);
            % hide the toolbar
            set(h,'menubar','none')
            % to hide the title
            set(h,'NumberTitle','off');
            
            
            colormap gray
            return
            % Window size:
            %             set(h,'Position',[360 80 560/sqrt(2) 560]);
            set(ax,'Position',[0 0 1 1]);   % location of the figure in the window
            set(ax,'PlotBoxAspectRatio',[1/sqrt(2) 1 1]);

            
            % Printing settings:
            set(h,'PaperUnits','centimeters');
            set(h,'PaperOrientation','landscape');
            set(h,'PaperType','A4');
            %set(h,'PaperSize',[21 29.7]);
            set(h,'PaperPosition',[0 0 21 29.7]);
            set(h,'PaperPositionMode','manual');
            set(h,'InvertHardcopy','off');   % keep the user background mode

            if isSaving
                fname = sprintf('%slambda_%.2gpx.tiff','..\',lambda);%*pxpercma4/pxpercm));
                fprintf('Saving \n%s\n',fname)
                %                 print('-dbmp','-r1000',fname);
                %print('-dtiff','-r2000',fname)
                %out2latex(fname,'r',1600,'type','tiff')
                imwrite(Iref,fname,"tiff","Compression","none","Resolution",100)
                fprintf('Done\n')
            end

            %% Input parser
            function parseInput(varargin)
                m = 1;
                items = varargin{:};
                for k=1:length(items)
                    switch items{m}
                        case {'lambda','wavelength'}
                            lambda = namevalue;
                        case {'rotation','rot'}
                            rotation = namevalue;
                        case 'ref'
                            image = namevalue;
                            if any(numel(image)==[2 4])
                                % Assume screesize argument
                                image = image(image~=1);
                                image = zeros(image(1),image(2));
                            end
                        case {'dot','checker','lattice','row','column'}
                            mode = items{m};
                        case 'save'
                            isSaving = true;
                        case {'resolution','-r'}
                            camRes = namevalue;
                        case {'size'}
                            paperSize = namevalue;
                        case 'density'
                            dotDensity = namevalue;
                        case 'dpi'
                            dpi = namevalue;
                        case ''
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
        function [B,Idef] = test_ss_pattern(iw,Iref,b,warpfactor)
            if nargin<4
                warpfactor = 1;
            end
            [fx,fz] = gradient(b);
            Idef    = WarpImage(Iref,fx*warpfactor,fz*warpfactor);

            % noise
            Idef = imnoise(Idef,"gaussian",1,.5);

            [u,w,B] = iw.ss_solver_fcd(Iref,Idef);

            if nargout
                return
            end
            
            tiledlayout(3,2)
            nexttile
            imagesc(Iref)
            title('Reference pattern')
            nexttile
            imagesc(Iref-Idef)
            title('Difference')

            
            nexttile
            imagesc(b)
            title('Input buoyancy field')
            colorbar
            nexttile
            imagesc(B)
            title('Reconstructed buoyancy field')
            subtitle(sprintf('%.2f',1/mean(b(:,:,1)./B,'all')))
            colorbar

            nexttile([1 2])
            imagesc(log10(abs(b./B)))
            title('log_{10}(|Input/Output|)')
            colorbar
            clear B
        end
    end
    %% NetCDF rendering
    methods
        function state = nc_create(iw,reset)
            %% Creates a NetCDF file in the savepath repository with the appropriate label
            % and generates the following variables:
            %   1. im-field (image)
            %   2. u-field  (horizontal velocity)
            %   3. v-field  (vertical velocity)
            %   4. b-field  (buoyancy)
            %   5. x-coord  (native)
            %   6. y-coord  (native)
            %   7. t-coord  (native)
            %   8. ds       (image resolution)
            %   9. dt       (time resolution)
            %   10. resize  (native resizing factor)
            %   11. rend    (rendered status)
            
            if nargin<2
                reset = false;
            end
            

            [~,fname] = nc_check(iw);

            % native dimensions
            [y,x]   = size(iw.im0);
            x2      = round(x/2);
            y2      = round(y/2);
            t       = iw.maxIndex;

            if reset
                fprintf('Deleting\n\t%s\n',fname)
                delete(fname)
            end

            % check for file
            if ~isfile(fname)||reset
%                 nccreate(fname,'im','Dimensions',{'y',y,'x',x,'t',Inf},'FillValue',nan)
                nccreate(fname,'u','Dimensions',{'y',y2,'x',x2,'t',Inf}, ...
                    'FillValue',nan,'Datatype','single')
                nccreate(fname,'v','Dimensions',{'y',y2,'x',x2,'t',Inf}, ...
                    'FillValue',nan,'Datatype','single')
%                 nccreate(fname,'b','Dimensions',{'y',y,'x',x,'t',Inf},'FillValue',nan)

                nccreate(fname,'t','Dimensions',{'t',t},'FillValue',nan)
                nccreate(fname,'rend','Dimensions',{'t',t},'FillValue',0)
                
                nccreate(fname,'ds','Dimensions',{'val',1})
                nccreate(fname,'dt','Dimensions',{'val',1});
                nccreate(fname,'resize','Dimensions',{'val',1});
                ncwrite(fname,'dt',milliseconds(mean(diff(iw.times))))
                ncwrite(fname,'ds',iw.ds)
                ncwrite(fname,'resize',.5)
                iw.hasNCFile = true;
            else
                iw.hasNCFile = true;
            end
            state = iw.hasNCFile;
            if ~nargout
                fprintf('Return state: %s\n',string(state))
                clear state
            end
        end
        function [state,fname] = nc_check(iw)
            try
                [~,str] = iw.savedLabel;
            catch
                state = false;
                fname = iw.savepath;
                return
            end
            str     = strcat('outputs_',str,'.nc');
            fname   = fullfile(iw.savepath,str);
            iw.ncfilepath = fname;
            state = isfile(fname);
        end
        function nc_render(iw,index,data,field)
            if ~iw.hasNCFile
                if ~iw.nc_create
                    error('Unable to render field. Cannot create netCDF file.')
                end
            end

            % resize 2D data to compressed dimensions
            ni  = ncinfo(iw.ncfilepath);
            s   = ni.Variables(strcmp({ni.Variables.Name},'u')).Size([1 2]);

            if any(strcmp(field,{'im','u','v','b'}))
                m=imresize(iw.ss_fcd.mask,s);
                data(~m)=nan;
                if all(size(data)~=s)
                    data = imresize(data,s);
                end
            end

            % write time-series data
            ncwrite(iw.ncfilepath,field,data,[1 1 index],[1 1 1])
            % write time index data
            ncwrite(iw.ncfilepath,'t',index,index)
            % write rendered status data
            ncwrite(iw.ncfilepath,'rend',1,index)
        end
        function [output,outputstruct] = nc_read(iw,index,field,isCheck)

            if nargin<4
                isCheck = false;
            end

            data1d = {'x','z','t','rend','ds','dt'};
            % handle 1D column data
            if nargin==2
                field = index;
                if any(strcmp(field,data1d))
                    output = ncread(iw.ncfilepath,field);
                    return
                end
                error('Not enough input arguments')
            end

            % get info about field dimensions
            ni      = ncinfo(iw.ncfilepath);

            Ustate = false;
            if strcmp(field,'U')
                field = 'u';
                Ustate = true;
            end
            target  = ni.Variables(strcmp({ni.Variables.Name},field));
            dim     = target.Size;
            fv      = target.FillValue;

            % handle non-unique indicies
            stride = 1;
            if numel(index)>1
                stride = round(mean(diff(index)));
            end

            switch numel(dim)
                case 1
                    indexing    = index(1);
                    startidx    = numel(index);
                    stride      = 1;
                case 2

                case 3
                    indexing = [1 1 index(1)];
                    startidx = [dim(1:2) numel(index)];
                    stride   = [1 1 stride];
            end

            try
                output  = ncread(iw.ncfilepath,field,indexing,startidx,stride);
                if Ustate
                    v  = ncread(iw.ncfilepath,'v',indexing,startidx,stride);
                    output = hypot(output,v);
                end
            catch
                output = nan;
            end
            output(output==fv)=nan;
            valid =  ~all(isnan(output),"all");

            if isCheck 
                output = valid;
            end

            outputstruct = [];

            if ~valid
                return
            end
            % build ss_fcd structure
            if Ustate
                iw.ss_fcd.U = output;
                return
            end
            if strcmp(field,'u')
                iw.ss_fcd.u = output;
            end
            if strcmp(field,'v')
                iw.ss_fcd.v = output;
            end
        end
    end
    % Batch rendering
    methods
        function renderImage(iw,index)
             getImageData(iw,index,'fcd','u',0,false);
%              iw.nc_render(index,iw.ss_fcd.im,'im')
             iw.nc_render(index,iw.ss_fcd.u,'u')
             iw.nc_render(index,iw.ss_fcd.v,'v')
        end
        function renderBatch(iw,start,step,stop,useParallel,waitbar_handle)
            if nargin<2;start = 1;end
            if nargin<3;step = 1;end
            if nargin<4;stop = iw.maxIndex;end
            if nargin<5;useParallel=gcp('nocreate');end
            if nargin<6;waitbar_handle=[];end


            if ischar(waitbar_handle)
                if any(strcmp(waitbar_handle,{'waitbar','wb'}))
                    waitbar_handle = waitbar(0,'Preallocating memory',...
                        'Name','Rendering experiment');
                end
            end
            loop = start:step:stop;
            
            flag = zeros(1,numel(loop));
            LOOP = 1:numel(loop);
            fprintf('Rendering\n')

            try
            if ~useParallel
                fprintf('Serial method\n')
                if isempty(waitbar_handle)
                    w = [];%waitbar(0,'Processing','Name','Serial computation');
                else
                    w = waitbar_handle;
                    try 
                        w.Title = 'Serial computation';
                    end
                end
                displayProgressparfor(w,numel(loop))
                for ii=LOOP
                    try
                        iw.renderImage(loop(ii))
                    catch ME
                        flag(ii) = true;
                    end
                    
                    displayProgressparfor(w)
                end
            else
                fprintf('Parallel method\n')
                q = parallel.pool.DataQueue;
                parfevalOnAll(@warning, 0,'off','all');
                afterEach(q,@displayProgressparfor)
                if isempty(waitbar_handle)
                    w = [];
                else
                    w = waitbar_handle;
                end
                displayProgressparfor(w,numel(loop))
                parfor ii=LOOP
                    try
                        iw.renderImage(loop(ii))
                        send(q,[])
                    catch
                        send(q,[])
                        continue
                    end
                end
            end
            catch ME
                close(waitbar_handle)
                rethrow(ME)
            end
        end
    end
    %% LabExperiment integration
    methods
        function loadLabExperimentData(iw)
            le = dir(fullfile(iw.savepath,'*experimentdata.mat'));
            if isempty(le)
                return
            end
            iw.experiment = load(fullfile(le.folder,le.name)).le;
            iw.N    = iw.experiment.N;
            iw.omega= unique(iw.N*iw.experiment.parameterspace.omegaN);
            iw.date = iw.experiment.dates;
        end
    end
    %% Saving / dating / labelling
    methods
        function [saveto,str,lbl]=savedLabel(iw)
            [root,imagefolder] = fileparts(iw.files.Folders);
            if isempty(imagefolder)
                [root,imagefolder] = fileparts(fileparts(iw.files.Files{1}));
            end
            if iscell(root)
                iw.savepath = unique(root);
            else
                iw.savepath = root;
            end
            str = char(join(imagefolder,'+'));
            lbl = sprintf('iwdata_%s.mat',str);
            saveto = fullfile(iw.savepath,lbl);%fullfile(iw.savepath,iw.label);
        end
        function saveDesign(iw)
            saveto=savedLabel(iw);
            if ~isfolder(fileparts(saveto))
                mkdir(fileparts(saveto))
            end
            if isfile(saveto)
                shortWarning(sprintf('\nOverwriting %s',char(saveto)))
            end
            
            % clear image outputs and non-essential data
            iw.resetOutputs

            save(char(saveto),'iw')
            iw.hasSavedFile = true;
        end
        function resetOutputs(iw)
            % Resets all non-essential properties to their default state
            % (for smaller saved .mat files)
            mask = iw.ss_fcd.mask;
            iw.ss_fcd       = struct('mask',mask);
            iw.lat          = struct;
            iw.im0          = [];
            iw.history      = struct;
            iw.VideoData    = struct;
        end
        function changeDate(iw,to)
            if nargin<2
                to = datetime("today");
            end
            iw.date     = to;
            oldlabel    = iw.label;
            iw.changeLabel(extractAfter(iw.label,'\'))
            iw.updateSavedFile(oldlabel)
        end
        function changeLabel(iw,to)
            if iw.internalCommand
                iw.internalCommand = false;
                return
            end
            oldlabel = iw.label;
            frontlabel = string(iw.date(1),'uuuu_MM_d\');
            if nargin<2
                to = 'untitled';
            end

            iw.internalCommand = true;
            iw.label = strcat(frontlabel,to);
            iw.updateSavedFile(oldlabel)
        end
        function updateSavedFile(iw,oldlabel)
            if ~iw.hasSavedFile
                return
            end
            from= fullfile(iw.savepath,strcat(oldlabel,'.mat'));
            to  = fullfile(iw.savepath,strcat(iw.label,'.mat'));
            shortWarning(sprintf('Changing saved filename from\n\t%s\n\tto\n\t%s',from,to))
            try
                movefile(from,to)
            catch ME
                if strcmp(ME.identifier,'MATLAB:MOVEFILE:ResourceNotFound')
                    mkdir(fileparts(to))
                    movefile(from,to)
                    rmdir(fileparts(from))
                end
            end
            %             if isempty(iw.oldpaths)
            iw.oldpaths{end+1}=from;
        end
    end
    methods (Hidden) % Configuration functions and initialization
        function initialize(iw) % configure defaults upon construction of class
            %% Synthetic schlieren solver (Digital image correlation)
            ss_dic = struct; %#ok<*PROP>
            ss_dic.u = [];
            ss_dic.v = [];
            ss_dic.controls = struct(...
                'fast',true,...
                'bs',4,...
                'bo',4,...
                'md',4,...
                'mq',0.01,...
                'nc',true,...
                'im0index',1,...
                'resize',50,...
                'roi',[], ...
                'fcdrolling',0,...
                'fcdrollingstep',1,...
                'rotation',0);

            %% Synthetic schlieren solver (Fast Checkerboard Demodulation)
            ss_fcd = struct;
            ss_fcd.u = [];
            ss_fcd.v = [];
            ss_fcd.mask = [];
            ss_fcd.fail = false;
            %% Configure
            iw.ss_dic = ss_dic;
            iw.ss_fcd = ss_fcd;
        end
        function runSetup(~) % interactive designer / tank+wave emulator
            %% Figure
            XLIM = 1.5;
            YLIM = 1.1;

            fig = uifigure('Name','Internal Wave Experiment Modeller',...
                'Resize','on');
            fig.Position = [fig.Position(1:2) 750 420];
            gl = uigridlayout('Parent',fig,'ColumnWidth',{'1x'},'RowHeight',{'1x',100});
            %% Axes
            ax = uiaxes('Parent',gl,...
                'XLim',[0 XLIM],...
                'YLim',[0 YLIM],...
                'Box','on');
            ax.Toolbar.Visible='off';
            ax.XLabel.String = 'Position (m)';
            ax.YLabel.String = 'Depth (m)';
            set(ax,'YDir','reverse')

            %% iwects
            ice     = drawrectangle('Position',[XLIM-.2 0 .2 YLIM],'Tag','ice','Parent',ax);
            iceth   = imdistline(ax,[ice.Position(1) sum(ice.Position([1 3]))],YLIM/3.*[1 1]);
            source  = drawline('Position',[0 YLIM/2;ice.Position(1)/2 YLIM/4],...
                'Label','',...
                'LabelAlpha',0,'Tag','source','Parent',ax);
            angle       = @(source) abs(rad2deg(atan(diff(source.Position(:,2))/diff(source.Position(:,1)))));
            radangle    = @(source) atan(diff(source.Position(:,2))/diff(source.Position(:,1)));

            %% Wavemaker
            h0  = 5e-2;
            wavemaker = [];
            %% Buttons
            gl2     = uigridlayout(gl,'ColumnWidth',repmat({80},[1 8]),'RowHeight',{30,30});
            Nbutton = uibutton(gl2,'Text','Get N');
            uilabel(gl2,'Text','N (rad/s)');
            N       = uieditfield(gl2,'numeric','Value',.3,'ValueChangedFcn',@updateangle);
            uilabel(gl2,'Text','Omega (rad/s)');
            omega   = uieditfield(gl2,'numeric','Value',N.Value*sin(deg2rad(angle(source))));
            uilabel(gl2,'Text','Fr');
            fr      = uieditfield(gl2,'numeric','Value',1);
            sim     = uibutton(gl2,'state','Value',false,'Text','Simulate','ValueChangedFcn',@runModel);
            uilabel(gl2,'Text','');
            uilabel(gl2,'Text','Ray tracing');
            rays2render = uispinner(gl2,'Limits',[0 10],'Value',1,'ValueChangedFcn',@drawRays);

            % copy axes data buttons
            uibutton(gl2,'Text','Copy axes','ButtonPushedFcn',@copyProtocol,'Tag','copy0');
            uibutton(gl2,'Text','Newtile','ButtonPushedFcn',@copyProtocol,'Tag','copy1','Visible','off');
            %% Rays
            negsource           = drawline(ax,'Position',source.Position,'Color',...
                'k','Tag','negsource',...
                'InteractionsAllowed','none');
            negsource.Position  = source.Position-[0 0;0 2*diff(source.Position(:,2))];
            for i=1:rays2render.Value
                rays_inline(i)  = drawline(ax,'Position',source.Position,'Color','k',...
                    'Tag',sprintf('inlinerays-%i',i),...
                    'InteractionsAllowed','none','Visible','off'); %#ok<*AGROW>
                rays_counter(i)  = drawline(ax,'Position',source.Position,'Color','k',...
                    'Tag',sprintf('counterrays-%i',i),...
                    'InteractionsAllowed','none','Visible','off');
            end
            %% Listeners
            addlistener(ice,'MovingROI',@(varargin)lineChangingFcn(ice));
            addlistener(source,'MovingROI',@(varargin)lineChangingFcn(source));
            %% Initiate
            runModel
            addColorbar('ax',ax,'cmap','balance','pivot',0)
            %             colorbar(ax)
            %% Functions
            function copyProtocol(evt,src)
                switch src.Source.Tag
                    case 'copy0'
                        shg,cla
                        tmpax   = gca;
                        tmpfig  = gcf;
                        [];
                    case 'copy1'

                end

            end
            function lineChangingFcn(iwect)
                switch iwect.Tag
                    case 'ice'
                        iceth.setPosition([ice.Position(1) YLIM/3;sum(ice.Position([1 3])) YLIM/3]);
                    case 'source'
                        iwect.Position(1,1)=0;
                end
                if source.Position(2,1)>=ice.Position(1)
                    source.Position(2,1)=ice.Position(1);
                end
                runModel
                updateWavemaker
                updateAngle
                updateRays

                function updateRays
                    clc
                    negsource.Position = source.Position-[0 0;0 2*diff(source.Position(:,2))];
                    rayTrace(source,rays_inline(1),'project')
                    rayTrace(negsource,rays_counter(1),'project')
                    if rays2render.Value==0
                        return
                    end
                    ch = get(ax,'Children');
                    delete(findall(ch,'Tag','intercept'))
                    for j=2:2*rays2render.Value-1
                        clr = 1-[1 1 1]*exp(-ceil(j/2)^2/(5*rays2render.Value));
                        rayTrace(rays_inline(j-1),rays_inline(j),'reflect',clr)
                        rayTrace(rays_counter(j-1),rays_counter(j),'reflect',clr)
                    end

                end
                function updateWavemaker
                    if isempty(wavemaker)
                        return
                    end
                    if ~sim.Value
                        declutterAxes('Patch')
                    end
                end
                function updateAngle
                    omega.Value = N.Value*sin(deg2rad(angle(source)));
                    source.Label = sprintf('%.1f deg (%.2fN)\nh=%.2f',...
                        rad2deg(radangle(source)),...
                        omega.Value/N.Value,...
                        source.Position(1,2));

                end
            end
            function drawRays(~,~)
                child   = get(ax,'Children');
                islines = cell2mat(arrayfun(@(x) contains(class(x),'Line'),child,'UniformOutput',false));
                arelines = child(islines);
                omitrays= cell2mat(arrayfun(@(x) contains(x.Tag,'ray'),arelines,'UniformOutput',false));
                delete(arelines(omitrays))
                count           = 2*rays2render.Value-1;
                if rays2render.Value==0
                    return
                end
                rays_inline     = giwects(1,count);
                rays_counter    = giwects(1,count);
                for k=1:count
                    rays_inline(k)  = drawline(ax,'Position',source.Position,'Color','k',...
                        'Tag',sprintf('inlinerays-%i',k),...
                        'InteractionsAllowed','none','Visible','off'); %#ok<*AGROW>
                    rays_counter(k)  = drawline(ax,'Position',source.Position,'Color','k',...
                        'Tag',sprintf('counterrays-%i',k),...
                        'InteractionsAllowed','none','Visible','off');
                end
                lineChangingFcn(source)
            end
            function runModel(~,~)
                %% Shield
                if ~sim.Value
                    declutterAxes('Image')
                    return
                end
                %% Parse model
                res = 1.5e2;
                X   = linspace(0,XLIM,res);
                Z   = linspace(0,YLIM,res);
                x0  = source.Position(1,2);
                H   = XLIM-iceth.getDistance;
                k0  = iw.k0;
                omegatransp = N.Value*cos(deg2rad(angle(source)));
                % Fix X
                x   = X;
                lim = find(x>H,1,'first');
                xx  = x(1:lim);
                % Fold Z to avoid periodic solution
                zfolds = 10;

                n_var = floor(N.Value/omega.Value);
                if n_var>3
                    n_var = 3;
                end

                [b,wavemaker] = iw_iterate('t0',...
                    'N',N.Value,...
                    'omega',omegatransp,...
                    'H',H,...
                    'L',zfolds.*[-YLIM YLIM],...
                    'k0',k0,...
                    'h0',h0,...
                    'x0',x0,...
                    'ut',fr.Value*N.Value/k0,...
                    'n',n_var,...
                    'nu',5e-5,...
                    'alpha0',2,...
                    'res',1/(zfolds*res+1),...
                    'zres',lim,...
                    'output','h');
                %% Folding z-axis
                zz          = linspace(zfolds*-YLIM,zfolds*YLIM,(zfolds*res));
                [~,regions] = unique(floor(zz/YLIM));
                regions(end)= regions(end)+1;
                regions     = regions(1:end);
                flippingstate = [0; repmat([0;1],[floor(numel(regions)/2) 1])];
                flippingstate = circshift(flippingstate,ceil(numel(regions)/2));

                B           = zeros([size(b,1) mean(diff(regions))+1 (numel(regions)-1)]);
                for ii=1:numel(regions)-1
                    space   = regions(ii):regions(ii+1);
                    in      = b(:,space);
                    if flippingstate(ii)
                        in  = fliplr(in);
                    end
                    B(:,:,ii) = in;
                end
                b           = fliplr(sum(B,3));
                %% Plotting
                declutterAxes({'Image','Patch'})
                hold(ax,'on')
                wavemaker = wavemaker(zz>0&zz<(YLIM));
                wavemaker([1 end])=0;
                fill(ax,wavemaker,zz(zz>0&zz<(YLIM)),'w')
                imagesc(ax,xx,Z,rot90(b,-1))
                hold(ax,'off')
                caxis(ax,[-1 1].*1e-3)
                set(ax,'Children',flipud(get(ax,'Children')))
            end
            function rayTrace(src,iwect,mode,clr)
                %%
                warning off
                iwect.Color = 'r';
                if nargin<4
                    clr = 'k';
                end
                findy   = @(x,dydx,origin) x*dydx+origin;
                findx   = @(y,dydx,origin) (y-origin)/dydx;
                switch mode
                    case 'project'
                        isprojection = true;
                    case 'reflect'
                        isprojection = false;
                end
                try
                    tmp = src.Position;
                    %                     if tmp
                    iwect.Position = calculateNewEndPoint(src,isprojection);
                    intercept = iwect.Position(:,1)==ice.Position(1);
                    if any(intercept)
                        line([0 ice.Position(1)],iwect.Position(intercept,2).*[1 1], ...
                            'Parent',ax,...
                            'LineStyle','--','Color','k','Tag','intercept')
                    end
                catch
                    %                     iwect.Position = calculateNewEndPoint(src,isprojection);
                    [];
                end
                iwect.Visible  = 'on';
                iwect.Color    = clr;
                %%
                function newpos = calculateNewEndPoint(src,isprojection)
                    grad    = diff(src.Position);
                    slope   = grad(2)/grad(1);
                    origin  = src.Position(1,:);
                    endpt   = src.Position(2,:);
                    tgt     = ice.Position(1);

                    if nargin<2
                        isprojection = false;
                    end
                    intercept = [tgt findy(tgt,slope,origin(2))];

                    if isprojection
                        if ~within(intercept(2),[0 YLIM])
                            ypt = YLIM;
                            if intercept(2)<0
                                ypt = 0;
                            end
                            intercept(1) = findx(ypt,slope,origin(2));
                            intercept(2) = ypt;
                        end
                        newpos = [endpt;intercept];
                        return
                    end

                    %% Incidence angle
                    if within(abs(radangle(src)),[0 pi])
                        quadrant = 1;
                    else
                        quadrant = 2;
                    end
                    if sign(radangle(src))<0
                        quadrant = 5-quadrant;
                    end

                    %% Boundary
                    ishorzsidewall = any(src.Position(:,1)==ice.Position(1)|src.Position(:,1)==0);

                    %% Reflection condition
                    slope   = -slope;

                    if ishorzsidewall
                        if any(quadrant==[1 4])
                            tgt = 0;
                        end
                    end

                    intercept = [tgt findy(tgt-endpt(1),slope,endpt(2))];
                    if ~within(intercept(2),[0 YLIM])
                        ypt = YLIM;
                        if intercept(2)<0
                            ypt = 0;
                        end
                        intercept(1) = findx(ypt,slope,endpt(2))+endpt(1);
                        intercept(2) = ypt;
                    end
                    %% Boundary intercept

                    newpos = [endpt;intercept];
                    if contains(src.Tag,'3')
                        [];
                    end
                end
            end
            function declutterAxes(classes)
                child = get(ax,'Children');
                omit = cell2mat(arrayfun(@(x) contains(class(x),classes),child,'UniformOutput',false));
                delete(child(omit))
            end
        end
    end

    %% Set functions and associated routines
    methods
        function set.files(iw,val)
            iw.files = val;
            manageFiles(iw)
        end
        function set.label(iw,val)
            changeLabel(iw,val)
        end
    end 
    methods 
        function manageFiles(iw)
            if iw.internalCommand
                return
            end
            loadLabExperimentData(iw)
            updateDatastore(iw)
            updateTimes(iw)
            updateSavePath(iw)
            updateLabel(iw)
        end
    end
    methods (Hidden)
        function updateTimes(iw,type)
            if iw.internalCommand
                return
            end
            if nargin<2
                type = 'date';
            end
            f   = iw.files;
            if isempty(f)
                return
            end
            switch class(iw.files)
                case 'imMerge'
                    % files should be sorted
                    iw.times = mean([f.dates1(1:f.maxIndex); f.dates2(1:f.maxIndex)]);
                    return
            end

             switch type
                 case 'index'
                    sortByIndex(iw)
                 case 'date'
                     sortByDate(iw)
             end
        end
        function sortByDate(iw)
            % files unsorted
            fold = iw.files.Folders;
            if isempty(fold)
                fold = fileparts(iw.files.Files{1});
            end
            tmp = arrayfun(@dir,(fullfile(fold,strcat('*.',iw.files.SupportedOutputFormats))),'UniformOutput',false);
            tmp = cell2mat(nonempty(tmp));
            tmp = tmp(~[tmp.isdir]);
            dt  = datetime({tmp.date});
            [dt,map] = sort(dt);
            %[~,fnames]=fileparts(iw.files.Files);
            %                     [~,map] = sort(cellfun(@(f) str2double(f(end-3:end)),fnames));
            iw.files.Files=iw.files.Files(map);
            iw.times = dt;
        end
        function sortByIndex(iw)
            file                = extractBefore(iw.files.Files,'.');
            try
                [~,files]=fileparts(file);
                [~,map]=sort(cellfun(@str2double,extractBetween(files,'C','T')));
                if isempty(map)
                    sortByDate(iw)
                    return
                end
            catch
                sortByDate(iw)
                return
                %[~,map]             = sort(cellfun(@(f) str2double(f(end-3:end)),file));
            end
            iw.files.Files      = iw.files.Files(map);
        end
        function validateDatastore(iw,imagePath)
            if nargin<2
                imagePath = iw.imagepath;
            end


            fold            = iw.files.Folders;
            if isempty(fold)
                fold = fileparts(iw.files.Files{1});
            end

            [root,fname]    = fileparts(fold);
            [root0,root1]   = fileparts(imagePath);

            if ~any(strcmp(root1,fname))
                root0 = imagePath;
            end

            iw.imagepath = imagePath;

            
            if all(strcmp(root0,root))
                return
            end
            
            newroot = fullfile(root0,fname);
            if all(isfolder(newroot))
                try
                    iw.files.Folders = newroot;
                catch
                    iw.files = imageDatastore(newroot);
                end
                laststate = iw.internalCommand;
                iw.internalCommand = false;
                updateDatastore(iw)
                iw.internalCommand = laststate;
                return
            end
        end
        function updateDatastore(iw)
            if iw.internalCommand
                return
            end
            fold = iw.files.Folders;

%             % check for valid folders (if source was moved)
%             valid = ~isfolder(fold);
%             if any(valid)
%                 validateDataStore(iw)
%             end

            if isempty(fold)
                fold = fileparts(iw.files.Files{1});
            end
            iw.internalCommand = true;
            switch class(iw.files)
                case 'imMerge'
                    try
                        ms = readstruct(fullfile(iw.fileLocation,'cache','mergeImageDefault.xml'));                        
                    end
                    iw.files = imMerge( ...
                        imageDatastore(fold{1}), ...
                        imageDatastore(fold{2}));
                    iw.files.mergeSettings = ms;
                otherwise
                    iw.files = imageDatastore(fold);
            end
            iw.internalCommand = false;

            try 
                iw.maxIndex = iw.files.maxIndex;
            catch
                iw.maxIndex = numel(iw.files.Files);
            end

            try
                iw.camValue = str2double(unique(extractBetween(iw.files.Files,'cam','_')));
            end
        end
        function updateSavePath(iw)
            if iw.internalCommand
                return
            end
            fold = iw.files.Folders;
            if isempty(fold)
                fold = fileparts(iw.files.Files{1});
            end
            fname = fileparts(fold);
            if isa(fname,'cell')
                fname = char(unique(fname));
            end
            iw.savepath = fname;
        end
        function updateLabel(iw)
            % Looks for file encoding
            fmt = iw.fnameformat;
            tmp = split(fmt,{'<','>'});
            tmp = nonempty(tmp(~contains(tmp,extractBetween(fmt,'<','>'))));
            f   = iw.files.Files;
            f   = unique(extractBetween(f,'cam','C'));
            if isempty(f)
                return
            end

            % keyword tags
            tags = {'omegaN','amp'};
            if ~contains(f{1},tags)
                return
            end
            iw.label = string(join(f,'+'));
        end
    end
    %% File handing
    methods
        function str = file_encodeInfo(~,omegaN,amp)
            omegaNs  = '';
            amps     = '';
            try
                if isa(omegaN,"double")
                    omegaN = num2str(omegaN);
                end
                omegaNs  = strrep(sprintf('%s',omegaN),'.','_');
            end
            try
                if isa(amp,"double")
                    amp = num2str(amp);
                end
                amps     = strrep(sprintf('%s',amp),'.','_');
            end
            str = sprintf('%s%s%s%s-','omegaN',omegaNs,'amp',amps);
        end
        function file_rename(iw,tags,values,range)
            %%FILE_RENAME renames all files in the imagedatastore that fit
            %%the user defined criteria. 
            % Arguements:
            %   tags    -   [cell,char] 
            %               {'omegaN','amp'} 
            %               or 'auto' (no values given, attempts to rewrite
            %               filenames if duplicates are present) 
            %   values  -   [array,double]
            %               [tag1,tag2]
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Example 1: 
            %   Append filenames to include missing information
            % filename = 'cam1_C1T0'
            % file_rename(iw,{'omegaN','amp'},[0.5,10],'all')
            % >> filename = 'cam1_omegaN0_5amp10C1T0'
            %
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Example 2: 
            %   Rewrite filenames to include missing information
            % filename = 'cam1_omegaN0_5amp10C1T0'
            % file_rename(iw,{'omegaN','amp'},[0.1,30],'all')
            % >> filename = 'cam1_omegaN0_1amp30C1T0'
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Guard-blocks
            isLoop = false;
            if nargin<2
                error('No tag arguments entered')
            end
            if nargin<3&&~contains(tags,'auto')
                error('No value arguments entered')
            end
            if nargin<4
                range = 1:iw.maxIndex;
            end
            if ~isa(tags,'cell')
                tags = {tags};
            end
            try
                if isa(range,'char')
                    if contains(range,'all')
                        range = 1:iw.maxIndex;
                    elseif contains(range,{'i','index','loop'})
                       isLoop = true;
                    else
                        error('Unknown range')
                    end
                end
            end

            if isLoop
                % loop through each index
                range = 1:iw.maxIndex;
                arrayfun(@(i) iw.file_rename(tags,values,i),range)
                home
                return
            end

            isAuto = strcmp(tags,'auto');

            f   = string(iw.files.Files(range));
            f0  = f;

            % Logical checks for expected format
            % cam<id|int>_omegaN<val>amp<val>C<count|int>T<time>.tiff
            [root,fnames]   = fileparts(f0);
            hasBasicFormat  = contains(fnames,{'cam','C','T'});
            hasInfoFmt      = contains(fnames,{'omegaN','amp'});
            hasRepeats      = [count(fnames,'omegaN')>1 count(fnames,'amp')>1];
            
            if any(~hasBasicFormat)
                fprintf('Adding basic info\n')
                Fnew = f(~hasBasicFormat);
                c   = 1:numel(Fnew);
                t   = arrayfun(@dir,Fnew);
                t   = datetime({t.date},'Format','uuuuMMdd''T''HHmmss');
                f   = arrayfun(@(idx) fullfile(root(1),sprintf('cam_C%i-%s',c(idx),t(idx))),1:numel(c))';
            end

            cam = @(f) string(cell2mat(arrayfun(@(camid) sprintf('cam%s',camid),extractBetween(f,'cam','_'),'UniformOutput',false)));

            if any(~hasInfoFmt)
                fprintf('Adding experiment info\n')
                [~,Fnew] = fileparts(f(~hasInfoFmt));
                Fnew = arrayfun(@(Fnew,cam) strrep(Fnew,cam,sprintf('%s%s',cam,getInfo(f))),Fnew,cam(f));
                f = fullfile(root(1),Fnew);
            end

            if any(hasRepeats,"all")
                fprintf('Removing repeated experiment info\n')
                Fnew0   = f(any(hasRepeats,2));
                [info,exptinfo] = getInfo(Fnew0);
                f       = replace(Fnew0,exptinfo,info);
            end
            
            f = replace(f,extractBetween(f,cam(f),'C'),getInfo(f));

            isSame = strcmp(f,f0);

            f0  = f0(~isSame);
            f   = f(~isSame);

            if isempty(f)
                fprintf('No changes\n')
                return
            end
            
            renameProcedure(iw,f0,f)

            function renameProcedure(iw,f,fnew)
                fprintf('Renaming files...')
                arrayfun(@(f,fnew) movefile(f,fnew),f,fnew)
                fprintf(' done\n')
                iw.updateDatastore
            end
            function [info,exptinfo] = getInfo(f)
                exptinfo    = extractBetween(f,cam(f),'C');
                Fnew        = replace(exptinfo,{'-'},'');
                if isAuto
                    omegaV  = join(nonempty(split(unique(extractBetween(Fnew,'omegaN','amp')),'_')),'.');
                    ampV    = unique(replace(extractBefore(extractAfter(Fnew,'amp'),'omegaN'),'_',''));
                    if ismissing(ampV)
                        ampV = replace(extractAfter(Fnew,'amp'),{'_','-'},'');
                    end
                    omegaV  = str2double(omegaV);
                    ampV    = str2double(ampV);
                    info    = ['_' iw.file_encodeInfo(omegaV,ampV)];
                else
                    info    = ['_' iw.file_encodeInfo( ...
                        values(strcmp(tags,'omegaN')), ...
                        values(strcmp(tags,'amp')))];
                end
            end
        end
    end
    %% IWDATA ARRAY functions
    % For use if multiple iwdata classes are called
    methods
        function getReferenceImages(iw)
            %%GETREFERENCEIMAGES searches for the blank background images
            %%in a given 'iwdata' class and copies them to a folder titled
            %%'allstreamsref'
            arrayfun(@(x) x.manageFiles,iw)
            tmp=iw([iw.camValue]==2);
            k = 1;
            arrayfun(@fcn,tmp)
            
            % Delete reference images in folder
            allstream = fullfile(iw.savepath,'allstreamsref');
            f = dir(allstream);
            if ~isempty(f)
                fprintf('Clearing files')
                f = fullfile(allstream,{f.name});
                delete(f)
            end

            arrayfun(@(x) x.copyRefImages,tmp)

            function fcn(selection)
                selection.ss_findreferenceImages(10)
                k = k+1;
            end
        end
    end

    %% ICEDATA functions
    methods
        function getBulkMeltingRate(iw,id)
            %%GETBULKMELTINGRATE locates the reference images near the
            %%beginning and end of the available time series and attempts
            %%to calculate the melting rate from the edge-detection
            %%algorithms nested in the 'icedata' class
            if numel(iw)>1
                iwIce=iw([iw.camValue]==2);
                figure
                dockfig
                tiledlayout(1,numel(iwIce))
                for i=1:numel(iwIce)
                    nexttile
                    iwIce(i).ss_dic.controls.rotation = -2;
                    iwIce(i).ss_dic.controls.roi=[254 79 514 2597];
                    try
                        iwIce(i).getBulkMeltingRate
                    catch
                        continue
                    end
                end
                return
            end
            f = iw.files.Files;
            if isempty(iw.refimages)
                [vS,lS] = iw.ss_findreferenceImages(10,[1 100],true);

                % last wave forcing index
                lw = find(contains(f,'amp-0'),1);
                if isempty(lw)
                    lw = iw.maxIndex;
                end
                [vE,lE] = iw.ss_findreferenceImages(10,lw-[100 0],true);
                index   = [(lS(vS)) (lE(vE))];
            end
            ims     = iw.loadimage(index);
            validS  = ismember(index,(lS(vS)));
            ims     = cat(3,mean(ims(:,:,validS),3),mean(ims(:,:,~validS),3));
            dt      = diff(iw.times([index(1) index(end)]));
            
            lbl = extractBetween(f(min(index)),sprintf('cam%i_',iw.camValue),'_C');
            
            imagesc(imfuse(ims(:,:,1),ims(:,:,2)))
            title(string(dt))
            subtitle(lbl,'interpreter','none')
            axis image
            return
            if nargin<2
                id = icedata('folder',Folder(iw.imagepath));
            end

            % Restrict files to indexes
            id.Folder.Files = id.Folder.Files([min(index) max(index)]);


        end
    end
end