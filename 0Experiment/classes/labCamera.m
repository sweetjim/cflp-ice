classdef labCamera<handle& matlab.mixin.SetGet
    %IWCAMERA connects to a camera

    properties
        operational logical = false
        camera 
        settings
        ds double = 1
        datastore matlab.io.datastore.ImageDatastore = imageDatastore(pwd)
    end
    properties
        label string = ''
        % default image data label          ['iw_label... .tiff']
        savepath string = '..\Data\ablation_images\'
        % image acquisition savepath
        interval duration = seconds(60)
        % image acquisition interval        [seconds]
        duration duration = seconds(10*60^2)
        % image acquisition duration        [seconds]
        elapsed duration = seconds(0)
        % image acquisition elapsed time    [datetime]
        count   double = 0
        % image acquisition image count
        starttime datetime
        % image acquisition start time      [datetime]
    end
    properties (Hidden)
        ImageTime struct = struct('recieve',struct('last',[],'av',nan),'send',struct('last',[],'av',nan))
        iw iwdata
        endtime datetime
        % image acquisition end time        [datetime]
        liveView
        camtimer timer
        camTimerState string = 'normal'
        fnamefmt function_handle = @(obj) sprintf( ...
            'cam%i_%sC%iT%s.tiff', ...
            obj.camera.DeviceID,...
            obj.label, ...
            obj.count, ...
            datetime('now','Format','uuuuMMdd''T''HHmmss'))
            %strrep(string(obj.elapsed,'hh:mm:ss'),':','_'))
        fnameformat string = 'cam<device_id|int>_<label|str>C<image_count|int>clk<datetime|uuuuMMddTHHmmss>'
                            % 'cam<device_id|int>_<label|str>C<image_count|int>T<elapsed_time|hh:mm:ss>'
        timeset struct = struct('Time',datetime,'From','Now','Value',seconds(0))
        UserData struct = struct( ...
            'Rot',0, ...
            'ROI',[1 1 3088 2064], ...
            'ds',nan, ...
            'command',@(camera_handle,userdata)rot90(imcrop(getsnapshot(camera_handle),userdata.ROI),userdata.Rot))
        UserDataPath string = '0Experiment\classes\camera_profiles\'
    end
    %% Constructor
    methods
        function c = labCamera(cam_id,savepath)
            % Establish uplink
            if nargin<1
                cam_id = 1;
            end
            
            if isa(cam_id,'videoinput')&&~isempty(cam_id)
                c.camera  = cam_id;
                cam_id    = cam_id.DeviceID;
            else
                try 
                    c.camera = videoinput("gentl", cam_id, "Mono8");
                catch
                    c.camera    = [];
                    cam_id      = nan;
                end
            end

            % Define camera resolutions
            switch cam_id
                case 1
                    c.ds = 0.0192*1e-2; % m/px
                case 2
                    c.ds = 3.9529e-04; % m/px
                otherwise
                    c.ds = []; % m/px
                    return
            end

            % look for camera userdata
            udpath = fullfile(c.UserDataPath,[c.camera.Name '_userdata.xml']);
            if isfile(udpath)
                ud = readstruct(udpath);
                try ud.command=str2func(ud.command);end
                try c.ds=ud.ds;end
                c.camera.UserData=ud;
                c.UserData=ud;
            end

            if nargin<3
                savepath = c.savepath;
            end
            c.savepath    = savepath;
            c.settings    = getselectedsource(c.camera);
        end
    end
    %% Settings
    methods
        function changeSettings(iwc)
            [~,~,imdata,hImage]=iwc.getLiveFeed(false,true);
            setappdata(hImage,'UpdatePreviewWindowFcn',@liveView); % livestream loop
            function liveView(~,event,~)
                % broadcast preview output to frontend ui
                iwc.liveView = event;
                try %#ok<*TRYNC>
                    imdata.CData = rot90(getsnapshot(iwc.camera),-1); % update CData (faster than calling imagesc)
                end
            end
        end
    end
    %% Acquisition
    methods
        function timelapse(iwc,ReadOutsOn)
            %% Set start time
            if nargin<2
                ReadOutsOn = true;
            end
            startnow = false;
            if iwc.starttime<datetime('now')
                startnow=true;
            end

            iwc.endtime     = iwc.starttime+seconds(iwc.duration);

            %% Read-outs
            if ReadOutsOn;iwc.status;end
            %% Open camera and timer
            % check if live-feeds are active
            openCamera(iwc,true)

            if isempty(iwc.savepath)
                iwc.setDirectory
            end
            warning off,delete(iwc.camtimer),warning on
            iwc.camtimer = timer(...
                "Period",seconds(iwc.interval),...
                'TasksToExecute',floor(iwc.duration/iwc.interval), ...
                'ExecutionMode','fixedRate',...
                'BusyMode','queue',...
                'StartFcn',@(~,~) windowsSound("notify"),...
                'TimerFcn',@acquisitionFcn,...
                'StopFcn',@stopTimeLapse,...
                'Tag','camtimer',...
                'UserData','normal');

            if startnow||iwc.starttime<datetime('now')
                start(iwc.camtimer)
            else
                startat(iwc.camtimer,iwc.starttime)
            end
            if ReadOutsOn;delayedMessage('Image acquisition started',iwc.starttime);end
            iwc.operational = false;
            
            %% Timer operation
            function acquisitionFcn(~,~)
               iwc.camTimerAcqFcn
            end
            function stopTimeLapse(~,~)
               iwc.camTimerStopFcn
               if ReadOutsOn
                   fprintf('\nImage acquisition stopped\n')
                   iwc.status
               end
            end
        end
        %% Timer functions
        function camTimerAcqFcn(iwc)
            if ~iwc.operational
                iwc.operational = true;
            end
            iwc.elapsed = datetime('now')-iwc.starttime;
            iwc.count   = iwc.count+1;
            iwc.writeImage(iwc.snapshot)
        end
        function camTimerStopFcn(iwc)
            switch iwc.camtimer.UserData
                case 'normal' % normal stop and delete routine
                    try windowsSound("notify"),end
                    iwc.operational = false;
                    iwc.closeCamera
                    delete(iwc.camtimer)
                case 'pause' % wait for start
                case 'restart'
                    iwc.camtimer.TasksToExecute = floor(iwc.duration/iwc.interval);
                    start(iwc.camtimer)
                    iwc.camtimer.UserData = 'normal';
            end
        end

        %% Acquisition functions
        function im = snapshot(iwc)
            t0          = tic;
            im = rot90(getsnapshot(iwc.camera),-1);
            iwc.ImageTime.recieve.last = toc(t0);
            iwc.ImageTime.recieve.av   = mean([iwc.ImageTime.recieve.av iwc.ImageTime.recieve.last],'omitnan');
        end
        function writeImage(iwc,im)
            t1          = tic;
            imwrite(im, ...
                fullfile( ...
                iwc.savepath, ...
                iwc.fnamefmt(iwc) ...
                ), ...
                'tiff', ...
                'Compression','none')
            iwc.ImageTime.send.last = toc(t1);
            iwc.ImageTime.send.av   = mean([iwc.ImageTime.send.av iwc.ImageTime.send.last],'omitnan');
        end
        function [hImage,camfig] = openCamera(iwc,checkIfOpen)
            if nargin<2
                checkIfOpen = false;
            end
            cam     = iwc.camera;
            vidRes  = cam.VideoResolution;
            nBands  = cam.NumberOfBands;
            tag     = sprintf('LiveView_%s_source',iwc.camera.Name);
            camfig  = findall(0,'Type','Figure','Tag',tag);

            if isempty(camfig)
                camfig  = figure('Name', 'LiveView: Preview', ...
                    'Visible','off', ...
                    'HandleVisibility','callback',...
                    'Tag',tag);
                axes(camfig)
            elseif checkIfOpen
                % ensure preview is active
                %...
                return
            end

            if checkIfOpen
                fprintf('Opening %s\n',iwc.camera.Name)
            end
            hImage  = image(zeros(vidRes(2), vidRes(1), nBands), ...
                'Parent',findall(camfig,'Type','Axes'));
            preview(cam, hImage);
            set(camfig,'Visible','off')
            
            if ~nargout
                clear hImage camfig
            end

        end
        function closeCamera(~)
            closepreview
            close(findall(0,'Type','Figure','Tag','LiveView'))
        end
        function setDirectory(iwc,path)
            if nargin<2
                path = iwc.savepath;
            end
            path = uigetdir(path,sprintf('Set image repository (%s)',iwc.camera.Name));
            if isa(path,'numeric')
                warning('Operation cancelled.')
                return
            end
            iwc.savepath = path;
        end
        function testConditions(iwc,starttime,noClear)

            if nargin<2
                starttime = 0;
            end
            if nargin<3
                noClear = false;
            end

            % defaults
            savepath0       = iwc.savepath;
            label0          = iwc.label;
            duration0       = iwc.duration;
            count0          = iwc.count;
            interval0       = iwc.interval;


            iwc.label       = 'test';
            iwc.savepath    = 'iwCam_test\';

            if ~isfolder(iwc.savepath)
                mkdir(iwc.savepath)
            end

            if ~noClear
                files = dir(fullfile(iwc.savepath,'*.tiff'));
                arrayfun(@(x) delete(fullfile(files(x).folder,files(x).name)),1:numel(files))
            end
            fprintf('TEST CONDITIONS\n')

            iwc.count       = 0;
            iwc.duration    = 10;
            iwc.interval    = 2;

            winopen(iwc.savepath)
            iwc.timelapse(starttime)
            delayedMessage('TEST CONCLUDED',iwc.duration)

            % revert changes
            delayedFunction(@revertChanges,iwc.duration)
            function revertChanges
                iwc.savepath    = savepath0;
                iwc.label       = label0;
                iwc.duration    = duration0;
                iwc.count       = count0;
                iwc.interval    = interval0;
            end
        end
    end
    %% Datastore
    methods
        function state = readDataStore(iwc)
            state = true;
            try
                iwc.datastore = imageDatastore(iwc.savepath);
            catch
                state = false;
            end
        end
        function openPath(iwc)
            winopen(iwc.savepath)
        end
    end
    %% Read-outs
    methods
        function status(iwc)
            %% Reports status of camera and image logging
            
            ImagesInPath = 0;
            if iwc.readDataStore
                ImagesInPath=numel(iwc.datastore.Files);
            end
            clear o s
            s{1}=sprintf('Camera "%s"',iwc.camera.Name);o{1}='';
%             if iwc.operational
%                 o{1} = sprintf(' (Running)');
%             else
%                 o{1} = sprintf(' (Not running)');
%             end
            s{end+1}='';o{end+1}='';
            if ~isempty(iwc.starttime)
                s{end+1}='Start-up:';
                o{end+1}=sprintf(' %s',string(iwc.starttime,'dd-MMM-uuuu hh:mm:ss'));
            end
            if ~isempty(iwc.endtime)
                s{end+1}='Shut-down:';
                o{end+1}=sprintf(' %s',string(iwc.endtime,'dd-MMM-uuuu hh:mm:ss'));
            end
            s{end+1}='Duration:';
            s{end+1}='Images:';
            s{end+1}='Interval:';
            s{end+1}= 'Average acquisition time:';
            s{end+1}= 'Average writing time:';
            o{end+1}=sprintf(' %s (HH:MM:SS)',string(iwc.duration,'hh:mm:ss'));
            o{end+1}=sprintf(' %i (taken) | %i (in path)',iwc.count,ImagesInPath);
            o{end+1}=sprintf(' %is',seconds(iwc.interval));
            o{end+1}=sprintf(' %.2fms',iwc.ImageTime.recieve.av*1e3);
            o{end+1}=sprintf(' %.2fms',iwc.ImageTime.send.av*1e3);
            try timerON = isvalid(iwc.camtimer);catch,timerON=true;end
            if timerON
                try
                    if ~timerEMPTY
                        s{end+1}='Execution Mode:';
                        o{end+1}=sprintf('%s',iwc.camtimer.ExecutionMode);
                    end
                end
            end
            s=[pad(s','left') pad(o','right')];
            s=arrayfun(@(x) cat(2,s{x,:}),1:size(s,1),'UniformOutput',false)';
            arrayfun(@(x) fprintf('%s\n',x{1}),s)
        end
    end
    %% Calibration
    methods
        function calibrate(iwc)

        end
    end
    %% Live preview
    methods
        function seeLiveView(iwc)
            [~,~,imdata,hImage] = iwc.getLiveFeed;
            setappdata(hImage,'UpdatePreviewWindowFcn',@liveView); % livestream loop
            function liveView(~,event,~)
                % broadcast preview output to frontend ui
                iwc.liveView = event;
                try %#ok<*TRYNC>
                    imdata.CData = rot90(getsnapshot(iwc.camera),-1); % update CData (faster than calling imagesc)
                end
            end
        end
        function seeSytheticSchilieren(iwc,closeLiveView)
            %% Preamble
            % close liveview command
            if nargin<2
                closeLiveView = false;
            end

            iwc.iw = iwdata;
            [f,ax,imdata,hImage,im0] = iwc.getLiveFeed(closeLiveView);
            setappdata(hImage,'UpdatePreviewWindowFcn',@liveView); % livestream loop

            %% UI controls
            poscon = [20 20 60 20];
            shift = @(x,y) poscon+[x 0 y 0];
            uicontrol(f(2), ...
                'String','New reference image', ...
                'Position',shift(0,90),...
                'Tooltip','Resets the reference image',...
                'Callback',@newRefIm);
            field = uicontrol(f(2),...
                'Style','popupmenu',...
                'Tooltip','Output velocity field',...
                'Position',shift(200,0),...
                'String',{'U','u','v'});
            resizer = uicontrol(f(2),...
                'Style','edit','String','.25',...
                'Position',shift(280,0),...
                'Tooltip','Resizing factor');
            avU = uicontrol(f(2),...
                'Style','edit','String','5',...
                'Position',shift(360,50),...
                'Tooltip','Averaging window');

            uicontrol('Style','text','Position',field.Position+[0 20 0 10],'Parent',f(2),'String','Velocity field');
            uicontrol('Style','text','Position',resizer.Position+[0 20 0 10],'Parent',f(2),'String','Resize (x100 %)');
            uicontrol('Style','text','Position',avU.Position+[0 20 0 10],'Parent',f(2),'String','Averaging window (frames)');

            % roi and mask
            [Y,X] = size(im0);
            roi     = [0 0 X Y]*iwc.ds*1e3;%[184,153,678,2686];%[1267,215,503,2368];
            pgon    = drawrectangle('Position',roi,'FaceAlpha',0,'Parent',ax,'Label','ROI');
            addlistener(pgon,'MovingROI',@roiCallback);

            % open output figure
            if isempty(findall(0,'Tag','playground'))
                figure('Tag','playground')
            else
                figs = findall(0,'Tag','playground');
                set(0,'CurrentFigure',figs(1))
            end
            clf
            tiledlayout(1,1,'Padding','compact')
            nexttile
            imagesc(im0,'Tag','output');
            cbar = addColorbar('location','southoutside','title','U');
            colormap gray
            caxis([0 1])
            set(gca,'XTick',[],'YTick',[])
            axis image

            %% Set broadcast variables
            U = []; % moving mean velocities
            V = [];
            cmap0 = 'gray';
            im0last = im0;

            %% Nested functions
            function liveView(~,event,~)
                % broadcast preview output to frontend ui
                iwc.liveView = event;
                try %#ok<*TRYNC>
                    im = rot90(getsnapshot(iwc.camera),-1);
                    imdata.CData = processimage(im); % update CData (faster than calling imagesc)
                end
            end
            function out = processimage(in)
                % Perform calculations on input image
                out = in; % ensure output is defined in case of catch errors
                id  = findall(0,'Tag','output','Type','Image');  % get output imagedata
                try delete(id(2:end));end        % remove additional Image types
                par = id.Parent;                 % get parent of imagedata (axes)

                tic

                % masking
                %                 mask = createMask(maskgon);
                imref = double(im0);
                imdef = double(in);

                % cropping and resizing
                imref = imcrop(imref,roi./(iwc.ds*1e3));
                imdef = imcrop(imdef,roi./(iwc.ds*1e3));

                %                 clc
                %                 resizer.String
                factor = str2double(resizer.String);
                imref = imresize(imref,factor);
                imdef = imresize(imdef,factor);

                par.Title.String = '';
                par.Title.Color = 'k';
                inputHandler(imdef)

                %% Velocity computation
                try
                    warning off
                    [u,v]           = iwc.iw.ss_solver_fcd(imref,imdef);
                    window = str2double(avU.String);
                    getMeanVelocity(u,v,window+1);

                    if window>0
                        umean = mean(U,3);
                        vmean = mean(V,3);
                        u = umean-u;
                        v = vmean-v;
                    end
                    switch field.String{field.Value}
                        case 'u'
                            out  = u;
                            cmap = 'balance';
                            str = "Horizontal velocity (normalized)";%"u'=\langleu\rangle-u";
                        case 'v'
                            out  = v;
                            cmap = 'balance';
                            str = "Vertical velocity (normalized)";%"v'=\langlev\rangle-v";
                        case 'U'
                            out = hypot(u,v);
                            cmap = 'gray';
                            str = "Speed (normalized)";%"U'=\langleU\rangle-U";
                    end
                    if isvalid(cbar)
                        cbar.Label.String = str;
                    end
                    inputHandler(out)
                    if ~strcmp(cmap0,cmap)||all(caxis(par)==[0 255])
                        cmap0 = cmap;
                        set(par,'Colormap',1-cmapCheck(cmap))
                        switch cmap
                            case 'gray'
                                caxis(par,[0 1])
                            case 'balance'
                                caxis(par,[-1 1])
                        end
                    end
                    tcomp = toc;
                    subtitle(par,sprintf('%.2f FPS %s%i%s',1/tcomp,'@',factor*1e2,'%'))
                    hold(par,'on')
                    delete(findall(par,'Type','Quiver'))
                    thresh = .1;%std2(hypot(u,v))*.5e2;
                    u(abs(u)<thresh/10)=0;
                    v(abs(v)<thresh/10)=0;
                    %vis_flow(u,v,'ax',par,'gx',50,'col','r','mag',thresh)
                    %pause(tcomp)
                    warning on
                catch ME
                    caxis(par,[0 255])
                    delete(findall(par,'Type','Quiver'))
                    set(par,'Colormap',cmapCheck('gray'))
                    par.Title.String = 'Velocity calculation error';
                    par.Title.Color = 'r';
                end
                out=in;

                %% Nested function
                function inputHandler(in)
                    try
                        id.CData=in;
                        XLIM = [0 size(id.CData,2)];
                        YLIM = [0 size(id.CData,1)];
                        if any(xlim(par)~=XLIM)
                            xlim(par,XLIM)
                        end
                        if any(ylim(par)~=YLIM)
                            ylim(par,YLIM)
                        end
                    catch
                        imagesc(par,in,'Tag','output')
                    end
                end
            end
            function newRefIm(~,~)
                % Callback for
                im0 = rot90(getsnapshot(iwc.camera),-1);
                t = timer('StartDelay',0, ...
                    'TimerFcn',@getMeanstate, ...
                    'Period',1, ...
                    'TasksToExecute',size(U,3),...
                    'StopFcn',@stoptimer,...
                    'ExecutionMode','fixedRate');
                k=1;
                ims = zeros([size(im0),t.Period*t.TasksToExecute]);
                start(t)
                Fdata = findall(0,'Tag','output');
                %                 par = Fdata.Parent;
                %                 par.Title.String = sprintf('Gathering new reference state (%is)',t.TasksToExecute-k-1);

                function getMeanstate(~,~)
                    %                     par.Title.String = sprintf('Gathering new reference state (%is)',t.TasksToExecute-k-1);
                    ims(:,:,k) =  rot90(getsnapshot(iwc.camera),-1);
                    k = k+1;
                end
                function stoptimer(timeriwc,~)
                    im0 = mean(ims,3);
                    %                     par.Title.String = '';
                    delete(timeriwc)
                end
            end
            function getMeanVelocity(u,v,steps)
                if nargin<3
                    steps = 20;
                end
                if isempty(U)
                    U = zeros([size(u) steps]);
                    V = U;
                end
                if any(size(U,[1 2])~=size(u))||size(U,3)~=steps
                    U = zeros([size(u) steps]);
                    V = U;
                end
                U(:,:,1) = u;
                V(:,:,1) = v;
                U = circshift(U,1,3);
                V = circshift(V,1,3);
            end
            function roiCallback(source,~)
                roi = source.Position;
            end
        end
        function [f,ax,imdata,hImage,im0] = getLiveFeed(iwc,closeLiveView,showControls)
            %% Preamble
            % Check if camera is connected
            if ~isvalid(iwc.camera)
                try
                    iwc.camera  = videoinput("gentl", 1, "Mono8");
                catch
                    error('Unable to connect to camera')
                end
            end
            % close liveview command
            if nargin<2
                closeLiveView = false;
            end
            if nargin<3
                showControls = false;
            end

            % check if liveview is already open (and restart if so)
            tag = sprintf('LiveView_%s',iwc.camera.Name);
            openfig = findall(0,'Type','Figure','Tag',tag);
            pos     = [];
            if (~isempty(openfig)||closeLiveView)&&~showControls
%                 f2pos = findall(0, ...
%                     'Type','Figure', ...
%                     'Tag',strcat(tag,'_source'));
%                 pos = f2pos.Position;
                if closeLiveView
                    close(openfig)
                    return
                end
                set(0,"CurrentFigure",openfig),shg
                imdata      = findall(openfig,'Type','Image');
                ax          = findall(openfig,'Type','Axes');
                [hImage,f]  = iwc.openCamera;
                return
            end

            % open live feed (preview) in background
            [hImage,f] = iwc.openCamera;
            uicontrol(f,'String', 'Close', 'Callback', 'close(gcf)');
            set(gcf,'Visible','off')

            % open roi/control figure
            figlabel = sprintf('LiveView: %s',iwc.camera.Name);
            if ~showControls
                f(2)=figure('Name',figlabel, ...
                    'Tag',tag, ...
                    'IntegerHandle','off',...
                    'HandleVisibility','callback',...
                    'CloseRequestFcn',@closeLiveViewFcn);
                ax = axes;
            else
                f(2)=uifigure('Name',figlabel, ...
                    'Tag',tag, ...
                    'IntegerHandle','off',...
                    'HandleVisibility','callback',...
                    'CloseRequestFcn',@closeLiveViewFcn);
                ax = buildControls(f(2));
            end

            if ~isempty(pos)
                f(2).Position = pos;
            end
            ax.Toolbar.Visible='off';
            im0    = rot90(getsnapshot(iwc.camera),-1);
            x      = 0:iwc.ds:size(im0,2)*iwc.ds;
            y      = 0:iwc.ds:size(im0,1)*iwc.ds;
            imdata = imagesc(ax,x*1e3,y*1e3,im0,'Tag','imdata');
            addlabels(ax,'x','Position (mm)','y','Depth (mm)','fs',10)
            caxis(ax,[0 255])
            colormap(ax,'gray')
%             set(ax,'XTick',[],'YTick',[])
            axis(ax,'image')

            %% Nested functions
            function closeLiveViewFcn(~,~)
                closepreview
                close(f(1),'force')
                delete(f(2))
            end
            function ax = buildControls(f)
                %%
                s   = iwc.settings;
                gl  = uigridlayout(f,'ColumnWidth',{200,'1x'},'RowHeight',{'1x'});
                p   = uipanel(gl,'Title','Settings','BorderType','none');
                gl2 = uigridlayout(p,'ColumnWidth',{100,'1x'},'RowHeight',repmat({30},1,12));
                glax = uigridlayout(gl,'ColumnWidth',{'1x'},'RowHeight',{'1x',50});
                ax  = uiaxes(glax);
                lbl = uilabel(glax,'Text',sprintf('%.1f FPS',1/(iwc.settings.ExposureTime*1e-6)));


                a1=createUIElement(sprintf('Exposure (%ss)',char(956)),'ExposureTime',s.ExposureTime);
                createUIElement('Gain','Gain',s.Gain);
                createUIElement('Gamma','Gamma',s.Gamma);
                a2=createUIElement('Digital Shift','DigitalShift',double(s.DigitalShift),'edit',[0 4]);
                set([a1 a2],ValueDisplayFormat='%i',RoundFractionalValues='on')
                createUIElement('Reverse X','ReverseX',strcmpi(s.ReverseX,'true'),'button');
                createUIElement('Reverse Y','ReverseY',strcmpi(s.ReverseY,'true'),'button');
                createUIElement([],[],[],'roi');
                createUIElement([],[],[],'resetroi');
                uilabel(gl2,'Text','');
                uilabel(gl2,'Text','');
                createUIElement('Set ds (m/px)','ds',0,'pushbutton',[],'on',@dsCallback);
                uilabel(gl2,'Text','');
                createUIElement('Distance (cm)','x',0,'edit',[0 inf],'off',@dsCallback);
                createUIElement('Accept','accept',[],'pushbutton',[],'off',@dsCallback);
                createUIElement('Cancel','cancel',[],'pushbutton',[],'off',@dsCallback);

                %% Nested function
                function el = createUIElement(label,tag,value,type,limits,visibility,functionCallback)
                    if nargin<4
                        type = 'edit';
                    end
                    if nargin<5
                        limits = [0 inf];
                    end
                    if nargin<7
                        functionCallback = @changeParameter;
                    end
                    if nargin<6
                        visibility='on';
                    end
                    switch type
                        case 'edit'
                            uilabel(gl2,'Text',label,...
                                'Tag',tag,...
                                'Visible',visibility);
                            el = uieditfield(gl2,'numeric', ...
                                'Limits',limits, ...
                                'Tag',tag,...
                                'Value',value, ...'
                                'ValueChangedFcn',functionCallback, ...
                                'ValueDisplayFormat','%.1f', ...
                                'RoundFractionalValues','off',...
                                'Visible',visibility);
                        case 'button'
                            el = uibutton(gl2,'state',...
                                'Text',label, ...
                                'Value',value,...
                                'Tag',tag,...
                                'ValueChangedFcn',functionCallback,...
                                'Visible',visibility);
                        case 'roi'
                            el = uibutton(gl2,'push',...
                                'Text','Set ROI',...
                                'Tag','setroi',...
                                'ButtonPushedFcn',functionCallback,...
                                'Visible',visibility);
                        case 'resetroi'
                            el = uibutton(gl2,'push',...
                                'Text','Reset ROI',...
                                'Tag','resetroi',...
                                'ButtonPushedFcn',functionCallback,...
                                'Visible',visibility);
                        case 'pushbutton'
                            el = uibutton(gl2,'push',...
                                'Text',label,...
                                'Tag',tag,...
                                'ButtonPushedFcn',functionCallback,...
                                'Visible',visibility);
                    end
                end
                function changeParameter(src,~)
                    imdata = findall(ax,'Type','Image','Tag','imdata');
                    warning off
                    switch src.Tag
                        case 'setroi'
                            roi = drawrectangle(ax,'FaceAlpha',0,'Label','ROI');
                        case 'resetroi'
                            set(iwc.camera,ROIPosition=[0 0 get(iwc.camera,'VideoResolution')])
                        case 'DigitalShift'
                            iwc.settings.DigitalShift = uint32(src.Value);
                        case {'ReverseX','ReverseY'}
                            val = 'False';
                            if src.Value
                                val = 'True';
                            end
                            command = sprintf('iwc.settings.%s="%s";',src.Tag,val);
                            eval(command)
                        otherwise
                            command = sprintf('iwc.settings.%s=%g;',src.Tag,src.Value);
                            try eval(command)
                            catch
                                iwc.settings.GainAuto = 'off';
                                iwc.settings.ExposureAuto = 'off';
                                eval(command)
                            end
                    end
                    lbl.Text = sprintf('%.1f FPS',1/(iwc.settings.ExposureTime*1e-6));
                    warning on
                end
                function dsCallback(src,~)
                    ref     = size(imdata.CData);
                    set(ax, ...
                        'XLim',[0 ref(2)], ...
                        'YLim',[0 ref(1)])
                    ylabel(ax,'Depth (px)')
                    xlabel(ax,'Position (px)')
                    set(imdata, ...
                        'XData',0:ref(2)-1, ...
                        'YData',0:ref(1)-1)
                    el      = findall(gl2,'-not','Type','uigridlayout');
                    disable = findall(el, ...
                        '-not','Tag','accept', ...
                        '-not','Tag','cancel',...
                        '-not','Tag','x');
                    keep    = findall(el, ...
                        'Tag','accept','-or', ...
                        'Tag','cancel','-or',...
                        'Tag','x');
                    %cellfun(@isempty,arrayfun(@(x) cell2mat(regexp(x.Tag,{'accept','cancel','x'})),el,'UniformOutput',false));
                    switch src.Tag
                        case 'x'
                            return
                        case 'ds'
                            set(disable,'Enable','off')
                            set(keep,'Visible','on')
                            ax.Toolbar.Visible='on';
                            imdistline(ax);
                            return
                        case 'accept'
                            Xpx=findall(ax,'Tag','imline');
                            Xcm=findall(el,'Tag','x','Type','uinumericeditfield').Value;
                            pospx = diff(...
                                cell2mat(arrayfun(@(x) [Xpx.Children(x).XData Xpx.Children(x).YData],[2 3],'UniformOutput',false)'));
                            Xpx=hypot(pospx(1),pospx(2));
                            iwc.ds = Xcm/Xpx*1e-2;
                            
                            ref     = size(imdata.CData);
                            maxDims = ref*iwc.ds;
                            set(imdata, ...
                                'XData',linspace(0,maxDims(2),ref(2)), ...
                                'YData',linspace(0,maxDims(1),ref(1)))
                            set(ax, ...
                                'XLim',[0 inf], ...
                                'YLim',[0 inf])
                            ylabel(ax,'Depth (m)')
                            xlabel(ax,'Position (m)')
                    end

                    delete(findall(ax,'Tag','imline'))
                    set(disable,'Enable','on')
                    set(keep,'Visible','off')
                    ax.Toolbar.Visible='on';
                end
            end
        end
    end
    %% Set/Get/Reset
    methods
        function setCamTimerState(iwc)
            try
                iwc.camtimer.UserData = iwc.camTimerState;
            catch ME
                if ~strcmp(ME.message,'Invalid or deleted object.')
                    rethrow(ME)
                end
            end
        end
        function reset(iwc)
            warning('Reseting count and elapsed time')
            iwc.count=0;
            iwc.elapsed=[];
            iwc.ImageTime = struct('recieve',struct('last',[],'av',nan),'send',struct('last',[],'av',nan));
        end
        function set.camTimerState(iwc,val)
            iwc.camTimerState = val;
            setCamTimerState(iwc)
        end
        function set.count(iwc,val)
            iwc.count=val;
        end
        function set.interval(iwc,val)
            iwc.interval=val;
        end
        function set.duration(iwc,val)
            iwc.duration=val;
        end
        function set.savepath(iwc,val)
            iwc.savepath = val;
        end
        function set.starttime(iwc,val)
            iwc.starttime=val;
        end
        function set.label(iwc,val)
            iwc.label = val;
        end
    end
end

