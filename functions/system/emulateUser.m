classdef emulateUser<handle
    %EMULATEUSER contains a list of functions that emulate user input or
    %utilise user input.

    properties
        data
        tmr
        state
    end
    properties (Hidden)
        formerWindowsButtonFcns
        rbt
    end
    %% Trace mouse pointer velocity
    methods
        function usr = emulateUser
            usr.rbt = java.awt.Robot;
        end
        function trackCursorOnAxes(usr,ax,onClick,period,timeout,ext_fcn)
            if nargin<2
                ax = gca;
            end
            if nargin<3
                onClick = false;
            end
            if nargin<4
                period = 0.01;
            end
            if nargin<5
                timeout = 5;
            end
            if nargin<6
                ext_fcn = [];
            end
            f = ancestor(ax,'Figure');
            idx = 1;
            t2e = timeout/period;
            log = nan([t2e 2]);

            [~,~,fBM,~,~] = checkWindowsButtonFcns(usr,f);


            if onClick
                t2e = 1e5;
                log = nan([1e5 2]);
                f.KeyPressFcn           = @keyActivate;
                f.KeyReleaseFcn         = @keyActivate;
                f.WindowButtonDownFcn   = @clickTimerOn;
                f.WindowButtonUpFcn     = @clickTimerOff;
                usr.data = {};
                clicks = 1;
            end
            ints = ax.Interactions;
            ax.Interactions = [];
            f.WindowButtonMotionFcn = @trackMouse;
            xbounds = ax.XLim;
            ybounds = ax.YLim;
            
            [dP,dS,dV] = initialiseTracePlot;
            t = buildTimer(period);
            usr.tmr = t;

            currentpoint    = [];
            usr.state       = false;
            afterClick      = false;
            keep            = [];
            if ~onClick
                start(t);
            end
            function parseFcn(ext_fcn)
                if isempty(ext_fcn)
                    return
                end
                if usr.state
                    ext_fcn();
                end
            end
            function clickTimerOn(~,~)
                if usr.state
                    t.StopFcn = @StopFcnClick;
                    xbounds = ax.XLim;
                    ybounds = ax.YLim;
                    try
                        start(t)
                    catch
                        t = buildTimer(period);
                        t.StopFcn = @StopFcnClick;
                        start(t)
                    end
                end
            end
            function clickTimerOff(~,~)
                try %#ok<*TRYNC>
                    if strcmp(t.Running,'off')
                        return
                    end
                    stop(t)
                    newpts = log(1:idx,:);
                    newpts = smooth2a(newpts,2,0);
                    usr.data{clicks} = newpts;
                    afterClick = true;
                    log = log*nan;
                    idx = 1;
                    clicks = clicks+1;
                    showPlot
                    afterClick = false;
                end
                parseFcn(ext_fcn)

            end
            function keyActivate(~,event)
                if ~strcmp(event.Key,'shift')
                    return
                end

                switch event.EventName
                    case 'KeyPress'
                        if isempty(keep)
                            keep = findall(ax,'Type','Images','Visible','on');
                            keep = cat(1,keep,findall(ax,'Type','Line','Visible','on'));
                        end
                        usr.state = true;
                        set(f,'Pointer','crosshair')
                        %                         if isempty(ext_fcn)
                        title(ax,'ready')
                        set(keep(arrayfun(@isvalid,keep)),'Visible','off')
                        %                         end
                    case 'KeyRelease'
                        usr.state = false;
                        set(f,'Pointer','arrow')
                        %                         if isempty(ext_fcn)
                        title(ax,'waiting')
                        set(keep(arrayfun(@isvalid,keep)),'Visible','on')
                        keep = [];
                        %                         end
                end
            end
            function trackMouse(~,~)
                currentpoint = get(ax,'CurrentPoint');
            end
            function StopFcnClick(~,~)
                %disp('Waiting for new click')
            end
            function StopFcnNormal(~,~)
                f.WindowButtonMotionFcn = fBM;
                disp('Timer finished')
                ax.Interactions = ints;
            end
            function ErrorFcn(~,~)
                ax.Interactions = ints;
                warning('Error occured')
                stop(t)
            end
            function timerFcn(~,~)
                if strcmp(get(f,'Pointer'),'arrow')
                    set(f,'Pointer','crosshair')
                end
                try
                    xval = currentpoint(1);
                    yval = currentpoint(3);
                catch
                    xval = nan;
                    yval = nan;
                end
                if ~within(xval,xbounds)
                    xval = nan;
                end
                if ~within(yval,ybounds)
                    yval = nan;
                end
                log(idx,:)=[xval yval];
                showPlot
                idx=idx+1;
            end
            function showPlot
                if afterClick
                    try
                        hold(ax,'on')
                        val = cell2mat(usr.data');
                        X = val(:,1);
                        Y = val(:,2);
                        plotToAxes(X,Y)
                        %                         cellfun(@(x) plotToAxes(x(:,1),x(:,2),usr.data))
                        hold(ax,'off')
                    end
                    return
                end
                plotToAxes(log(:,1),log(:,2))
                function plotToAxes(X,Y)
                    v = gradient([X Y]')';
                    v = log10(sqrt(v(:,1).^2+v(:,2).^2));
                    v(isinf(v))=0;

                    if ~isvalid(dP)||~isvalid(dS)||~isvalid(dV)
                        [dP,dS,dV]=initialiseTracePlot;
                    end
                    dP.XData = X;
                    dP.YData = Y;
                    dS.XData = X;
                    dS.YData = Y;
                    dS.CData = v;
                    dV.XData = [X X];
                    dV.YData = [Y Y];
                    dV.CData = [v v];
                    warning off
                    dV.ZData = [v v]*0;
                    warning on
                end
            end
            function t = buildTimer(period)
                t = timer('ExecutionMode', 'fixedRate', ...
                    'Period', period, ...
                    'TasksToExecute', t2e, ...
                    'TimerFcn',@timerFcn,...
                    'StopFcn',@StopFcnNormal,...
                    'ErrorFcn',@ErrorFcn,...
                    'Tag','trace');
            end
            function [dP,dS,dV]=initialiseTracePlot
                hold(ax,'on')
                dP = plot(ax,log(:,1),log(:,2),'Color','k','Tag','trace');
                dS = scatter(ax,log(:,1),log(:,2),10,zeros(1,size(log,1)),'filled','Tag','trace');
                dV = coloredLine(log(:,1),log(:,2),zeros(1,size(log,1)),ax);
                dV.Tag = 'trace';
            end

        end
        function [fBD,fBU,fBM,fKP,fKR,KP,KR] = checkWindowsButtonFcns(usr,f)
            fBD=f.WindowButtonDownFcn;
            fBU=f.WindowButtonUpFcn;
            fBM=f.WindowButtonMotionFcn;
            fKP=f.WindowKeyPressFcn;
            fKR=f.WindowKeyReleaseFcn;
            KP=f.KeyPressFcn;
            KR=f.KeyReleaseFcn;
            usr.formerWindowsButtonFcns = {f,fBD,fBU,fBM,fKP,fKR,KP,KR};
        end
        function resetWindowsButtonFcns(usr)
            fcns = usr.formerWindowsButtonFcns;
            f       = fcns{1};
            fcns    = fcns(2:end);
            f.WindowButtonDownFcn   = fcns{1};
            f.WindowButtonUpFcn     = fcns{2};
            f.WindowButtonMotionFcn = fcns{3};
            f.WindowKeyPressFcn     = fcns{4};
            f.WindowKeyReleaseFcn   = fcns{5};
            f.KeyPressFcn           = fcns{6};
            f.KeyReleaseFcn         = fcns{7};
        end
        function clearData(usr)
            usr.data = [];
        end
    end
    %% (old trace methods)
    methods (Hidden)
        function watchUser(~,timeout,period)
            if nargin<2
                timeout = 5;
            end
            if nargin<3
                period = 0.01;
            end
            screensize=get(0,'ScreenSize');
            log = nan([timeout/period 2]);
            idx = 1;
            t = timer('ExecutionMode', 'fixedRate', ...
                'Period', period, ...
                'TasksToExecute', timeout/period, ...
                'TimerFcn',@mouseLocation,...
                'StopFcn',@stopfcn);
            start(t);
            clf
            plot(nan,nan)
            hold(gca,'on')
            function stopfcn(~,~)
                disp('Emulator timeout')
                hold(gca,'off')
                showPlot
                delete(t)
            end
            function mouseLocation(~,~)
                showPlot
                loc = get(0, 'PointerLocation');
                log(idx,:)=loc;
                idx=idx+1;
            end
            function showPlot
                plot(log(:,1),log(:,2))
                v = gradient(log')';
                v = log10(sqrt(v(:,1).^2+v(:,2).^2));
                v(isinf(v))=0;
                %%
                hold on
                plot(log(:,1),log(:,2),'Color',[0 0 0])
                scatter(log(:,1),log(:,2),10,...
                    rescale(v,0,1),'filled')
                hold off
                xlim([0 screensize(3)])
                ylim([0 screensize(4)])
            end
            function screenGrab
                % Take screen capture
                robot = usr.rbt;
                pos = screensize; % [left top width height]
                rect = java.awt.Rectangle(pos(1),pos(2),pos(3),pos(4));
                cap = robot.createScreenCapture(rect);
                % Convert to an RGB image
                rgb = typecast(cap.getRGB(0,0,cap.getWidth,cap.getHeight,[],0,cap.getWidth),'uint8');
                imgData = zeros(cap.getHeight,cap.getWidth,3,'uint8');
                imgData(:,:,1) = reshape(rgb(3:4:end),cap.getWidth,[])';
                imgData(:,:,2) = reshape(rgb(2:4:end),cap.getWidth,[])';
                imgData(:,:,3) = reshape(rgb(1:4:end),cap.getWidth,[])';
                % Show or save to file
                im = imgData;
                delete(robot)
            end
        end
    end
end

