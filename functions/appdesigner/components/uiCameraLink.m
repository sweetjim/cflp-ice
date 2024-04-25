classdef uiCameraLink < matlab.ui.componentcontainer.ComponentContainer

    % Properties that correspond to underlying components
    properties (Access = private, Transient, NonCopyable)
        GridLayout            matlab.ui.container.GridLayout
        NotificationsLabel    matlab.ui.control.Label
        GridLayout_3          matlab.ui.container.GridLayout
        PreviewButton         matlab.ui.control.Button
        ConnectButton         matlab.ui.control.Button
        DeviceIDListBox       matlab.ui.control.ListBox
        DeviceListBoxLabel    matlab.ui.control.Label
        GridLayout_2          matlab.ui.container.GridLayout
        AdaptorListBox        matlab.ui.control.ListBox
        AdaptorListBoxLabel   matlab.ui.control.Label
        CameraSettingsPanel   matlab.ui.container.Panel
        GridLayout2           matlab.ui.container.GridLayout
        timeSheduleLabel      matlab.ui.control.Label
        ScheduleDate          matlab.ui.control.DatePicker
        ScheduleDropDown      matlab.ui.control.DropDown
        BreaktimeLabel        matlab.ui.control.Label
        ScheduleLabel         matlab.ui.control.Label
        ImageintervalLabel    matlab.ui.control.Label
        ChangesettingsButton  matlab.ui.control.Button
        IntervalTime          uiDuration
        ScheduleTime          uiDuration
        IterationsTime        uiDuration
    end

    % Events with associated public callbacks
    events (HasCallbackProperty, NotifyAccess = private)
        SelectionChanged
    end

    properties (Access = private)
        adaptor
        device
        sheduleTimer timer = timer;
    end
    
    properties (Access = public)
        BreakTime duration = minutes(10);
        Interval duration = seconds(0);
        Iterations double = 0;
        SheduledTime datetime = datetime('now');
        camera
        labCamera
        test
        time {mustBeA(time, "duration")} = seconds(0);
    end
    
    
    methods (Access = public)
        
        function connect(comp)
            try
                comp.NotificationsLabel.Text = 'Attempting to connect to camera...';
                comp.camera = videoinput(comp.adaptor,comp.device);
                comp.NotificationsLabel.Text = 'Camera is active';
            catch
                comp.NotificationsLabel.Text = 'Failed to connect to camera';
            end
            
        end
    end
    
    methods (Access = private)
        function validCam = validcamera(comp)
            ids         = str2double(comp.DeviceIDListBox.Value);
            validids    = arrayfun(@(lc) lc.camera.DeviceID,comp.labCamera);
            validCam    = comp.labCamera(validids==ids);
        end
        function str = atimeReadOut(comp)
            fmt     = 'dd-MMM-uuuu HH:mm:ss';
            settime = comp.SheduledTime;
            endtime = settime;
            str = sprintf( ...
                'Set for %s\nEnds at %s\nElapsed time: %s\n%Image count: %i', ...
                string(settime,fmt), ...
                string(endtime,fmt), ...
                string(elapsedtime), ...
                imagecount);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startUpFcn(comp)
            warning off
            info = imaqhwinfo;
            warning on
            comp.AdaptorListBox.Items    = info.InstalledAdaptors;
            comp.NotificationsLabel.Text = '';
            comp.timeSheduleLabel.Text = '';
            comp.ConnectButton.UserData = false;

            comp.sheduleTimer = timer(Period = 1, ...
                ExecutionMode='fixedRate',...
                TasksToExecute=inf,...
                TimerFcn = @updateSheduledTimeDisplay, ...
                Tag = 'uiCameraLinkTimer');

            function updateSheduledTimeDisplay(~,~)
                comp.SheduledTime = string(datetime('now')+comp.ScheduleTime.SumTime);
                comp.timeSheduleLabel.Text = timeReadOut(comp);
            end
        end

        % Value changed function: AdaptorListBox, DeviceIDListBox
        function changeAttribute(comp, event)
            if comp.AdaptorListBox==event.Source
                if ~any(strcmp(comp.AdaptorListBox.Value,event.PreviousValue))||isempty(comp.DeviceIDListBox.Items)
                    % update device id list
                    warning off
                    info = imaqhwinfo(event.Value);
                    comp.DeviceIDListBox.Items = {};
                    warning on
                    if isempty(info)
                        return
                    end
                    if contains(lastwarn,'No devices were detected for the')&&isempty(info.DeviceIDs)
                        comp.NotificationsLabel.Text = 'No devices were detected';
                        return
                    end
                    comp.DeviceIDListBox.Items = arrayfun(@(x) sprintf('%i',x),1:numel(info.DeviceIDs),'UniformOutput',false);
                    comp.NotificationsLabel.Text = '';
                end
                comp.adaptor = comp.AdaptorListBox.Value;
                return
            end
            comp.adaptor = comp.AdaptorListBox.Value;
            comp.device  = comp.DeviceIDListBox.Value;

            comp.ChangesettingsButton.Enable = 'on';
            if numel(comp.device)>1
                % multi-select on
                comp.ChangesettingsButton.Enable = 'off';
            end
            
            notify(comp, 'SelectionChanged');
        end

        % Callback function
        function toggleSettings(comp, event)
            
            switch event.Source
                case comp.openButton
                    newsize = {'1x',20,'1x'};
                case comp.closeButton
                    newsize={'1x',20,0};
            end
            
            comp.GridLayout.ColumnWidth = newsize;
        end

        % Callback function
        function changeTime(comp, event)
            comp.time = event.Source.time;
        end

        % Button pushed function: ConnectButton
        function connectCallback(comp, event)
            switch comp.ConnectButton.UserData
                case true
                    warning off;imaqreset,warning on
                    comp.ConnectButton.Text = 'Connect';
                    comp.ConnectButton.UserData = false;
                    comp.AdaptorListBox.Enable = 'on';
                    comp.DeviceIDListBox.Items = comp.DeviceIDListBox.UserData;
                otherwise
                    comp.ConnectButton.Text = 'Disconnect';
                    comp.ConnectButton.UserData = true;
                    warning off
                    cameras         = cellfun(@(ids)videoinput(comp.adaptor,ids),comp.device,'UniformOutput',false);
                    comp.camera     = [cameras{:}];
                    comp.labCamera  = arrayfun(@(idx) labCamera(comp.camera(idx)),1:size(comp.camera,2)); %#ok<CPROPLC> 
                    warning on
                    comp.DeviceIDListBox.UserData = comp.DeviceIDListBox.Items;
                    comp.DeviceIDListBox.Items = comp.DeviceIDListBox.Value;
                    comp.AdaptorListBox.Enable = 'off';
            end


        end

        % Callback function: ScheduleDate, ScheduleDropDown, ScheduleTime
        function scheduleChanged(comp, event)
            switch comp.ScheduleDropDown.Value
                case 'From now'
                    comp.ScheduleDate.Enable = "off";
                    comp.SheduledTime=datetime('now')+comp.ScheduleTime.SumTime;
                    if strcmp(comp.sheduleTimer.Running,'off')
                        start(comp.sheduleTimer)
                    end
                case 'At time'
                    stop(comp.sheduleTimer)
                    comp.ScheduleDate.Enable = "on";
                    comp.SheduledTime=comp.ScheduleDate.Value+comp.ScheduleTime.SumTime;
                    comp.timeSheduleLabel.Text = sprintf('Set for %s',string(comp.SheduledTime,'dd-MMM-uuuu HH:mm:ss'));
            end
        end

        % Button pushed function: ChangesettingsButton
        function openCameraSettings(comp, event)
            validcamera(comp).changeSettings;
        end

        % Button pushed function: PreviewButton
        function previewCallback(comp, event)
            arrayfun(@(cams) cams.seeLiveView,validcamera(comp))
        end

        % Callback function: IntervalTime
        function IntervalTimeChanged(comp, event)
            
        end
    end

    methods (Access = protected)
        
        % Code that executes when the value of a public property is changed
        function update(comp)
            % Use this function to update the underlying components
            
        end

        % Create the underlying components
        function setup(comp)

            comp.Position = [1 1 608 398];
            comp.BackgroundColor = [0.94 0.94 0.94];

            % Create GridLayout
            comp.GridLayout = uigridlayout(comp);
            comp.GridLayout.ColumnWidth = {136, '1x'};
            comp.GridLayout.RowHeight = {100, 130, '1x', '1x', 20};

            % Create CameraSettingsPanel
            comp.CameraSettingsPanel = uipanel(comp.GridLayout);
            comp.CameraSettingsPanel.BorderType = 'none';
            comp.CameraSettingsPanel.Title = 'Camera Settings';
            comp.CameraSettingsPanel.BackgroundColor = [0.902 0.902 0.902];
            comp.CameraSettingsPanel.Layout.Row = [1 4];
            comp.CameraSettingsPanel.Layout.Column = 2;

            % Create GridLayout2
            comp.GridLayout2 = uigridlayout(comp.CameraSettingsPanel);
            comp.GridLayout2.ColumnWidth = {'1x', '1x', '1x'};
            comp.GridLayout2.RowHeight = {20, 50, 50, 50, 20, 50, 48};
            comp.GridLayout2.Padding = [0 0 0 10];

            % Create IterationsTime
            comp.IterationsTime = uiDuration(comp.GridLayout2);
            comp.IterationsTime.MinutesValue = 10;
            comp.IterationsTime.Layout.Row = 5;
            comp.IterationsTime.Layout.Column = 2;

            % Create ScheduleTime
            comp.ScheduleTime = uiDuration(comp.GridLayout2);
            comp.ScheduleTime.SelectionChangedFcn = matlab.apps.createCallbackFcn(comp, @scheduleChanged, true);
            comp.ScheduleTime.Layout.Row = 1;
            comp.ScheduleTime.Layout.Column = 3;

            % Create IntervalTime
            comp.IntervalTime = uiDuration(comp.GridLayout2);
            comp.IntervalTime.SelectionChangedFcn = matlab.apps.createCallbackFcn(comp, @IntervalTimeChanged, true);
            comp.IntervalTime.Layout.Row = 1;
            comp.IntervalTime.Layout.Column = 2;

            % Create ChangesettingsButton
            comp.ChangesettingsButton = uibutton(comp.GridLayout2, 'push');
            comp.ChangesettingsButton.ButtonPushedFcn = matlab.apps.createCallbackFcn(comp, @openCameraSettings, true);
            comp.ChangesettingsButton.Layout.Row = 1;
            comp.ChangesettingsButton.Layout.Column = 1;
            comp.ChangesettingsButton.Text = 'Change settings';

            % Create ImageintervalLabel
            comp.ImageintervalLabel = uilabel(comp.GridLayout2);
            comp.ImageintervalLabel.Layout.Row = 2;
            comp.ImageintervalLabel.Layout.Column = 1;
            comp.ImageintervalLabel.Text = 'Image interval';

            % Create ScheduleLabel
            comp.ScheduleLabel = uilabel(comp.GridLayout2);
            comp.ScheduleLabel.Layout.Row = 5;
            comp.ScheduleLabel.Layout.Column = 1;
            comp.ScheduleLabel.Text = 'Schedule';

            % Create BreaktimeLabel
            comp.BreaktimeLabel = uilabel(comp.GridLayout2);
            comp.BreaktimeLabel.Layout.Row = 4;
            comp.BreaktimeLabel.Layout.Column = 1;
            comp.BreaktimeLabel.Text = 'Break time';

            % Create ScheduleDropDown
            comp.ScheduleDropDown = uidropdown(comp.GridLayout2);
            comp.ScheduleDropDown.Items = {'From now', 'At time'};
            comp.ScheduleDropDown.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @scheduleChanged, true);
            comp.ScheduleDropDown.Layout.Row = 5;
            comp.ScheduleDropDown.Layout.Column = 2;
            comp.ScheduleDropDown.Value = 'From now';

            % Create ScheduleDate
            comp.ScheduleDate = uidatepicker(comp.GridLayout2);
            comp.ScheduleDate.DisplayFormat = 'uuuu-MM-dd';
            comp.ScheduleDate.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @scheduleChanged, true);
            comp.ScheduleDate.Enable = 'off';
            comp.ScheduleDate.Layout.Row = 5;
            comp.ScheduleDate.Layout.Column = 3;

            % Create timeSheduleLabel
            comp.timeSheduleLabel = uilabel(comp.GridLayout2);
            comp.timeSheduleLabel.HorizontalAlignment = 'right';
            comp.timeSheduleLabel.VerticalAlignment = 'top';
            comp.timeSheduleLabel.Layout.Row = 7;
            comp.timeSheduleLabel.Layout.Column = [2 3];
            comp.timeSheduleLabel.Text = 'timeSheduleLabel';

            % Create GridLayout_2
            comp.GridLayout_2 = uigridlayout(comp.GridLayout);
            comp.GridLayout_2.RowHeight = {'1x', '1x', '1x', '1x'};
            comp.GridLayout_2.Padding = [0 0 0 0];
            comp.GridLayout_2.Layout.Row = 1;
            comp.GridLayout_2.Layout.Column = 1;

            % Create AdaptorListBoxLabel
            comp.AdaptorListBoxLabel = uilabel(comp.GridLayout_2);
            comp.AdaptorListBoxLabel.HorizontalAlignment = 'right';
            comp.AdaptorListBoxLabel.Layout.Row = 1;
            comp.AdaptorListBoxLabel.Layout.Column = 1;
            comp.AdaptorListBoxLabel.Text = 'Adaptor';

            % Create AdaptorListBox
            comp.AdaptorListBox = uilistbox(comp.GridLayout_2);
            comp.AdaptorListBox.Items = {};
            comp.AdaptorListBox.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @changeAttribute, true);
            comp.AdaptorListBox.Layout.Row = [1 4];
            comp.AdaptorListBox.Layout.Column = 2;
            comp.AdaptorListBox.Value = {};

            % Create GridLayout_3
            comp.GridLayout_3 = uigridlayout(comp.GridLayout);
            comp.GridLayout_3.RowHeight = {30, 60, '1x', 30, 30};
            comp.GridLayout_3.Padding = [0 0 0 0];
            comp.GridLayout_3.Layout.Row = [2 4];
            comp.GridLayout_3.Layout.Column = 1;

            % Create DeviceListBoxLabel
            comp.DeviceListBoxLabel = uilabel(comp.GridLayout_3);
            comp.DeviceListBoxLabel.HorizontalAlignment = 'right';
            comp.DeviceListBoxLabel.Layout.Row = 1;
            comp.DeviceListBoxLabel.Layout.Column = 1;
            comp.DeviceListBoxLabel.Text = 'Device';

            % Create DeviceIDListBox
            comp.DeviceIDListBox = uilistbox(comp.GridLayout_3);
            comp.DeviceIDListBox.Items = {};
            comp.DeviceIDListBox.Multiselect = 'on';
            comp.DeviceIDListBox.ValueChangedFcn = matlab.apps.createCallbackFcn(comp, @changeAttribute, true);
            comp.DeviceIDListBox.Layout.Row = [1 2];
            comp.DeviceIDListBox.Layout.Column = 2;
            comp.DeviceIDListBox.Value = {};

            % Create ConnectButton
            comp.ConnectButton = uibutton(comp.GridLayout_3, 'push');
            comp.ConnectButton.ButtonPushedFcn = matlab.apps.createCallbackFcn(comp, @connectCallback, true);
            comp.ConnectButton.Layout.Row = 4;
            comp.ConnectButton.Layout.Column = [1 2];
            comp.ConnectButton.Text = 'Connect';

            % Create PreviewButton
            comp.PreviewButton = uibutton(comp.GridLayout_3, 'push');
            comp.PreviewButton.ButtonPushedFcn = matlab.apps.createCallbackFcn(comp, @previewCallback, true);
            comp.PreviewButton.Layout.Row = 5;
            comp.PreviewButton.Layout.Column = [1 2];
            comp.PreviewButton.Text = 'Preview';

            % Create NotificationsLabel
            comp.NotificationsLabel = uilabel(comp.GridLayout);
            comp.NotificationsLabel.Layout.Row = 5;
            comp.NotificationsLabel.Layout.Column = [1 2];
            comp.NotificationsLabel.Text = 'notifications';
            
            % Execute the startup function
            startUpFcn(comp)
        end
    end
end