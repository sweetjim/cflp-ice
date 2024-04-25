function outfig = imagepv_old(images,cmap)
%IMAGEPV generates an interactive figure where a time-series can be
%viewed for playback. Features include: start/stop functionality, refresh
%rate of image playback (period), c-axis limits (input same as normal), and
%drawing a region of interest (ROI).
%
% Inputs:
%   images  -   [3D matrix, double] time series array
%   cmap    -   [str] colormap (OPTIONAL)

%% Construction of figure elements
fig = uifigure('Name','Image time series preview', ...
    'CloseRequestFcn',@closeFig);
gl = uigridlayout(fig, ...
    "ColumnWidth",{'1x'}, ...
    'RowHeight',{'1x',75});
ax = uiaxes(gl);
gl2 = uigridlayout(gl, ...
    'ColumnWidth',{'1x','1x','1x','1x','1x','1x'}, ...
    'RowHeight',{'1x','1x'});
startbut=uibutton(gl2,'state', ...
    'Value',0, ...
    'Text','Start', ...
    'ValueChangedFcn',@startTimer);
uilabel(gl2,'Text','Period (s)', ...
    'HorizontalAlignment','right');
period = uieditfield(gl2, ...
    'numeric','Value',0.1, ...
    'Limits',[1e-6 inf],...
    'ValueChangedFcn',@changePeriod);
uilabel(gl2,'Text','C-axis', ...
    'HorizontalAlignment','right');
CLIM = uieditfield(gl2, ...
    'Text','Value','[-inf inf]', ...
    'Tooltip',{'Written as if command: caxis(input).',...
    'Can also use "images" as the time-series variable for eval input.',...
    'e.g., [-1 1]*std(images,[],"all")'},...
    'ValueChangedFcn',@caxisCallback);
uibutton(gl2, ...
    'Text','Draw ROI', ...
    'ButtonPushedFcn',@roiCallback);
uidropdown(gl2,'Value','normal','Items',{'image','normal'},'ValueChangedFcn',@axContextMenu);
uilabel(gl2,'Text','Rolling (vertical)', ...
    'HorizontalAlignment','right');
shift = uieditfield(gl2, ...
    'numeric','Value',0, ...
    'Limits',[-inf inf]);
shiftcount = 0;

%% Initialization of image and timer
imdata = imagesc(ax,images(:,:,1));
axis(ax,'xy','tight')
if nargin<2
    cmap = 'gray';
end
colormap(ax,cmapCheck(cmap))
caxis(ax,[-inf inf])
roi = [];
k   = 1;
tr = timer('TimerFcn',@timerOperation,...
    'StopFcn',@stopFcn,...
    'TasksToExecute',inf,...
    'ExecutionMode','fixedRate',...
    'Period',.1);
roi_count = 1;

if nargout==1
    outfig = fig;
end
%% Nested functions
    function closeFig(~,~)
        stop(tr)
        delete(tr)
        delete(fig)
    end
    function timerOperation(~,~)
        k = k+1;
        if k>size(images,3)
            k = 1;
        end
        im = images(:,:,k);
        if shift.Value~=0
            im = circshift(im,shiftcount,1);
            shiftcount = shiftcount+shift.Value;
            if shiftcount>size(im,1)
                shiftcount = shiftcount-size(im,1);
            end
        end
        imdata.CData    = im;
        ax.Title.String = sprintf('%i/%i',k,size(images,3));
    end
    function roiCallback(~,~)
        %delete(findall(ax,'Type','images.roi.Polygon'))
        roi_info = struct('count',roi_count);
        roi_count = roi_count+1;
        drawpolygon(ax,'FaceAlpha',0,'UserData',roi_info);
        roi = findall(ax,'Type','images.roi.Polygon');
        tmp=arrayfun(@(x) createMask(x),roi,'UniformOutput',false);
        mask = sum(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),[]),3);
        assignin('base',"roi",roi)
        assignin('base',"mask",mask)
    end
    function startTimer(~,event)
        switch event.Value
            case 1 
                event.Source.Text = 'Stop';
                start(tr)
            case 0
                stop(tr)
                event.Source.Text = 'Start';
        end
    end
    function changePeriod(~,~)
        stop(tr)
        tr.Period = period.Value;
        startbut.Text = 'Stop';
        startbut.Value = 1;
        start(tr)
    end
    function caxisCallback(~,~)
        val = CLIM.Value;
        try
            eval(sprintf('caxis(ax,%s);',val))
        catch
            val = str2double(split(extractBetween(val,'[',']'),{' ',','}));
            if isempty(val)
                return
            end
            val(isnan(val)) = inf;
            caxis(ax,val')
        end
        if strcmp(cmap,'balance')
            cmocean('balance','pivot',0,ax)
        end
    end
    function stopFcn(~,~)
        imdata.CData = images(:,:,k);
    end
    function axContextMenu(~,event)
        axis(ax,event.Value)
    end
end