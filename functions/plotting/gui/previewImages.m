function [outfig,outax] = previewImages(images,cmap,clim)
%PREVIEWIMAGES generates an interactive figure where a time-series can be
%viewed for playback. Features include: start/stop functionality, refresh
%rate of image playback (period), c-axis limits (input same as normal), and
%drawing a region of interest (ROI).
%
% Inputs:
%   images  -   [3D matrix, double] time series array
%   cmap    -   [str] colormap (OPTIONAL) OR [1x2 double] caxis

%% Construction of figure elements
fig = uifigure('Name','Image time series preview', ...
    'CloseRequestFcn',@closeFig);
gl = uigridlayout(fig, ...
    "ColumnWidth",{'1x'}, ...
    'RowHeight',{'1x',100});
ax = uiaxes(gl);
gl2 = uigridlayout(gl, ...
    'ColumnWidth',{'1x','1x','1x','1x','1x','1x','1x'}, ...
    'RowHeight',{20,'1x','1x'});
imlabel  = uilabel(gl2,"Text",'Index (1)');
imslider = uislider(gl2, ...
    "Value",1, ...
    'Limits',[1 size(images,3)], ...
    'ValueChangedFcn',@sliderCallback,...
    'ValueChangingFcn',@sliderCallback);
imslider.Layout.Column = [2 numel(gl2.ColumnWidth)];
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
cm = uicontextmenu(fig);
uimenu(cm,"Text",'Rectangle','Checked','on','MenuSelectedFcn',@cmCallback)
uimenu(cm,"Text",'Polyline','Checked','off','MenuSelectedFcn',@cmCallback)
uimenu(cm,"Text",'Polygon','Checked','off','MenuSelectedFcn',@cmCallback)
uibutton(gl2, ...
    'Text','Draw ROI', ...
    'ButtonPushedFcn',@roiCallback,'ContextMenu',cm);
cropBut = uibutton(gl2,'state','Value',0, ...
    'Text','Activate', ...
    'Enable','off', ...
    'ValueChangedFcn',@roiCallback);

uidropdown(gl2,'Value','normal','Items',{'image','normal'},'ValueChangedFcn',@axContextMenu);
uilabel(gl2,'Text','Rolling (vert.)', ...
    'HorizontalAlignment','right');
shifty = uieditfield(gl2, ...
    'numeric','Value',0, ...
    'Limits',[-inf inf],'RoundFractionalValues','on');
shiftcountY = 0;
uilabel(gl2,'Text','Rolling (horz.)', ...
    'HorizontalAlignment','right');
shiftx = uieditfield(gl2, ...
    'numeric','Value',0, ...
    'Limits',[-inf inf],'RoundFractionalValues','on');
shiftcountX = 0;
rolling = uibutton(gl2,'state','Text','Toggle rolling','Value',0);
%% Initialization of image and timer
imdata = imagesc(ax,images(:,:,1));
axis(ax,'xy','tight')
if nargin<2
    cmap = 'gray';
end
if nargin<3
    clim = [-inf inf];
end
if isa(cmap,"double")
    clim = cmap;
    cmap = 'gray';
end
CLIM.Value = sprintf('[%i %i]',clim(1),clim(2));
colormap(ax,cmapCheck(cmap))
caxis(ax,clim)
roi = [];
k   = 1;
tr = timer('TimerFcn',@timerOperation,...
    'StopFcn',@stopFcn,...
    'TasksToExecute',inf,...
    'ExecutionMode','fixedRate',...
    'Period',.1);
roi_count = 1;
gbl = struct;
if nargout>=1
    outfig = fig;
    outax  = ax;
end
%% Nested functions (front-end)
    function timerOperation(~,~)
        k = k+1;
        if k>size(images,3)
            k = 1;
        end
        im = images(:,:,k);
        
        try %#ok<*TRYNC> 
            roitmp = findall(ax,'Type','Images');
            switch class(roitmp(1))
                case 'images.roi.Rectangle'
                    if cropBut.Value
                        im = imcrop(im,roitmp(1).Position);
                    end
                case 'images.roi.Polygon'
                    if cropBut.Value
                        im = imwarp(imcrop(im,gbl.crop),gbl.map);
                    end
            end
        end
        
        if rolling.Value
            try
                if shifty.Value~=0
                    im = circshift(im,shiftcountY,1);
                    shiftcountY = shiftcountY+shifty.Value;
                    if shiftcountY>size(im,1)
                        shiftcountY = shiftcountY-size(im,1);
                    end
                end
            end
            try
                if shiftx.Value~=0
                    im = circshift(im,shiftcountX,2);
                    shiftcountX = shiftcountX+shiftx.Value;
                    if shiftcountX>size(im,2)
                        shiftcountX = shiftcountX-size(im,2);
                    end
                end
            end
        end
%         im = imgaussfilt(im,1);
        imdata.CData    = im;
        ax.Title.String = sprintf('%i/%i',k,size(images,3));
    end
    function roiCallback(~,event)
        %delete(findall(ax,'Type','images.roi.Polygon'))
        if event.Source==cropBut
            im = images(:,:,k);
            set(roi,'Visible','on')
            if cropBut.Value
                roitmp = findall(ax,'Type','Images');
                switch class(roitmp(1))
                    case 'images.roi.Rectangle'
                        im = imcrop(im,roitmp(1).Position);
                    case 'images.roi.Polygon'
                        [map,crop]  = warpCoordinates(imdata.CData,roitmp(1));
                        delete(findall(ax,'Tag','cont'))
                        gbl.map     = map;
                        gbl.crop    = crop;
                        im          = imwarp(imcrop(im,crop),map);
                end
                set(roitmp,'Visible','off')
            end
            imdata.CData = im;
            return
        end
        roi_info        = struct('count',roi_count);
        roi_count       = roi_count+1;
        cropBut.Enable  = 'off';
        cropBut.Value   = 0;
        switch cm.Children([cm.Children.Checked]).Text
            case 'Rectangle'
                drawrectangle(ax,'FaceAlpha',0,'UserData',roi_info,'Rotatable',true);
                cropBut.Enable = 'on';
            case 'Polygon'
                pgon = drawpolygon(ax,'FaceAlpha',0,'UserData',roi_info);
                cropBut.Enable = 'on';
                addlistener(pgon,'ROIMoved',@pgonTransform);
            case 'Polyline'
                drawpolyline(ax,'UserData',roi_info);
        end
        roi = findall(ax,'Type','Images');
        tmp=arrayfun(@(x) createMask(x),roi,'UniformOutput',false);
        mask = sum(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),[]),3);
        assignin('base',"roi",roi)
        assignin('base',"mask",mask)
    end
    function pgonTransform(pgon,~)
        if cropBut.Value
            [map,crop,wf]  = warpCoordinates(imdata.CData,pgon,false);
            gbl.map     = map;
            gbl.crop    = crop;
            gbl.wp      = wf;
        else
            [mask,crop]  = warpCoordinates(imdata.CData,pgon,true);
            gbl.mask    = mask;
            gbl.crop    = crop;
            clim = caxis(ax);
            delete(findall(ax,'Tag','cont'))
            hold(ax,'on')
            contour(ax,mask,'r','Tag','cont')
            hold(ax,'off')
            caxis(ax,clim)
        end
    end
    %% Back-end
    function sliderCallback(~,event)
        value           = round(event.Value);
        imlabel.Text    = sprintf('Index (%i)',value);
        imdata.CData    = images(:,:,value);
        switch event.EventName
            case 'ValueChanged'
                event.Source.Value = value;
        end
    end
    function closeFig(~,~)
        stop(tr)
        delete(tr)
        delete(fig)
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
    function cmCallback(~,event)
        set(cm.Children,'Checked','off')
        event.Source.Checked='on';
    end
end