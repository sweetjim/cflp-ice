function app_waveTankEmulator % interactive designer / tank+wave emulator
%% Figure
XLIM = 1.5;
YLIM = 1.1;

fig = uifigure('Name','Internal Wave Experiment Modeller',...
    'Resize','on');
fig.Position = [fig.Position(1:2) 750 420];
gl = uigridlayout('Parent',fig,'ColumnWidth',{'1x'},'RowHeight',{'1x',40});
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
h0          = 5e-2;
m0          = 8.5e-2;
Z           = linspace(0,YLIM,1e3);
h           = @(m0,z0) h0*exp(-((Z-z0)/m0).^2);
patch(ax,h(m0,YLIM/2),Z,'k','Tag','wavemaker');
m0line      = drawline('Parent',ax,'Label','','Tag','m0','Position',[0 YLIM/2+m0;0 YLIM/2-m0]);
%% Buttons
gl2     = uigridlayout(gl,'ColumnWidth',repmat({80},[1 8]),'RowHeight',{30});
uilabel(gl2,'Text','N (rad/s)');
N       = uieditfield(gl2,'numeric','Value',.3,'ValueChangedFcn',@updateLabel);
uilabel(gl2,'Text','Omega (rad/s)');
omega   = uieditfield(gl2,'numeric','Value',N.Value*sin(deg2rad(angle(source))));
uilabel(gl2,'Text','Ray length (m)');
raylength = uieditfield(gl2,'numeric','Limits',[0 100],'Value',0,'ValueChangedFcn',@drawRays,'UserData',false);
uilabel(gl2,'Text','','Tag','raylength');
harmonics = uibutton(gl2,'state','Text','Harmonics','Value',0,'ValueChangedFcn',@drawRays);
%% Rays
negsource           = drawline(ax,'Position',source.Position,'Color',...
    'k','Tag','negsource',...
    'InteractionsAllowed','none');
negsource.Position  = source.Position-[0 0;0 2*diff(source.Position(:,2))];
% cbar = addColorbar(ax,'cmap','thermal','limits',[0 100],'levels',10,'title','Wave amplitude (%)');

line(ax,[0 XLIM],[0 0],'Tag','floor-line','LineStyle','--','Color','b','SelectionHighlight','off')
line(ax,[0 XLIM],[YLIM YLIM],'Tag','ceil-line','LineStyle','--','Color','b','SelectionHighlight','off')
rayfloor    = drawpoint(ax,'Position',[XLIM/2 0],'Tag','floor');
rayceil     = drawpoint(ax,'Position',[XLIM/2 YLIM],'Tag','ceil');
%% Listeners
addlistener(rayfloor,'MovingROI',@(varargin)lineChangingFcn(rayfloor));
addlistener(rayceil,'MovingROI',@(varargin)lineChangingFcn(rayceil));
addlistener(m0line,'MovingROI',@(varargin)waveMakerChangingFcn(m0line));
addlistener(ice,'MovingROI',@(varargin)lineChangingFcn(ice));
addlistener(source,'MovingROI',@(varargin)lineChangingFcn(source));
%% Initiate
k       = @(theta,m) m/tan(deg2rad(theta));%sqrt(m.^2/(sec(theta)-1));
k0      = k(deg2rad(angle(source)),m0);

calculateRayTrajectory
%% Functions
    function waveMakerChangingFcn(src)
        wm                  = findall(ax,'Tag','wavemaker');
        src.Position(:,1)   = 0;
        dif                 = abs(src.Position(:,2)-source.Position(1,2));

        switch sign(diff(dif))
            case -1
                m0   = min(dif);
            case 1
                m0   = max(dif);
        end
        k0                      = k(deg2rad(angle(source)),m0);
        src.Position(:,2)       = m0.*[1 -1]+source.Position(1,2);
        wm.Vertices(:,1)        = h(m0,source.Position(1,2));
        wm.Vertices([1 end],1)  = 0;
        m0                      = 1/m0;

        updateLabel
        src.LabelAlpha          = 0;
        src.Label               = sprintf('%.2f cm',1/m0*1e2);
        src.LabelTextColor      = 'r';
        cl                      = @() set(src,'Label','');
        delayedFunction(cl,5)
    end
    function lineChangingFcn(iwect)
        switch iwect.Tag
            case 'ice'
                iceth.setPosition([ice.Position(1) YLIM/3;sum(ice.Position([1 3])) YLIM/3]);
            case 'source'
                iwect.Position(1,1)=0;
                fixWavemaker
            case {'floor','ceil'}
                rayBounds
        end
        if source.Position(2,1)>=ice.Position(1)
            source.Position(2,1)=ice.Position(1);
        end
        k0      = k(deg2rad(angle(source)),m0);
        updateLabel
        updateRays

        function fixWavemaker
            wm = findall(ax,'Tag','wavemaker');
            M0 = findall(ax,'Tag','m0');
            zdiff   = abs(mean(M0.Position(:,2))-M0.Position(1,2));
            wm.Vertices(:,1)  = h(zdiff,source.Position(1,2));
            M0.Position(:,2) = source.Position(1,2)+[-1;1]*zdiff;
        end
        function updateRays
            negsource.Position = source.Position-[0 0;0 2*diff(source.Position(:,2))];
            ch = get(ax,'Children');
            delete(findall(ch,'Tag','intercept'))
            drawRays
        end
        function rayBounds
            src = source.Position(1,2);
            pos = iwect.Position(2);
            switch iwect.Tag
                case 'floor'
                    if pos>src
                        iwect.Position(2) = src;
                    end
                case 'ceil'
                    if pos<src
                        iwect.Position(2) = src;
                    end
            end
            
            set(findall(ax,'Tag',sprintf('%s-line',iwect.Tag)),YData = iwect.Position(2)*[1 1])
        end
    end
    function updateLabel(~,~)
        omega.Value = N.Value*sin(deg2rad(angle(source)));
        nu      = 1e-6;
        K       = sqrt(k0^2+m0^2);
        Alpha   = @(theta) nu/(2*N.Value)*K^3/sin(pi-theta);%nu/2*(k0./omega.Value).^3./sqrt(N.Value.^2.*(1-omega.Value.^2));
        Eloss   = @(theta,distance) exp(-Alpha(theta)*distance);

        n = 2;
        omega_n = omega.Value;
        while true
            if n>numel(raylength.UserData)-1
                break
            end
            if omega.Value/N.Value<1/n
                omega_n(end+1) = omega.Value*n; %#ok<AGROW> 
            else
                break
            end
            n = n+1;
        end

        omega_n = asin(omega_n./N.Value);
        try
            dists   = raylength.UserData(2:1+numel(omega_n));
            Eloss   = arrayfun(@(x,y) Eloss(x,y),omega_n,dists);
            Eloss_str = arrayfun(@(x) sprintf('%.2f',x),Eloss*1e2,'UniformOutput',false);
            Eloss_str = join(Eloss_str,',');
            source.Label = sprintf('%.1f deg (%.2fN) h=%.2f A=(%s)%s',...
                rad2deg(radangle(source)),...
                omega.Value/N.Value,...
                source.Position(1,2), ...
                Eloss_str{1},'%');
        catch
            source.Label = sprintf('%.1f deg (%.2fN) h=%.2f',...
                rad2deg(radangle(source)),...
                omega.Value/N.Value,...
                source.Position(1,2));
        end
    end
    function drawRays(~,~)
        raylength.Tag = 'normal';
        if raylength.Value==0
            raylength.Tag = 'intercept';
        end
        calculateRayTrajectory
        updateLabel
        % lineChangingFcn(source)
    end
    function calculateRayTrajectory
        %%
        theta   = deg2rad(angle(source));
        pts     = 2e5;
        cutoff  = raylength.Value;
        if strcmp(raylength.Tag,'intercept')
            cutoff = 10;
        end
        maxX    = XLIM-ice.Position(3);
        origin  = source.Position(1,2);
        
        raylength.UserData = cutoff;
        if theta==pi/2
            return
        end
        [X0,Y0,dis] = getTrajectory(theta);

        lbl = findall(fig,'Tag','raylength');
        lbl.Text = '';
        if strcmp(raylength.Tag,'intercept')
            set(lbl,'Text',sprintf('%.2fm',max(dis)))
        end
        
        omega.Value = N.Value*sin(deg2rad(angle(source)));
        nu      = 1e-6;
        K       = sqrt(k0^2+m0^2);
        Alpha   = @(theta) nu/(2*N.Value)*K^3/sin(pi-theta);
        Eloss   = exp(-Alpha(deg2rad(angle(source)))*dis)*1e2;

        warning off
        rayplot     = findall(ax,'Tag','ray');
%         dissipation = findall(ax,'Tag','loss');
        if isempty(rayplot)
%             if isempty(X0)
%                 X0      = source.Position(:,1)';
%                 Y0      = [source.Position(:,2)'; negsource.Position(:,2)']';
%                 dis     = hypot( ...
%                 mean(diff(X0))*cumsum(abs(gradient(X0,mean(diff(X0))))), ...
%                 mean(diff(Y0(:,1)))*cumsum(abs(gradient(Y0(:,1),mean(diff(Y0(:,1))))))');
%                 Eloss   = exp(-Alpha(deg2rad(angle(source)))*dis)*1e2;
%             end
%             sur(1)=coloredLine(X0,Y0(:,1),Eloss,ax);
%             sur(2)=coloredLine(X0,Y0(:,2),Eloss,ax);
%             set(sur,SelectionHighlight='off',Tag='loss');
            hold(ax,'on')
            plot(ax,X0,Y0,':r','Tag','ray','SelectionHighlight','off')
            plot(ax,X0,Y0,'--r','Tag','ray','SelectionHighlight','off')
            plot(ax,X0,Y0,'-r','Tag','ray','SelectionHighlight','off')
            hold(ax,'off')
        else
            %set(dissipation(1),'XData',[X0;X0]','YData',[Y0(:,1) Y0(:,1)],'CData',[Eloss; Eloss]','Zdata',[X0;X0]'*0)
            %set(dissipation(2),'XData',[X0;X0]','YData',[Y0(:,2) Y0(:,2)],'CData',[Eloss; Eloss]','Zdata',[X0;X0]'*0)
            set(rayplot(1),XData=X0,YData=Y0(:,1))
            set(rayplot(2),XData=X0,YData=Y0(:,2))
            set(rayplot(3:end),'Visible','off')
            
            if harmonics.Value
                if omega.Value/N.Value<.5
                    [X,Y] = getTrajectory(asin(2*omega.Value/N.Value));
                    set(rayplot(3),XData=X,YData=Y(:,1))
                    set(rayplot(4),XData=X,YData=Y(:,2))
                    set(rayplot(3:4),'Visible','on')
                end
                if omega.Value/N.Value<1/3
                    [X,Y] = getTrajectory(asin(3*omega.Value/N.Value));
                    set(rayplot(5),XData=X,YData=Y(:,1))
                    set(rayplot(6),XData=X,YData=Y(:,2))
                    set(rayplot(4:6),'Visible','on')
                end
            end
           
        end 
        warning on

        % wavelength shading
%         Ys1 = Y0(:,1)+1/k0.*[-1 1];
%         Ys2 = Y0(:,2)+1/k0.*[-1 1];
%         Ys1(Ys1>YLIM)=YLIM;Ys1(Ys1<0)=0;
%         Ys2(Ys2>YLIM)=YLIM;Ys2(Ys2<0)=0;
%         delete(findall(ax,'Tag','shaded'))
        % p(1) = shadedErrorBarPoly(X0,Ys1(:,1),Ys1(:,2),'x','k',ax);
        % p(2) = shadedErrorBarPoly(X0,Ys2(:,1),Ys2(:,2),'x','k',ax);
        % set(p,'EdgeAlpha',0,'SelectionHighlight','off','Tag','shaded')

        children = ax.Children;
        [~,csort]=sort(arrayfun(@(x) (isequal(x,findall(children,'Type','Images','Tag','source'))),children),'descend');
        ax.Children = ax.Children(csort);

        function [X,Y,S] = getTrajectory(theta)
            %%
            minY    = rayfloor.Position(2);
            maxY    = rayceil.Position(2);
            X       = linspace(0,XLIM*100,pts);
            Z       = linspace(0,YLIM*100,pts);
            Xfreq   = tan(YLIM/XLIM)/(maxX*tan(theta));%1.125
            Yfreq   = YLIM^2/(maxY-minY);%*1.41*(YLIM/maxX)^-1.7;
            
            Xs      = rescale(sawtooth(X*Xfreq*2*pi,.5),0,maxX);
            Ys      = rescale(sawtooth(Z*Yfreq*2*pi,.5),minY,maxY);
            Ysn     = -Ys;
            idx     = find(Ys>=origin,1);
            idx2    = find(Ys>=maxY-origin,1);
            Xs      = Xs(1:end-idx);
            Ys      = Ys(1+idx:end);
            Ysn     = Ysn(1+idx2:end)+maxY+minY;
            Ysn     = Ysn(1+find(Ysn<origin,1):end);

            if numel(Ysn)>numel(Ys)
                Ysn     = Ysn(1:numel(Ys));
            else
                Ys     = Ys(1:numel(Ysn));
                Xs     = Xs(1:numel(Ysn));
            end
            % hide beyond cutoff length
%             s               = hypot(cumsum(Xs)*mean(diff(X)),cumsum(Ys)*mean(diff(Z)));
            s=hypot( ...
                mean(diff(X))*cumsum(abs(gradient(Xs,mean(diff(X))))), ...
                mean(diff(Z))*cumsum(abs(gradient(Ys,mean(diff(Z))))));

            if strcmp(raylength.Tag,'intercept')
                cutoff =  s(find(Xs>maxX*(1-.01),1,'first'));
            end
            try %#ok<*TRYNC> 
                raylength.UserData(end+1) = cutoff;
                idx             = find(s>cutoff,1);
                Xs(idx:end)     = [];
                Ys(idx:end)     = [];
                Ysn(idx:end)    = [];
            end
            S = s(s<cutoff);
            X = Xs;
            Y = [Ys' Ysn'];
%             cutoff
        end
    end
end
