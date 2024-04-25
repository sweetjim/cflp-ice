function [f,ax] = iw_dissipation_distance
nu      = 1e-6;
Le      = 8e-2;
k0      = pi/(2*Le);
N       = linspace(0,.3,2e2+1)';
Omega   = linspace(0,1,2e2);
Alpha   =@(k0) nu/2*(k0./Omega).^3./sqrt(N.^2.*(1-Omega.^2));
x       = linspace(0,1e3,1e3);
[~,~,X] = meshgrid(Omega,N,x);

cutoff  = .9;
dis = calculateDistance(k0,cutoff);
%%
f   = uifigure('Name','Dissipation distance calculator');
f.Position = f.Position.*[1 1 0 1]+[0 0 725 0];
gl  = uigridlayout(f,'ColumnWidth',{700},'RowHeight',{'1x',70});
ax  = uiaxes(gl);
gl2 = uigridlayout(gl,'ColumnWidth',{100,200,150,200},'RowHeight',{50});

klabel      = uilabel(gl2,'Text','E-length (cm)');
kslider     = uislider(gl2,'Limits',[0 20],'Value',Le*1e2,'ValueChangingFcn',@updateLabel,'MajorTicksMode','auto','ValueChangedFcn',@updateAxes);
clabel      = uilabel(gl2,'Text','Attenutation cutoff (%)');
cutoffslider= uislider(gl2,'Limits',[0 1],'Value',cutoff,'ValueChangingFcn',@updateLabel,'ValueChangedFcn',@updateAxes);
imdata      = pcolor(ax,Omega,N,dis/1.5);
shading(ax,"interp")
addlabels(ax,'x','$\Omega$','y','$N$ (rad/s)','latex')
addColorbar(ax,'title','Attenuation distance (m)','cmap','thermal','levels',5)

ylim(ax,[0 max(N)])
xlim(ax,[0 max(Omega)])
line(ax,xlim(ax),0.26.*[1 1],'LineStyle','--','Color','r')
line(ax,xlim(ax),0.23.*[1 1],'LineStyle','--','Color','r')
arrayfun(@(x) line(ax,x.*[1 1],ylim(ax),'LineStyle','--','Color','r'),[.3 .4 .5 .7 .8])
% hold on
% contour(Omega,N,dis/1.5,'k','LevelStep',1)
caxis(ax,[0 5])
% imagesc(Omega,N,dis,'AlphaData',~isnan(dis))

    function dis = calculateDistance(k0,cutoff)
        dis     = NaN(numel(N),numel(Omega));
        alpha   = Alpha(k0);
        for i=1:numel(N)
            for j=1:numel(Omega)
                val=x(find(exp(-alpha(i,j)*x)<cutoff,1,'first'));
                if isempty(val)
                    continue
                end
                dis(i,j)=val;
            end
        end
    end
    function updateAxes(~,~)
        dis = calculateDistance(pi/(2*kslider.Value*1e-2),cutoffslider.Value);
        imdata.CData = dis;
    end
    function updateLabel(source,event)
        val = event.Source.Value;
        switch source
            case kslider
                klabel.Text   = sprintf('%s\n(%.2g)','E-length (cm)',val);
            case cutoffslider
                clabel.Text   = sprintf('%s\n(%.2g)','Attenutation cutoff (%)',val);
        end
    end

end