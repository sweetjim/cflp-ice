function appPlotting(ax,varargin)

input       = 'analysis';
inputstruct = [];
cmap        = 'parula';
parseInput(varargin)
persistent time_units

switch input
    case 'analysis'
        
        H       = inputstruct.H;
        ds      = inputstruct.ds;
        t       = inputstruct.t;
        Z       = inputstruct.z;
        dt_str  = lower(inputstruct.dt_str);
        cam     = inputstruct.cam;
        
        if isempty(time_units)
            time_units  = dt_str;
            tplot       = t;
        else
            tplot = t;
            if ~strcmp(lower(time_units),lower(dt_str))
                tplot = converttime(t,time_units,dt_str);
            end
        end 
        
        switch lower(inputstruct.colormode)
            case {'height','thickness'}
                C       = H*ds;
                clabel  = '$h$ (m)';
            case {'gradient','dht'}
                [C,~]   = gradient(H*ds,tplot,Z);
                clabel  = sprintf('%s (m/%s)','$\partial_t h$',dt_str);
            case 'dhtt'
                [ht,~]  = gradient(H*ds,tplot,Z);
                [C,~]   = gradient(ht);
                clabel  = sprintf('%s (m/%s)','$\partial_{tt} h$',dt_str);
            case 'laplacian'
                C       = smooth2a(del2(H*ds,tplot,Z),5,1);
                clabel  = '$\nabla^2h$';
        end
        
        switch lower(inputstruct.mode)
            case 'waterfall'
                waterfall(ax,Z,tplot,H'*ds,C')
                view(ax,cam(1),cam(2))
            case 'surf'
                surf(ax,Z,tplot,H'*ds,C')
                shading(ax,'interp')
                view(ax,cam(1),cam(2))
            case {'contour','contoured'}
                %%
                clim = [0 max(H*ds)];
                surf(ax,Z,tplot,H'*ds,C')
                hold(ax,'on')
                caxis(ax,max(caxis(ax)).*[-1 1])
                %         cmocean(ax,'balance','pivot',0,2^5)
                
                
                contour3(ax,Z,tplot,H'*ds,'k',...
                    'LevelList',0:0.01:clim(2))
                hold(ax,'off')
                %         view(90,90)
                view(ax,cam(1),cam(2))
                shading(ax,'interp')
        end
        
        std2 = 2*std(C,[],'all');
        switch lower(inputstruct.colormode)
            case {'height','thickness'}
                clim = [0 Inf];
%             case 'laplacian'
%                 clim = 1e1.*[-1 1]
            otherwise
                clim = std2.*[-1 1];
        end
        
        addColorbar('ax',ax,'latex','label',clabel,'fs',15,...
            'cmap',cmap,'limits',clim,'location','northoutside')
        
        addlabels('ax',ax,'latex','fs',15,...
            'x','$z$ (m)',...
            'z','$h$ (m)',...
            'y',sprintf('%s (%s)','$t$',dt_str))
        
        axis(ax,'tight')
    case 'melt'
        
end

    function parseInput(varargin)
        m = 1;
        items = varargin{:};
        for k=1:length(items)
            switch items{m}
                %% Name arguments
                case inputclasses
                    input       = items{m};
                    inputstruct = items{m+1}; 
                    m = m+1;
                case colormapclasses
                    cmap        = items{m};
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

