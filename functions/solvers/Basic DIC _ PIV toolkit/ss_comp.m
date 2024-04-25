function ss_comp(files,rec,im00)
index = 1;
if nargin<2
   rec = [];
   im00 = [];
end
if nargin<3
   im00 = [];
end

%% Generation
fig = uifigure('Name','Synthetic Schlieren',...
    'WindowScrollWheelFcn',@scroll);
gl = uigridlayout(fig,'RowHeight',...
    {'1x','1x','1x','1x','1x','1x','1x','1x'},...
    'ColumnWidth',{'1x','1x'});

uilabel(gl,'Text','index')
uislider(gl,'Limits',[1 numel(files.Files)],...
    'ValueChangedFcn',@updateIndex,...
    'Tag','scrollable')
blockspacingtt = 'resolution';
borderoverlaptt = 'blocksize = blockspacing + 2*borderoverlap';
maxdisptt = 'sets size of search window = blocksize + 2*maxdisp';
minqualitytt = 'minimum quality of the correlation to be a valid vector, NaN otherwise';
normcorrtt = 'use normalized cross correlation? relatively slow but robust';

uilabel(gl,'Text','resize')
imqual = uislider(gl,'Limits',[0 1],...
    'Value',1,'ValueChangingFcn',@parseSolver);
uilabel(gl,'Text','blockspacing')
bspacing = uislider(gl,'Limits',[0 30],...
    'Value',4,'ValueChangedFcn',@parseSolver,'Tooltip',blockspacingtt);
uilabel(gl,'Text','borderoverlap')
boverlap = uislider(gl,'Limits',[0 5],...
    'Value',4,'ValueChangedFcn',@parseSolver,'Tooltip',borderoverlaptt);
uilabel(gl,'Text','maxdisp')
mdisp = uislider(gl,'Limits',[0 30],...
    'Value',4,'ValueChangedFcn',@parseSolver,'Tooltip',maxdisptt);
uilabel(gl,'Text','minquality')
mqual = uislider(gl,'Limits',[0 0.01],...
    'Value',0.01,'ValueChangedFcn',@parseSolver,'Tooltip',minqualitytt);
normcorr = uibutton(gl,'state','Value',0,...
    'Text','normcorr','ValueChangedFcn',@parseSolver,'Tooltip',normcorrtt);

imagebut = uibutton(gl,'state',...
    'Text','Image','ValueChangedFcn',@imageMode);
gpubut = uibutton(gl,'state',...
    'Text','GPU','ValueChangedFcn',@parseSolver);


if isempty(im00)
    im00 = readimage(files,1);
end
if ~isempty(rec)
    im00 = imcrop(im00,rec);
end
im0 = im00;
imagesc(im00)
% dockfig
%%
    function updateIndex(~,source)
        index = round(source.Value);
%         disp('index active')
        parseSolver
    end
    function imageMode(~,~)
        if ~imagebut.Value
           parseSolver 
        end
        im = readimage(files,index);
        if ~isempty(rec)
            im  = imcrop(im,rec);
        end
        im = imresize(im,imqual.Value);
        imagesc(im)
        clc
    end
    function scroll(~,event)
        verticalScrollCount = event.VerticalScrollCount;
        hFound              = pointer2object(fig,{'scrollable'});
        if isempty(hFound)
            return
        end
        switch hFound.Type
            case 'uislider'
                val     = hFound.Value - verticalScrollCount;
                
                if val>=hFound.Limits(1) && val<=hFound.Limits(2)
                    hFound.Value = hFound.Value - verticalScrollCount;
                end
                index = round(hFound.Value);
                parseSolver
        end
        
    end
    function parseSolver(~,~)
%         disp('solver active')
        disp('Processing')
        if imagebut.Value
           imageMode
           return
        end
        
        im = readimage(files,index);
        if ~isempty(rec)
%         im0 = imcrop(im00,rec);
        im  = imcrop(im,rec);
        end
        im0 = imresize(im00,round(imqual.Value,1));
        im  = imresize(im,round(imqual.Value,1));
        
        %%
        if gpubut.Value
            im0 = gpuArray(im0);
            im = gpuArray(im);
        end
       
        [u,v,Xq,Yq,q] = solver(im0,im,...
            'blockspacing',round(bspacing.Value),...
            'blockoverlap',round(boverlap.Value),...
            'maxdisp',round(mdisp.Value),...
            'minquality',mqual.Value,...
            'normcorr',normcorr.Value);
        %%

%         imagesc(q)
%         return
        %%
%         tiledlayout(2,1)
%         nexttile
%         imagesc(im)
%         colormap gray
%         hold on
        u_mag = sqrt(u.^2+v.^2);
        
        
        U = smooth2a(gather(u),ceil(size(u,1)/200),ceil(size(u,2)/200));
        imagesc(q)
%         caxis([-10 10]/5)
%         cmocean('balance','pivot',0),colorbar
        
        
        return
        
        
        V = smooth2a(gather(v),ceil(size(u,1)/200),ceil(size(u,2)/200));
%         streamslice(U,V,2)
%         vis_flow(u,v,100, 1, 5, 'm',1);
%         quiver(Xq, Yq, u, v, 0, 'y')
%%      
        cla
        hold on
        %~isnan(vorticity(U,V)))
        cutoff = 2;
        imagesc(rescale(im,-1,1),'AlphaData',sqrt(u.^2+v.^2)<cutoff)
        imagesc(vorticity(U,V),'AlphaData',sqrt(u.^2+v.^2)>cutoff)
        hold off
        caxis(1.5.*[-1 1])
%         cmocean('balance','pivot',0),colorbar
%         hold off
%         caxis([0 10])
%         nexttile
%         imagesc(readimage(files,index))
%         drawrectangle(gca,'Position',rec,'InteractionsAllowed','none');
%         cmocean('balance','pivot',0)
%         caxis(30.*[-1 1])
    end
    function [u,v,Xq,Yq,q] = solver(Iref,Idef,varargin)
        %SS_DECOMP Summary of this function goes here
        %   Detailed explanation goes here
        
        bs = 4;   % resolution
        bo = 4;  % blocksize = blockspacing + 2*borderoverlap
        md = 4;        % sets size of search window = blocksize + 2*maxdisp
        mq = 0.01;  % minimum quality of the correlation to be a valid vector, NaN otherwise
        nc = 1;       % use normalized cross correlation? relatively slow but robust
        parseInput(varargin)
        
        if mq==0
            mq=nan;
        end
        
        switch class(Iref)
            case 'gpuArray'
                warpmethod = 'linear';
            otherwise
                warpmethod = 'cubic';
        end
        %%
        tic
        [u, v, br, bc, q] = dic_dispfield(Iref, Idef, bs, bo, md, [], mq, nc);
        flag = false;
        switch class(Iref)
            case 'gpuArray'
                try
                    u = medfilt2(u);
                    v = medfilt2(v);
                catch
                    u = medfilt2(gather(u),'symmetric');
                    v = medfilt2(gather(v),'symmetric');
                end
            otherwise
                u = medfilt2(u,'symmetric');
                v = medfilt2(v,'symmetric');
        end        
        
        % warp cycle(s) if required
        Iref_w = interpimwarp(Iref, u, v, bc, br, warpmethod);
        smalldisp = 2; % allow small displacement around current solution
        [du, dv, ~, ~, q] = dic_dispfield(Iref_w, Idef, bs, bo, md, smalldisp, mq, nc);
        u = u + du;
        v = v + dv;
        
        switch class(Iref)
            case 'gpuArray'
                try
                    u = medfilt2(u);
                    v = medfilt2(v);
                catch
                    flag = true;
                    u = medfilt2(gather(u),'symmetric');
                    v = medfilt2(gather(v),'symmetric');
                end
            otherwise
                u = medfilt2(u,'symmetric');
                v = medfilt2(v,'symmetric');
        end   
        
        % final interpolation for display
        [Xq, Yq] = meshgrid(1:size(Iref,2), 1:size(Iref,1));
        % u = interp2(bc,br,u,Xq,Yq,'linear',0);
        % v = interp2(bc,br,v,Xq,Yq,'linear',0);
        [Iref_w, u, v] = interpimwarp(Iref, u, v, bc, br, warpmethod);
        switch class(Iref)
            case 'gpuArray'
                [du, dv] = of_dispfield(Iref_w,Idef, .2);
            otherwise
                roi = ~isnan(u);
                [du, dv] = of_dispfield(Iref_w,Idef, .1, roi);
        end
        u = u + du;
        v = v + dv;
        clc
        if flag
            disp('Used CPU')
        end
        toc
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

