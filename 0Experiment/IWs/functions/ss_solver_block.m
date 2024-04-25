function [u,v,Xq,Yq,q]=ss_solver_block(Iref,Idef,bs,bo,md,mq,nc,warpmethod,switcher)
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
    function filtering(switcher)
        if switcher
            try %#ok<*TRYNC>
                u = medfilt2(u);
                v = medfilt2(v);
                return
            end
        end
        u = medfilt2(gather(u),'symmetric');
        v = medfilt2(gather(v),'symmetric');
    end
end
