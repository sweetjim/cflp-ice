function showStreamfunction(u,v)
U = u;
V = v;
fac = max([max(abs(U),[],'all')...
    max(abs(V),[],'all')]);
U = U/fac;
V = V/fac;
[dUdx,dUdz]=gradient(U);
[dVdx,dVdz]=gradient(V);
clf
[~,psi]=flowfun(U,V);


tiledlayout(2,3)
nexttile(1)
imagesc(U)
addColorbar('cmap','balance','pivot',0,'title','U','location','southoutside')
nexttile(4)
imagesc(V)
addColorbar('cmap','balance','pivot',0,'title','V','location','southoutside')

nexttile(2,[2 1])
% imagesc(psi)
% hold on
contourf(psi)
axis ij
% vis_flow(U,V,'gx',100)
addColorbar('cmap','balance','pivot',0,'title','\Psi','location','southoutside')

nexttile
[dxpsi,dzpsi] = gradient(psi);
ur = -dzpsi;
vr = dxpsi;
imagesc(ur)
addColorbar('cmap','balance','pivot',0,'title','U','location','southoutside')
nexttile
imagesc(vr)
addColorbar('cmap','balance','pivot',0,'title','V','location','southoutside')
end