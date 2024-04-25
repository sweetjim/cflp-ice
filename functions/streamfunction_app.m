function streamfunction_app(u,v)
% igure;
clf
tiledlayout(1,3,'TileSpacing','compact');
for i=1:3;n(i)=nexttile;end
im(1)=imagesc(n(1),u./std(u,[],'all'));
im(2)=imagesc(n(2),v./std(v,[],'all'));
im(3)=imagesc(n(3),u);
addColorbar(n,'cmap','balance','pivot',0,'limits',[-1 1]*5,'latex')%,'title',{'$u/\sigma_u$','$v/\sigma_v$','$\Psi/\sigma_\Psi$'},
pgon = drawrectangle('Parent',n(1),'FaceAlpha',0);
pgonv = drawrectangle('Parent',n(2),'FaceAlpha',0,'Position',pgon.Position,'InteractionsAllowed','none');
addlistener(pgon,'MovingROI',@movingROI);
    function movingROI(pgon,~)
        pgonv.Position = pgon.Position;
        val = cumsum(u.*createMask(pgon),'reverse');
%         val = flowfun(u.*createMask(pgon),v.*createMask(pgon));  

        set(im(3),'CData',val.*createMask(pgon))
        caxis(n(3),[-1 1]*1e-2)%std(val,[],'all'))
    end
end

