function Iref = ss_pattern(k0,savefile)
%SS_PATTERN Summary of this function goes here
%   Detailed explanation goes here
% s = diam * 72/25.4;   % diameter, in points units (1 point = 1/72 inch = 25.4/72 mm)
%a4 @ 600dpi = [4961 7016] px = [21 29.7] cm
[x,y] = meshgrid(0:4961-1,0:7016-1);

% desired resolution 2pi/k = 3.4px  in sensor
% 50mm lens, 3088x2064px, 3.5m dist, ~66 px/cm
pxpercm     = 66;
pxpercma4   = 4961/21;

lambdarec = 6;
if nargin<1
    k0 = 2*pi/(lambdarec*pxpercma4/pxpercm);
end
th = (pi/4)*.2;
xp = cos(th)*x + sin(th)*y;
yp = -sin(th)*x + cos(th)*y;
Iref = .5+(cos(k0*xp) + cos(k0*yp))/4;

if ~nargout
    h=gcf;
    % Window size:
    set(h,'Position',[360 80 560/sqrt(2) 560]);
    set(gca,'Position',[0 0 1 1]);   % location of the figure in the window
    set(gca,'PlotBoxAspectRatio',[1/sqrt(2) 1 1]);

    % Printing settings:
    set(h,'PaperUnits','centimeters');
    set(h,'PaperOrientation','portrait');
    set(h,'PaperType','A4');
    %set(h,'PaperSize',[21 29.7]);
    set(h,'PaperPosition',[0 0 21 29.7]);
    set(h,'PaperPositionMode','manual');
    set(h,'InvertHardcopy','off');   % keep the user background mode
    imshow(Iref),shg,colormap gray
    clear Iref
end
if nargin<2
    savefile = false;
end
if savefile
    print('-dtiff','-r600','ss_pattern');
end


