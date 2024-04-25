function boundary(S,T,z)
%% Boundary properties
% Generates a map of the boundary layer thickness 'eta' from
%
% eta = alpha(T_b-T_w)/(-beta dS0/dz)
%
% where alpha is the thermal expansion coefficient (function of T and S),
% beta the haline contraction coefficient (as before), and dS0/dz is the
% ambient solute gradient.
% -------------------------------------------------------------------------
% %  Singular arguments:
% -------------------------------------------------------------------------
%   S   :   [struct]    -   [g/kg]
%       'i'             -   ice salinity            [double]
%       'b'             -   interface salinity      [vec, double]
%       'w' or 'inf'    -   ambient salinity        [vec, double]
%
%   T   :   [struct]    -   [degC]
%       'i'             -   ice temperature         [vec, double]
%       'b'             -   interface temperature   ![double]!
%       'w' or 'inf'    -   ambient temperature     [vec, double]
%
%   z   :   [vec]       -   depth profile           [vec]
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Initialization
eta     = eta_map(S,T,z);
%% Plotting eta
pcolor(S.w,T.w,(eta))
shading interp
hold on
contour(S.w,T.w,abs(eta),'k','ShowText','on','LabelSpacing',1e3)
[C,h]=contour(S.w,T.w,eta,'b','ShowText','on','LevelList',[0],'LabelSpacing',5e2);
clabel(C,h,'Color','b');
contour(S.w,T.w,abs(eta),[0.1 .5],'k','ShowText','on','LabelSpacing',1e3)
hold off
axis xy

grid on

addColorbar('cmap','balance','pivot',0,...
    'latex','fs',15,...
    'title','$|\eta|$ (m)')
addlabels('latex',...
    'y','$T$ ($^\circ$C)',...
    'x','$d_zS_0$ (g/kg m)',...
    'fs',15)
set(gca,'XScale','log','layer','top')

% Tb tick mark
newtick(gca,T.b,'$T_b\rightarrow$','y')

% dS0/dz cr mark
idx = find(T.w>=T.b,1,'first');
Scr = S.w(find(sum(eta(idx.*[1 1]+[-1 1],:)',2)>0,1,'first'));
str = '\begin{tabular}{c} $\uparrow$\\$d_zS_{0,cr}$\end{tabular}';
newtick(gca,Scr,str,'x')
latexformat('ratio',[1 .75])
% out2latex('eta_map')
shg
end

