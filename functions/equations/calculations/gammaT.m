function GammaT = gammaT(S,T)
%% Heat transfer coefficient GammaT
% Calculates the heat transfer coefficient from McConnochie & Kerr (2015)'s
% empirical, homogeneous theory 
% 
% gamma_T = Nu/Ra^(1/3) = 1/(k*Delta T) * (kappaT*nu)/(g alpha*Delta T)^(1/3)
% 
% Noting that alpha = alpha(Sw,Tw)
% -------------------------------------------------------------------------
% %  Singular arguments:
% -------------------------------------------------------------------------
%   S   :   [struct]    -   [g/kg]
%       'i'             -   ice salinity            [double]
%       'w' or 'inf'    -   ambient salinity        [vec, double]
%
%   T   :   [struct]    -   [degC]
%       'i'             -   ice temperature         [vec, double]
%       'b'             -   interface temperature   [vec, double]
%       'w' or 'inf'    -   ambient temperature     [vec, double]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Initialization
alpha   = sw_ALPHA(S.w,T.w);
beta    = sw_BETA(S.w,T.w);
deltaT  = (ones(size(T.w)).*T.b-T.w);

Tin     = mean([T.w T.b]);
kappaT  = property('water','kappaT','T',Tin);
nu      = property('water','nu','T',Tin);
k       = property('water','k','T',Tin);
g       = 9.8;

clear GammaT
GammaT  = 1./(k.*deltaT).*(kappaT.*nu./(g*alpha.*(deltaT))).^(1/3);


if numel(S.w)>1&&numel(T.w)==1
    %% S-vector
        plot(S.w,real(GammaT))
        addlabels('latex','fs',15,'x','$S (g/kg)$','y','Re$(\gamma_T)/q_T$')
elseif numel(S.w)==1&&numel(T.w)>1
    %% T-vector
        plot(T.w,real(GammaT))
        addlabels('latex','fs',15,'x','$S (g/kg)$','y','Re$(\gamma_T)/q_T$')
else
    %% ST-matrix or value
%     for i=1:length(T.w)
%         deltaT      = T.b-T.w(i);
%         GammaT(i,:)  = 1./(k.*deltaT).*(kappaT.*nu./(g*alpha.*(deltaT))).^(1/3);
%     end
        GammaT  = GammaT';
        Re_sign = sign(real(GammaT));
        Im_sign = sign(imag(GammaT));
        Re_g    = real(GammaT);
        clf
        pcolor(S.w,T.w,(Re_g.*Re_sign)),shading interp
        hold on
        contour(S.w,T.w,(Re_g.*Re_sign),...
            10.^[-6:-2],'k','ShowText','on','LabelSpacing',1e3);
        hold off
        addlabels('latex','fs',15,'x','$S$ (g/kg)','y','$T$ ($^\circ$ C)')
        addColorbar('latex','fs',15,...
            'title','Re$(\gamma_T)/q_T$',...
            'cmap','thermal','log','flip','limits',[1e-5 10])
        axis xy
        set(gca,'Layer','top','XScale','log')
        grid on
        % Tb tick mark
        newtick(gca,T.b,'$T_b\rightarrow$','y')

end
latexformat('ratio',[1 .75])
% out2latex('gammaT_map')
end

