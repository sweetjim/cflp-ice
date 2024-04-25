%% have a look over a range of T/S
s=10:50;
t=0:30;
[S,T]=ndgrid(s,t);

figure(1)
contourf(t,s,sw_dens(S,T,100),30);
colorbar
xlabel('T (deg. C)');
ylabel('S (psu)');
title('seawater density (kg/m^3)');

fname='../../Latex/cabbeling/sw_prop2.eps';
printeps(21,21,500,fname);

% pretty much linear with S, but not T
figure(2)
plot(t,sw_dens(S,T,100))

%% quadratic fits for the cabbeling problem

% mean pressure in db (evaluate at z=-H/2)
mp=100; % 100m

% reference salinity
S0=35; %psu

dT=0.05;
t=4:dT:23;
dS=0.1;
s=25:dS:45;
[T,S]=meshgrid(t,s);
rho=sw_dens(S,T,mp);

drhodT=interp2(T(:,1:end-1)+dT/2,S(:,1:end-1),diff(rho,1,2)/dT,T,S,'cubic',0);
d2rhodT2=interp2(T(:,1:end-1)+dT/2,S(:,1:end-1),diff(drhodT,1,2)/dT,T,S,'cubic',0);
drhodS=interp2(T(1:end-1,:),S(1:end-1,:)+dS/2,diff(rho,1,1)/dS,T,S,'cubic',0);


% plot
it=5:length(t)-4;
is=5:length(s)-4;

figure(3)

subplot(3,1,1)
contourf(T(:,it),S(:,it),d2rhodT2(:,it)./rho(:,it),30);
colorbar
title('cabbeling parameter c');
xlabel('T');
ylabel('S');

subplot(3,1,2)
contourf(T(:,it),S(:,it),drhodT(:,it)./rho(:,it),30);
colorbar
title('thermal expansion gammaT');
xlabel('T');
ylabel('S');

subplot(3,1,3)
contourf(T(is,:),S(is,:),drhodS(is,:)./rho(is,:),30);
colorbar
title('saline contraction gammaS');
xlabel('T');
ylabel('S');

%fname='../../Latex/cabbeling/sw_prop1.eps';
%printeps(21,21,500,fname);


% reference temps/salinity
% NP - Kuroshio
clear NP 
NP.T0=15; % 13 to 17
NP.S0=34.3;
NP.N2=1e-4;
NP.dT=4;
NP.f=9e-5;

% Indian -SAF
clear IN
IN.T0=6.5;
IN.S0=34.1;
IN.N2=3e-5;
IN.dT=5;
IN.f=1e-4;

% Gulf Stream
clear GS
GS.T0=18;%17.5;
GS.S0=36.05;%35.7;
GS.N2=4e-5;
GS.dT=5;
GS.f=9e-5;

% alpha*Kh
minaKh=0.01;
maxaKh=10;
g=9.81;


% For the North Pacific
NP.rho0=sw_dens(NP.S0,NP.T0,mp);
NP.gammaT=interp2(T,S,drhodT,NP.T0,NP.S0)/NP.rho0;
NP.gammaS=interp2(T,S,drhodS,NP.T0,NP.S0)/NP.rho0;
NP.c=interp2(T,S,d2rhodT2,NP.T0,NP.S0)/NP.rho0;

NP.maxWMT=sqrt(2*maxaKh*NP.f)*g*-NP.c*NP.dT^2/(8*NP.N2)*2*sqrt(2/pi);
NP.minWMT=sqrt(2*minaKh*NP.f)*g*-NP.c*NP.dT^2/(8*NP.N2)*2*sqrt(2/pi);

NP


% For the Indian
IN.rho0=sw_dens(IN.S0,IN.T0,mp);
IN.gammaT=interp2(T,S,drhodT,IN.T0,IN.S0)/IN.rho0;
IN.gammaS=interp2(T,S,drhodS,IN.T0,IN.S0)/IN.rho0;
IN.c=interp2(T,S,d2rhodT2,IN.T0,IN.S0)/IN.rho0;

IN.maxWMT=sqrt(2*maxaKh*IN.f)*g*-IN.c*IN.dT^2/(8*IN.N2)*2*sqrt(2/pi);
IN.minWMT=sqrt(2*minaKh*IN.f)*g*-IN.c*IN.dT^2/(8*IN.N2)*2*sqrt(2/pi);

IN

% For the Gulf Stream
GS.rho0=sw_dens(GS.S0,GS.T0,mp);
GS.gammaT=interp2(T,S,drhodT,GS.T0,GS.S0)/GS.rho0;
GS.gammaS=interp2(T,S,drhodS,GS.T0,GS.S0)/GS.rho0;
GS.c=interp2(T,S,d2rhodT2,GS.T0,GS.S0)/GS.rho0;

GS.maxWMT=sqrt(2*maxaKh*GS.f)*g*-GS.c*GS.dT^2/(8*GS.N2)*2*sqrt(2/pi);
GS.minWMT=sqrt(2*minaKh*GS.f)*g*-GS.c*GS.dT^2/(8*GS.N2)*2*sqrt(2/pi);

GS

save cab_data NP GS IN 

% rho min (at base of domain is)
%GS.rho0*(1+GS.N2/g*200)

% fit
%f = fittype('rho0+a*(x-b)^n','problem','n','options',s);
%[c2,gof2] = fit(cdate,pop,f,'problem',2)



