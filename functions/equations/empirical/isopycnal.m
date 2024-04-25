function [s,t] = isopycnal(rho0)
%ISOPYCNAL is a function that returns the range of salinities and
%temperatures that satisfy the input target density.
%%
sx = linspace(-5,50,1e3);
ty = linspace(-10,50,1e3);
[S,T] = meshgrid(sx,ty);
rho = density(S,T);
fig = figure('Visible','off');


if numel(rho0)>1
    si = zeros([numel(sx) numel(rho0)]);
    ti = si;
   for i=1:numel(rho0)
       [stmp,ttmp] = isopycnal(rho0(i));
       si(1:numel(stmp),i) = stmp;
       ti(1:numel(stmp),i) = ttmp;
   end
   si(si==0)=nan;
   s = si;
   t = ti;
   return
end
C = contourf(sx,ty,(real(rho)),'ShowText','on','LevelList',rho0);
s = C(1,2:end);
t = C(2,2:end);
t(s<0) = [];
s(s<0) = [];
close(fig)
end

