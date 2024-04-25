function [ux,uy,vor,ux_horn,uy_horn,error1]=OpticalFlowPhysics_fun(I1,I2,lambda_1,lambda_2,method)
if nargin<5
    method = 'ls';
end
%%
% Horn's solution as an initial approximation of u and v
[Ix, Iy, It] = computeDerivatives(I1,I2,method);
%%
maxnum_1 = 500;
tol_1   = 10^(-12);
% [u0,v0]   = HS(I1,I2,lambda_1,maxnum_1);
[u,v] = horn_schunk_estimator(Ix, Iy, It, lambda_1, tol_1, maxnum_1);
ux_horn = v;
uy_horn = u;


% new model
Dm=1*10^(-3);
f=Dm*opt_laplacian(I1,1);

maxnum  = 60;
tol     = 1e-8;
% lambda_2 = 4000; % 2000 for dI0

dx=1;
dt=1; % unit time

% [u1,v1,error1] = liu_shen_estimator(I1, I2,...
%     f, dx, dt,...
%     lambda_2, tol, maxnum,...
%     u0, v0);

[u,v,error1] = liu_shen_estimator(I1, I2,...
    f, dx, dt,...
    lambda_2, tol, maxnum,...
    u, v);

ux=v;
uy=u;

vor=vorticity(ux, uy);
end


