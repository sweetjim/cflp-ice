function out = velocity_field(imref,imdef,method,alpha,iterations,lambda,illumination,filter,diffmethod)
%% Velocity field calculator from Optical Flow methods
% Parameters:
%   imref   -   Reference image (2d double matrix)
%   imdef   -   Deformed image  (2d double matrix)
%   method  -   "fast" (Horn-Schunck alg.) or "slow" (Liu-Shen alg.)
%   alpha   -   a parameter that reflects the influence of the smoothness term.
%   iterations - (Horn-Shunck) number of iterations for velocity calculation
%
% Output:
%   out     -   structure containing:
%           im1     -   the reference image
%           im2     -   the deformed image
%           u       -   u velocity field
%           v       -   v velocity field
%           U       -   sqrt(u^2+v^2) field
%           vor     -   vorticity field (not great)
%% Default input arguments
if nargin<3
    method          = 'fast';
end
if nargin<4
    alpha           = 20;
end
if nargin<5
    iterations      = 2e2;
end
if nargin<6
    lambda          = 10;
end
if nargin<7
    illumination    = 0;
end
if nargin<8
    filter          = 1;
end
if nargin<9
    diffmethod      = 'hs';
end

if isempty(imref)
    error('No reference image data found')
end
if isempty(imdef)
    error('No deformed image data found')
end
%% Parse method
switch method
    case 'fast'
        [u,v] = HS(imref,imdef,alpha,iterations,diffmethod);
    case 'slow'
        [u,v] = LS(imref,imdef,alpha,lambda,filter,illumination);
    case 'lc'
        [u,v] = LucasKanadeOpticalFlow(imref,imdef);
end
%% Conversions
% calculate the velocity magnitude
u_mag   = sqrt(u.^2+v.^2);
u_max   = max(max(u_mag));
u_mag   = u_mag/u_max;
%%
% Spurious value omission
umean10 = mean(log10(u_mag),'all');
if umean10<-5
    u5std10 = 2*std2(log10(u_mag));
    omit = log10(u_mag)>umean10+u5std10;
    u_mag(omit)=nan;
    u(omit)=nan;
    v(omit)=nan;
end
% U=u;
% V=v;

% imagesc(U),colorbar
%%


% calculate vorticity
vor     = vorticity(u,v);
vor_max = max(max(abs(vor)));
vor     = vor/vor_max;

% calculate the 2nd invariant
Q       = invariant2_factor(gather(u),gather(v),1,1);

%% Output structure
out     = struct;
out.im1 = imref;
out.im2 = imdef;
out.u   = u;
out.v   = v;
out.U   = u_mag*u_max;
out.Q   = Q;
out.vor = vor;
end