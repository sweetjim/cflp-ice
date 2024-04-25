function [u,v]  = LS(I1,I2,alpha,lambda,filter,illumination,kernel_iterations,usePIV)
%% https://doi.org/10.1017/S0022112008003273
%% DEFAULTS
if nargin<3
    error('Not enough input arguments.')
end
if nargin<4
    alpha = 20;
end
if nargin<5
    lambda = 2e3;
end
if nargin<6
    filter  = 6;
end
if nargin<7
    illumination = 0;
end
if nargin<8
    kernel_iterations = 1;
end
if nargin<9
    usePIV = false;
end


%%
% Set the Parameters for Optical Flow Computation

% Set the lagrange multipleirs in optical computation
lambda_1 = alpha;  % the Horn_schunck estimator for initial field
lambda_2 = lambda; % the Liu-Shen estimator for refined estimation

%% Preamble 
no_iteration    = kernel_iterations;    % Number of iterations in the coarse-to-fine iterative process from initial estimation, "1" means 1 iteration
scale_im        = 1;                    % Initial coarse field estimation in the coarse-to-fine iterative process, scale_im is a scale factor for down-sizing of images
size_average    = illumination;         % [pixels] For local illumination intensity adjustment, To bypass it, set size_average = 0
size_filter     = filter;               % [pixels] Gaussian filter size for removing random noise in images

% correcting the global and local intensity change in images
[m1,n1]         = size(I1);
window_shifting = [1;n1;1;m1]; % [x1,x2,y1,y2] deines a rectangular window for global intensity correction
[I1,I2]         = correction_illumination(I1,I2,window_shifting,size_average);
edge_width      = 1; % [pixels] cleaning the left and upper edges since some data near the edges are corrupted due to interperlation

% pre-processing for reducing random noise,
% and downsampling images if displacements are large
if size_filter>0
    [I1,I2] = pre_processing_a(I1,I2,scale_im,size_filter);
end

% initial correlation calculation for a coarse-grained velocity field (ux0,uy0)
% ux is the velocity (pixels/unit time) in the image x-coordinate (from the left-up corner to right)
% uy is the velocity (pixels/unit time) in the image y-coordinate (from the left-up corner to bottom)

%% FFT cross-correlation algorithm
if usePIV
    pivPar.iaSizeX = [64 16 8];         % size of interrogation area in X
    pivPar.iaStepX = [32 8 4];          % grid spacing of velocity vectors in X

    pivPar.ccMethod = 'fft';
    warning off
    [pivData1] = pivAnalyzeImagePair(I1,I2,pivPar);
    warning on

    ux0=pivData1.U;
    uy0=pivData1.V;

    % re-size the initial velocity field (u0, v0)
    [n0,m0] = size(ux0);
    [n1,m1] = size(I1);
    scale   = round((n1*m1/(n0*m0))^0.5);
    ux0     = imresize(ux0,scale);
    uy0     = imresize(uy0,scale);
else
    [ux0,uy0]=OpticalFlowPhysics_fun(I1,I2,lambda_1,lambda_2);
end
%% Iterative correction
% generate the shifted image from Im1 based on the initial coarse-grained velocity field (ux0, uy0),
% and then calculate velocity difference for iterative correction
% estimate the displacement vector and make correction in iterations

ux  = ux0;
uy  = uy0;

k=1;
while k<=no_iteration
    [Im1_shift,uxI,uyI]=shift_image_fun_refine_1(ux,uy,I1,I2);
    I1=double(Im1_shift);
    I2=double(I2);
    % calculation of correction of the optical flow
    [dux,duy]=OpticalFlowPhysics_fun(I1,I2,lambda_1,lambda_2);

    % refined optical flow
    ux_corr=uxI+dux;
    uy_corr=uyI+duy;


    k=k+1;
end

% refined velocity field
ux = ux_corr;    %%%%%
uy = uy_corr;    %%%%%

% clean up the edges
ux(:,1:edge_width)=ux(:,(edge_width+1):(2*edge_width));
uy(:,1:edge_width)=uy(:,(edge_width+1):(2*edge_width));

ux(1:edge_width,:)=ux((edge_width+1):(2*edge_width),:);
uy(1:edge_width,:)=uy((edge_width+1):(2*edge_width),:);

u = ux;
v = uy;
end