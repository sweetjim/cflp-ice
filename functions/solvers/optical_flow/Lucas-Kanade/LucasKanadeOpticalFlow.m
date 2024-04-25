function [u,v] = LucasKanadeOpticalFlow(I1,I2,WindowSize,MaxIter,NumLevels)
% This	function computes optical flow using Iterative Pyramid Lucas-Kanade Algorithm	
% Input:    I1,I2 - images (I2 is back-projected to I1)	
%           WindowSize - size of local	neighbourhood	around	the	pixel
%           MaxIter - maximum iterations that we allow (for	each level of the pyramid)	
%           NumLevels - number of levels in	the	image	pyramid	
% Output:	[u,v] – Optical	flow warp parameters (velocity	fields)	

% ww = 45;
% w = round(ww/2);
% 
% % Lucas Kanade Here
% % for each point, calculate I_x, I_y, I_t
% Ix_m = conv2(I1,[-1 1; -1 1], 'valid'); % partial on x
% Iy_m = conv2(I1, [-1 -1; 1 1], 'valid'); % partial on y
% It_m = conv2(I1, ones(2), 'valid') + conv2(I2, -ones(2), 'valid'); % partial on t
% u = zeros(size(I1));
% v = zeros(size(I2));
% 
% % within window ww * ww
% for i = w+1:size(Ix_m,1)-w
%    for j = w+1:size(Ix_m,2)-w
%       Ix = Ix_m(i-w:i+w, j-w:j+w);
%       Iy = Iy_m(i-w:i+w, j-w:j+w);
%       It = It_m(i-w:i+w, j-w:j+w);
% 
%       Ix = Ix(:);
%       Iy = Iy(:);
%       b = -It(:); % get b here
% 
%       A = [Ix Iy]; % get A here
%       nu = pinv(A)*b; % get velocity here
% 
%       u(i,j)=nu(1);
%       v(i,j)=nu(2);
%    end
% end
% 
% % downsize u and v
% % u_deci = u(1:10:end, 1:10:end);
% % v_deci = v(1:10:end, 1:10:end);
% % % get coordinate for u and v in the original frame
% % [m, n] = size(im1t);
% % [X,Y] = meshgrid(1:n, 1:m);
% % X_deci = X(1:20:end, 1:20:end);
% % Y_deci = Y(1:20:end, 1:20:end);
% return

%%
I1=im2double(I1);
I2=im2double(I2);


% Decompose the image to gaussian pyramid
P1 = cell(NumLevels,1);
P2 = cell(NumLevels,1);
P1{1}= I1;
P2{1}= I2;
for lvl=2:NumLevels
    P1{lvl}=impyramid(P1{lvl-1},'reduce');
    P2{lvl}=impyramid(P2{lvl-1},'reduce');
end

% Initiazlie optical flow array
u=zeros(size(P1{NumLevels}));
v=zeros(size(P1{NumLevels}));

% Compute optical flow for each lvl pyramid,
% use the lower level as initial guess for the higer one.
for lvl=NumLevels:-1:1
    for j=1:MaxIter

        % Apply initial guess
        P2_wrapped = WarpImage(P2{lvl},u,v);
        
        % Compute flow field
        [du, dv] = LucasKanadeStep(P1{lvl},P2_wrapped, WindowSize);
        u = u + du; 
        v = v + dv;
    end
    
    % Upscale flow field
    if (lvl~=1)
        u = u * size(P1{lvl-1},1) / size(P1{lvl},1);
        v = v * size(P1{lvl-1},2) / size(P1{lvl},2);
        u = imresize(u,size(P1{lvl-1})); 
        v = imresize(v,size(P1{lvl-1}));
    end
end
end