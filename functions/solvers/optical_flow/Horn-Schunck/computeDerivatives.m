function [fx, fy, ft] = computeDerivatives(im1, im2,method)

if nargin<3
    method = 'hs';
end

if size(im2,1)==0
    im2=zeros(size(im1));
end
useconv     = @(im1,im2,kern) conv2(im1,kern,'same') + conv2(im2, kern,'same');
useimfilt   = @(im1,im2,kern) imfilter((im1+im2)/2, kern, 'symmetric',  'same');
kernt       = 1/4*ones(2);
switch method
    case 'hs'
        % Horn-Schunck original method
        kern = 1/4*...
               [-1 1
                -1 1];
        fx = useconv(im1,im2,kern);
        fy = useconv(im1,im2,kern');
        ft = conv2(im1,kernt,'same') + conv2(im2,-kernt,'same');
    case 'barron'
        % derivatives as in Barron (not advisable)
        kern = 1/12*...
            [-1 8 0 -1 1];
        fx= conv2(im1,kern,'same');
        fy= conv2(im1,kern','same');
        ft = conv2(im1,kernt,'same') + conv2(im2,-kernt,'same');
%         fx=-fx;
%         fy=-fy;
    case 'diff'
        % An alternative way to compute the spatiotemporal derivatives is to use simple finite difference masks.
        fx = conv2(im1,[1 -1],'same');
        fy = conv2(im1,[1; -1],'same');
        ft= im2-im1;
    case 'ls'
        % derivates as in Liu-Shen
        kern = 1/2*...
             [0  0  0
              0 -1 -1
              0  1  1];
        kernt = 1/4*...
             [0  0  0
              0  1  1
              0  1  1];
        %
        fx = useimfilt(im1,im2,kern');
        fy = useimfilt(im1,im2,kern);
        ft = imfilter(im2-im1, kernt, 'symmetric',  'same');
end
return
%%
% if nargin<3
%     method = 'hs';
% end
% 
% if size(im2,1)==0
%     im2=zeros(size(im1));
% end
% 
% switch method
%     case 'hs'
%         % Horn-Schunck original method
%         fx = conv2(im1,0.25* [-1 1; -1 1],'same') + conv2(im2, 0.25*[-1 1; -1 1],'same');
%         fy = conv2(im1, 0.25*[-1 -1; 1 1], 'same') + conv2(im2, 0.25*[-1 -1; 1 1], 'same');
%         ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
%     case 'barron'
%         % derivatives as in Barron
%         fx= conv2(im1,(1/12)*[-1 8 0 -8 1],'same');
%         fy= conv2(im1,(1/12)*[-1 8 0 -8 1]','same');
%         ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
%         fx=-fx;fy=-fy;
%     case 'diff'
%         % An alternative way to compute the spatiotemporal derivatives is to use simple finite difference masks.
%         fx = conv2(im1,[1 -1],'same');
%         fy = conv2(im1,[1; -1],'same');
%         ft= im2-im1;
%     case 'ls'
%         D1 = [0, 0, 0; 0,-1,-1;0,1,1]/2;
%         F1 = [0, 0, 0; 0,1,1;0,1,1]/4;
%         %
%         Ix = imfilter((I1+I2)/2, D1, 'symmetric',  'same');
%         Iy = imfilter((I1+I2)/2, D1', 'symmetric',  'same');
%         It = imfilter(I2-I1, F1, 'symmetric',  'same');
% 
% end