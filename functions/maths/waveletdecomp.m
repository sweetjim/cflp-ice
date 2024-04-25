function [A,H,V,D] = waveletdecomp(in,levels,method)
%% Wavelet decomposition
n       = levels; 
w       = method; 
[C,S]   = wavedec2(in,n,w);

[A,H,V,D] = waveletLevel(n);

% figure(1)
% tiledlayout(2,2)
% sgtitle(sprintf('Wavelet decomposition (n = %i, w = %s)',n,w))
% nexttile
% imagesc(A)
% title('Approximation')
% colorbar
% nexttile
% imagesc(H)
% title('Horizontal details')
% colorbar
% nexttile
% imagesc(V)
% title('Vertical details')
% colorbar
% nexttile
% imagesc(D)
% title('Diagonal details')
% colorbar

    function [dec,Ai,Hi,Vi,Di] = waveletLevel(level)
        Ai = appcoef2(C,S,w,level); % Approximation
        [Hi,Vi,Di] = detcoef2('a',C,S,level); % Details
        Ai = wcodemat(Ai,128);
        Hi = wcodemat(Hi,128);
        Vi = wcodemat(Vi,128);
        Di = wcodemat(Di,128);
        dec = [Ai Hi;Vi Di];
    end
end

