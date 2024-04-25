function [V,Xstd,score,s2]=getPCA(x,y)
%% GETPCA 

%%
x       = check1d(x);
y       = check1d(y);
X       = [x;y];
npts    = numel(x);
Xavg    = mean(X,2); % Compute mean
B       = X - Xavg*ones(1,npts); % Mean-subtracted Data
[U,S,V] = svd(B/sqrt(npts),'econ'); % PCA via SVD
theta   = (0:.01:1)*2*pi;
Xstd    = U*S*[cos(theta); sin(theta)];
plot(x,y,'.')
hold(gca,'on')
for ij  = 1:3
    plot(gca,Xavg(1)+ij*Xstd(1,:),Xavg(2) + ij*Xstd(2,:),'r-')
end
hold(gca,'off')
[V,score,s2]=pca(X);

    function out = check1d(in)
        flatten = @(TwoDArray) reshape(TwoDArray',[1 size(TwoDArray,1)*size(TwoDArray,2)]);
        out = in;
        if ~any(size(in)-1==0)
            out = flatten(in);
        end
    end
end