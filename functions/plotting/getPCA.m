function [V,score,s2] = getPCA(x1,x2)
X       = [flatten(x1);flatten(x2)];
npts    = numel(x1);
Xavg    = mean(X,2,'omitnan'); % Compute mean
B       = X - Xavg*ones(1,npts); % Mean-subtracted Data
[U,S,V] = svd(B/sqrt(npts),'econ'); % PCA via SVD
theta   = (0:.01:1)*2*pi;
Xstd    = U*S*[cos(theta); sin(theta)];
for ij=1:3
    plot(Xavg(1)+ij*Xstd(1,:),Xavg(2) + ij*Xstd(2,:),'r-')
end
[V,score,s2]=pca(X);
end