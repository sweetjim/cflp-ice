function dist_map = getLayerThickness(in,varargin)

[normalize,window,method,inverse,hypersample] = parseInput(varargin);

if inverse
    in(isnan(in)) = 0;
    in = ~in; 
end

%%
dist_map = zeros(size(in));
for k=1:size(in,2)
    kmin = k-floor(window/2);
    kmax = k+floor(window/2);
    
    if kmin<1;kmin=1;end
    if kmax>size(in,2);kmax=size(in,2);end
    
    if sum(~in(:,kmax))==0
        kmin = k;
        kmax = k;
    end
    
    dist_map(:,kmin:kmax) = bwdist(in(:,kmin:kmax),method);
        
    if (normalize)&&(sum(dist_map(:,kmin:kmax)>0,[1 2])>0)
        normalizer = nanmax(dist_map(:,kmin:kmax),[],'all');
        dist_map(:,kmin:kmax) =  dist_map(:,kmin:kmax)/normalizer;
    end
    
end

% imagesc(dist_map)
% pause(.1)
%%
function [normalize,window,method,inverse,hypersample] = parseInput(varargin)
        
        normalize   = 0;
        window      = 1;
        method      = 'euclidean';
        inverse     = 0;
        hypersample = 0;
        
        m = 1;
        items = varargin{:};
        for K=1:length(items)
            switch items{m}
                case {'euclidean','cityblock','chessboard','quasi-euclidean'}
                    method = items{m};
                case 'normalize'
                    normalize = 1;
                case {'window','sampling'}
                    window = namevalue;
                    if floor(size(in,2)/2)<window
                        window = floor(size(in,2)/2);
                    end
                case {'hypersample','oversampling','interp'}
                    hypersample = round(namevalue);
                    if hypersample<0;hypersample=0;end
                case {'invert','inverse'}
                    inverse = 1;
            end
            m = m+1;
            if m>length(items);break;end
        end
        function out = namevalue
            out = items{m+1};
            m   = m+1;
        end
    end
end

