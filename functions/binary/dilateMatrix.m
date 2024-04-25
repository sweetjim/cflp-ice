function matrix = dilateMatrix(type,n,varargin)
%% Dilate Matrix
% Generate a logical array.
% -------------------------------------------------------------------------
% %  Parameters:
% -------------------------------------------------------------------------
% 
% type: [char] (Required)
%   'line'
%   'cross'
%   'star'
%   'rect'
%   'circ'
%   'diamond'
% 
% size: [int, array] (Required)
%   Neighborhood size.  
% 
% feature: [int, double] (Optional)
%   Inclination, or approximation.
% 
% -------------------------------------------------------------------------
% % Outputs: (Optional)
% -------------------------------------------------------------------------
% 
% matrix: [logical, array]
%   The desired dilation matrix.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%

if nargin>2
    feature = varargin{1};
else 
    feature = 0;
end

type = lower(type);
switch type
    case 'none'
        matrix = 1;
        return
    case 'line'
        obj   = strel('line',n,feature);
        matrix = obj.Neighborhood;
        if mod(feature,90)==0&&length(matrix)~=n
            switch size(matrix,2)
                case 1
                    matrix = logical(ones(1,n))';
                otherwise
                    matrix = logical(ones(1,n));
            end
        return        
        end
    case 'cross'
        line1 = strel('line',n,feature);
        line2 = strel('line',n,-feature);
        matrix = imbinarize(line1.Neighborhood+line2.Neighborhood);
        return
    case 'star'
        line1   = strel('line',n,feature);
        line2   = strel('line',n,-feature);
        crosses = imbinarize(line1.Neighborhood+line2.Neighborhood);
        [H,L]   = size(crosses);
        space   = zeros(H,L);
        space(ceil(median(1:H)),:) = 1;
        space(:,ceil(median(1:L))) = 1;
        matrix  = imbinarize(crosses + space);
        return
    case {'rect','rectangle'}
%         if numel(feature)==1
%            feature = feature.*[1 1]; 
%         end
        obj   = strel('rectangle',[n feature]);
    case {'circ','circle'}
        try
            obj   = strel('disk',n,feature);
        catch
            obj   = strel('disk',n,0);
        end
    case 'diamond'
        obj   = strel('diamond',n);
end
matrix = obj.Neighborhood;

end

