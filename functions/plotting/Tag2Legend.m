function L = Tag2Legend(ax,varargin)
%TAG2LEGEND Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
    ax = gca;
end
c           = get(ax,'Children');
[~,idx]     = unique(arrayfun(@(c) c.Tag,c,'UniformOutput',false),'last');
tags        = arrayfun(@(c) c.Tag,c(idx),'UniformOutput',false);
C           = c(idx);
last        = cellfun(@isempty,arrayfun(@(c) c.UserData,C,'UniformOutput',false));
L           = legend(C(last),tags(last),varargin{:});
if ~nargout
    clear L
end
end

