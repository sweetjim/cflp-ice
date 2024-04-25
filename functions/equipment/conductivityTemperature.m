function Tout = conductivityTemperature(TV)
if ismatrix(TV)&&all([~iscolumn(TV) ~isvector(TV)])
    Tout = cell2mat(arrayfun(@(idx) conductivityTemperature(TV(:,idx)),1:size(TV,2),'UniformOutput',false));
    return
end
infmap = isinf(TV);
nanmap = isnan(TV);
TV(isinf(TV))=0;
TV(isnan(TV))=0;
cv = [-0.0976    0.7225   -4.7043   11.3297];%[-3.1948   10.6287];
Tout = polyval(cv,TV);
Tout(nanmap)=nan;
Tout(infmap)=nan;
end

