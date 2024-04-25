function Sout = conductivitySalinity(SV)
if ismatrix(SV)&&all([~iscolumn(SV) ~isvector(SV)])
    Sout = cell2mat(arrayfun(@(idx) conductivityTemperature(SV(:,idx)),1:size(SV,2),'UniformOutput',false));
    return
end
infmap = isinf(SV);
nanmap = isnan(SV);
SV(isinf(SV))=0;
SV(isnan(SV))=0;
cv = [-0.0468    0.4660    8.7402   26.2790];%[0.0001   -0.0063    0.3203   -4.9768]; % at 22.71pm0.27 degC
Sout = polyval(cv,SV);
Sout(nanmap)=nan;
Sout(infmap)=nan;
end
