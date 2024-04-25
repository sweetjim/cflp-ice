function varargout = dependencies(filename)
%% Dependencies
% Shows all child functions of the file entered.
%
[fList,pList] = matlab.codetools.requiredFilesAndProducts(filename);
if nargout==0
fList'
{pList.Name}'
return
end

varargout{1} = fList;
varargout{2} = pList;
end

