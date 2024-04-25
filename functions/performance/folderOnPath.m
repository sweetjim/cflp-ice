%% FOLDERONPATH 
% Attempts to convert string representations of folders or files paths into
% strings relative to the working directory. 
%   folderOnPath(root) - <root> [string,char,isFolder,isFile]
function root = folderOnPath(root)
if isempty(root)
    return
end
if iscell(root)
    root = cellfun(@folderOnPath,root,UniformOutput=false);
    return
end
if all(~[isa(root,'string') isa(root,'char')])
    error('Arguments must be string or char type')
end
if ~(isfolder(root)||isfile(root))
    error('Arguments must be valid folders or files')
end
label = split(root,'\');
onPWD = ismember(split(pwd,'\'),label);
if all(onPWD)
    root = join(label(~ismember(label,split(pwd,'\'))),'\');
end
if iscell(root)
    root=root{1};
end
end