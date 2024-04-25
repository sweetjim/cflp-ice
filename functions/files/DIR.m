function dirinfo = DIR(path)
%DIR is identical to dir but it has no files or directory tags ('.','..')
dirinfo = dir(path);
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
tf = ismember( {dirinfo.name}, {'.', '..'});
dirinfo(tf) = [];  %remove current and parent directory.
end

