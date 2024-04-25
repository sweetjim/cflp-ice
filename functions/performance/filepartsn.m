function [filepath,name,ext] = filepartsn(filename,n,command)
%FILEPARTSN is a recursive fileparts function.
%
persistent keepfilename
if isempty(keepfilename)
    keepfilename = '';
end

if nargin<3
    command = 'reverse';
end

switch command
    case {'reverse','keep'}
        [filepath,name,ext]=fileparts(filename);
        if strcmp(command,'keep')
            keepfilename = fullfile(name,keepfilename);
        end
        if n>0
            [filepath,name,ext]=filepartsn(filepath,n-1,command);
            if strcmp(command,'keep')&&n>1
                keepfilename = fullfile(name,keepfilename);
            end
        else
            return
        end
    case 'forward'
        filepath = extractAfter(filename,filesep);
        name = '';
        ext = '';
        if n>0
            [filepath,name,ext]=filepartsn(filepath,n-1,'forward');
        else
            return
        end
end
if strcmp(command,'keep')
    name = keepfilename;
end

keepfilename = [];
end

