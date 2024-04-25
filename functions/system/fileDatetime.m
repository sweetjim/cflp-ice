function [time,files,map] = fileDatetime(filename,mode)
%FILEDATETIME returns the datetime of the input file(s).
%   fileDatetime(filename) - Returns the datetime of the last modification 
%
%   fileDatetime(filename,'created') - Returns the datetime of the file
%                                      creation 
if ~any([isfolder(filename) isfile(filename)])
    error('Must be a valid file')
end

if ~iscell(filename)
    filename = {filename};
end
if nargin<2
    mode = 'modified';
end
if numel(filename)==1&&isfolder(filename)
    if strcmp(mode,'modified')
        folder = dir(filename{:});
        time = datetime(folder(1).date);
        return
    end
end
switch mode
    case 'created'
        if numel(filename)>1
            time = cellfun(@(filename) fileDatetime(filename,'created'),filename);
            return
        end
        d = System.IO.File.GetCreationTime(filename{:});
        time = datetime(d.Year, d.Month, d.Day, d.Hour, d.Minute, d.Second,d.Millisecond);
    otherwise
        folders = filename(isfolder(filename));
        if nargin<2
            mode = 'files';
        end
        if contains(mode,{'children','Children','files'})
            files = dir(folders{:});
            time  = datetime({files(~[files.isdir]).date});
            files = fullfile(files(1).folder,{files(~[files.isdir]).name});
            [time,map]  = sort(time);
            files       = files(map);
            return
        end

        files   = cellfun(@dir,filename(~isfolder(filename)));
        if ~isempty(files)
            filetime    = datetime({files.date});
        end
end
end

