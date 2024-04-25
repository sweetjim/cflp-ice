classdef watchDirectory<handle
    %WATCHDIRECTORY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FileWatcher
        ValidEvents
        LastEvent
        LastTime
    end
    
    methods
        function wd = watchDirectory(root,filter,watchFor)
            if nargin<1
                root = uigetdir(pwd);
            end
            if nargin<2
                filter = '';
            end
            if nargin<3
                watchFor = 'all';
            end
            wd.FileWatcher = System.IO.FileSystemWatcher(root);
            wd.FileWatcher.Filter = filter;
            wd.FileWatcher.EnableRaisingEvents = true;
            switch watchFor
                case 'all'
                    wd.ValidEvents = events(wd.FileWatcher);
                    cellfun(@(event) addlistener(wd.FileWatcher,event,@eventhandlerChanged),wd.ValidEvents)
                otherwise
                    cellfun(@(event) addlistener(wd.FileWatcher,event,@eventhandlerChanged),watchFor)
            end
            function eventhandlerChanged(~,event)
                wd.LastTime     = datetime;
                wd.LastEvent    = event;
            end
        end
        
        
    end
end

