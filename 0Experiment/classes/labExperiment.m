classdef labExperiment<handle
    %LABEXPERIMENT is a dynamic (handle) superclass that stores information
    %about a laboratory experiment.
    
    %% Experiment properties
    properties
        root string {mustBeText}=pwd
        dates datetime 
        S   (1,2) double {mustBeNumeric(S)}     = [0 0];
        rho (1,2) double {mustBeNumeric(rho)}   = density(0,0).*[1 1];
        N   (1,1) double {mustBeNumeric(N)}     = 0;
        T   (1,1) double {mustBeNumeric(T)}     = 0;
        eta (1,1) double {mustBeNumeric(eta)}   = 0;
    end
    properties 
        comments        string = '';
        savefile        string {mustBeText}
        inputs          struct
        parameterspace  table
        mixinglog       struct
        probelog        struct
        camera          struct
        driveSaveFile   string {mustBeText}
        driveUpdate     logical = false
    end
    properties (Constant,Hidden)
        savetag char = 'experimentdata'
        drivePath {mustBeFolder(drivePath)} = '0Experiment\experimentdata'
    end
    %% Internal properties
    properties (Access=private,Hidden)

    end

    %% Constructor
    methods
        function le = labExperiment(file2load)
            if nargin>=1
                le.loadExperiment(file2load)
            end
            le.savefile         = fullfile(le.root,le.savetag);
            le.driveSaveFile    = fullfile(getCloneFolder(le),le.savetag);
        end 
    end
    %% Recording scheduling
    methods (Hidden)
        function createSchedule(~)
            
        end
    end
    %% Saving
    methods
        function saveExperiment(le)
            % primary path
            save(le.savefile)
            % clone path (matlab drive)
            clonefolder = getCloneFolder(le);
            if ~isfolder(clonefolder)
                mkdir(clonefolder)
            end
            le.driveSaveFile    = fullfile(clonefolder,le.savetag);
            save(le.driveSaveFile)
        end
    end
    %% Loading
    methods
        function le=loadExperiment(le,file2load)
            switch class(file2load)
                case 'struct'
                    if all([isfield(file2load,'name') isfield(file2load,'folder')])
                        file2load = fullfile(file2load.folder,file2load.name);
                    end
                case {'string','char'}
                otherwise
                    return
            end
            if ~isfile(file2load)
                warning('Input is not a valid file')
                return
            end
            load(file2load,'le');
        end
        function load_labviewlog(le)

        end
    end
    %% Set properties
    methods
        function set.dates(le,val)
            if numel(val)==1
                val=[val val];
            end
            le.dates=val;
        end
        function set.S(le,val)
            le.S=val;
        end
        function set.rho(le,val)
            le.rho=val;
        end
        function set.N(le,val)
            le.N=val;
        end
        function set.T(le,val)
            le.T=val;
        end
        function set.root(le,val)
            le.root=val;
            updateSaveFile(le)
        end
    end
    %% Private methods
    methods (Access=private)
        function clonefolder = getCloneFolder(le)
            [~,rootfolder]      = fileparts(le.root);
            clonefolder         = fullfile(le.drivePath,rootfolder);
        end
        function updateSaveFile(le)
            le.savefile         = fullfile(le.root,le.savetag);
            le.driveSaveFile    = fullfile(getCloneFolder(le),le.savetag);
        end
    end
end

