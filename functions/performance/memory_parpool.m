function memory_parpool(initialize)
persistent memory_gcp
%MEMORY_PARPOOL Summary of this function goes here
%   Detailed explanation goes here
if ~nargin
    initialize = false;
end

p = gcp('nocreate');

if initialize
    if ~isempty(p)
        delete(p)
    end
    memory_gcp = struct('without',readMemory,'new',[],'now',[]);
    parpool;
    memory_gcp.new = readMemory;
    memory_gcp.now = readMemory;
    outputs
    return
end

if isempty(p)
    return
end
memory_gcp.now = readMemory;
outputs

if ~nargin
    if memory_gcp.now>0.7
        fprintf('Memory leakage. Restarting parallel pool\n')
        delete(p)
        memory_gcp.without = readMemory;
        parpool;
        memory_gcp.now = readMemory;
    end
end

    function outputs
        fprintf('Memory allocation:\n\t%.2f %s (without parpool)',memory_gcp.without*1e2,'%')
        fprintf('\n\t%.2f %s (new parpool)',memory_gcp.new*1e2,'%')
        fprintf('\n\t%.2f %s (current parpool)\n',memory_gcp.now*1e2,'%')
    end
    function usage = readMemory
        tmp = memory;
        usage = tmp.MemUsedMATLAB/tmp.MaxPossibleArrayBytes;
    end
end

