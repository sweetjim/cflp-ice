function [s,t,Sstate] = conductivityConvert(data,progress)
if nargin<2
    progress = false;
end

Sstate = true;
if ~isa(data,'struct')
    error('Input argument must be a structure containing "S" and "T" fields.')
end
if isempty(fieldnames(data))
    error('Input structure is empty')
end
if ~(any(contains(fieldnames(data),{'S','T'})))
    error('Input structure must contain both "S" and "T" fields')
end

load('labview_logs\probe\calibration\conductivityLookupTable.mat','lookup')
sref = lookup.S;
tref = lookup.T;
csref = lookup.CS;


T = data.T;
t = conductivityTemperature(T);

S = data.S;
% assume data is of (:,n) size
samples = 1:size(S,2);

Slu = zeros(size(S));
% locate indicies of C_S
for i=samples
    si      = S(:,i);
    ti      = t(:,i);
    % get T index
    tidx = arrayfun(@(idx) find(tref>ti(idx),1),1:numel(ti),'UniformOutput',false);
    sidx = arrayfun(@(idx) find(csref(:,tidx{idx})>si(idx),1),1:numel(ti),'UniformOutput',false,'ErrorHandler',@(~,~) {nan});
    SI = cellfun(@(sidx) sref(sidx),sidx,'UniformOutput',false);
    SI(cellfun(@isempty,SI))={nan};
    Slu(:,i) = cell2mat(SI);

    try %#ok<TRYNC>
        if isvalid(progress)
            progress.Notification = sprintf('Converting data (%.0f%s)',i/numel(samples)*1e2,'%');
            pause(1e-5)
        end
    end
end
s = Slu;
if all(isnan(Slu))
    s = S;
    Sstate = false;
end
end

