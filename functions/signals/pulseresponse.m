function [Y,fall,adjust,saturation] = pulseresponse(x,fs,y)
%PULSERESPONSE returns the average response time of the input signal to a
%pulsed change

if nargin<2||isempty(fs)
    fs = 1:numel(x);
end

pp = pulseperiod(x,fs);
tidx0   = [1; arrayfun(@(pt) find(fs>=pt,1,'first'),cumsum(pp))];
tidxf   = [tidx0(2:end)-1; numel(fs)];
Y = arrayfun(@(idx) y(tidx0(idx):tidxf(idx)),1:numel(tidx0),'UniformOutput',false);

maxlength = max(cellfun(@(Y) max(size(Y)),Y)); 

Y = cellfun(@(Y) interp1(fs, ...
                rescale(Y,0,1), ...
                linspace(fs(1),tidxf(1),maxlength))', ...
                Y,'UniformOutput',false);
Y = cell2mat(Y);
Y = mean(Y,2,'omitnan');

[~,lt,ut] = falltime(Y);

fs0         = fs(1:maxlength);
adjust      = fs0(find(fs0>=ut,1,'first'));
saturation  = fs0(find(fs0>=lt,1,'first'));
fall        = ut;

if ~nargout
    falltime(Y,fs0)
    ylabel('Normalized signal level')
    clear Y fall adjust saturation
end
end

