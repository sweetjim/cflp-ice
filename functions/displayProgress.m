function varargout = displayProgress(string,val,start_val,end_val,varargin)
% DISPLAYPROGRESS is a command line waitbar as a percentage
% 
% Inputs:
%   string      - [str] Prefix before percentage output
%   val         - [int] Index of loop
%   start_val   - [int] Starting value of loop
%   end_val     - [int] End value of loop
%   'delete'    - (OPTIONAL) deletes the printed output
%
% Example:
% loop = 40:100;
% for i=loop
%   pause(.1) % example calculation time
%   displayProgress('Progress:',i,loop(1),loop(end),'delete')
% end
%%

if ~nargout
    if val==start_val
        fprintf('\r%s:\t',string)
    end
end
abdiff = @(i) diff([i end_val]);
frac = round((abdiff(start_val)-abdiff(val))/abdiff(start_val)*100);

if nargout>0
    varargout{1} = sprintf('%s = %i%%',string,frac);
    return
end

if frac == 0; frac=1;end
if val>start_val
    for j=0:log10(frac*1e3)
        fprintf('\b'); 
    end
end
lines = fprintf(' %i %%', frac);
pause(1e-5); 

if val==end_val
    if nargin>4
        if strcmp(varargin{1},'delete')
            fprintf(repmat('\b',1, lines+length(string)+1))
             fprintf('\r');
        else
            fprintf('\r');
        end
    else
        fprintf('\r');
    end
end
end

