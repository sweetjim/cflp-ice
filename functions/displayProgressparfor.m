function displayProgressparfor(wb,iterations,est_time)
%DISPLAYPROGRESSPARFOR
persistent count N h timeest tinit
 
if nargin==2
    h = wb;
    count = 0;
    tinit = tic;
    N = iterations;
    timeest = [];
elseif nargin==3
    h = wb;
    count = 0;
    tinit = tic;
    N = iterations;
    timeest = est_time; 
else
    count = count+1;
    if isempty(h)
        displayProgress('Progress',count,1,N,'delete')
        return
    end
    telapsed    = toc(tinit);
    progress    = count/N;
    rate        = telapsed/progress;
    esttime     = rate*(1-progress);

    str = sprintf('Processing (%s elapsed, %s remaining)', ...
        datestr(seconds(telapsed),'MM:SS'),...
        datestr(seconds(esttime),'MM:SS')); %#ok<DATST> 
    if isvalid(h)
        switch class(h)
            case 'matlab.ui.dialog.ProgressDialog'
                h.Value     = progress;
                h.Message   = str;
            otherwise
                waitbar(progress,h,str)
        end
    end
    try %#ok<TRYNC> 
        if h.CancelRequested
            close(h)
        end
    end
        
    if count==N
        delete(h)
    end
end

