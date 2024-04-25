function delayedFunction(fcn_handle,delay,tag)
%DELAYEDMESSAGE executes a function after a delay.
% Fcn       - [function handle]
% Delay     - [double,datetime]
if nargin<3
    tag = 'messageTimer';
end
try stop(timerfindall('Tag',tag));end
switch class(delay)
    case 'datetime'
        delay = seconds(delay-datetime('now'));
    case 'duration'
        delay = seconds;
end
if delay<0
    delay = 0;
end
warning off
start(timer(StartDelay=delay, ...
    BusyMode="queue",...
    TasksToExecute=1, ...
    TimerFcn=@timerExecuteFcn, ...
    StopFcn=@(t,~) delete(t),...
    Tag=tag))
warning on

    function timerExecuteFcn(~,~)
        fcn_handle();
    end
end