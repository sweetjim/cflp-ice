function delayedMessage(message,delay)
%DELAYEDMESSAGE prints a message to the command window after a delay.
% Message   - [string]
% Delay     - [double,datetime]
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
    TasksToExecute=1, ...
    TimerFcn=@(~,~)fprintf('%s\r',message), ...
    StopFcn=@(t,~) delete(t)))
warning on

