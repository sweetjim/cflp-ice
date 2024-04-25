function t = scheduleTime(dd,hh,mm,ss)
%SCHEDULETIME creates a datetime variable relative to today.
% Inputs:
%   dd  -   [char,numeric] relative day from today  (default is "now")
%   hh  -   [numeric] relative hour of the day      (default is 0)
%   mm  -   [numeric] relative minute of the hour   (default is 0)
%   ss  -   [numeric] relative second of the minute (default is 0)

if nargin<1
    dd = datetime('now');
end
if nargin<2
    hh=0;
end
if nargin<3
    mm=0;
end
if nargin<4
    ss=0;
end
if isa(dd,"numeric")
    dd = datetime('today')+day(dd);
end
t=datetime(dd)+duration(hh,mm,ss);


