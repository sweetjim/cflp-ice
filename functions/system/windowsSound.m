function y=windowsSound(label)
if ~nargin
    label = 'warning';
end

f = 2e4;
switch label
    case 'notify'
        y=audioread("C:\Windows\Media\Windows Background.wav");
    case 'ding'
        y=audioread("C:\Windows\Media\Windows Ding.wav");
    case 'error'
        y=audioread("C:\Windows\Media\Windows Critical Stop.wav");
        f= 5e4;
    case 'warning'
        y=audioread("C:\Windows\Media\Windows Error.wav");
        f= 5e4;
end
sound(y,f);
if ~nargout
    clear y
end
end

