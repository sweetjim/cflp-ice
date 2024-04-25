function errorHandler(ME)
% Catch new inputs that saved files do not already have

for ii=1:length(ME.stack)
    errorLine   = ME.stack(ii).line;
    errorFile   = ME.stack(ii).file;
    file        = splitlines(fileread(errorFile));
    errorCode   = strrep(file{errorLine},' ','');
    
    if contains(errorFile,'PhD')
        break
    end
end
matlab.desktop.editor.openAndGoToLine(errorFile,errorLine);
error('%s\nLine: %i\nCode: %s\n\n',ME.message,errorLine,errorCode)
end