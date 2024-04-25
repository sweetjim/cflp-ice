function changebackground(imagepath,style)
%CHANGEBACKGROUND changes the desktop background to the selected imagepath
%('.bmp' extension)
pause(1)
if ~isfile(imagepath)
    error('Cannot locate file')
end
[fpath,fname,ext] = fileparts(imagepath);
if ~strcmp(ext,'.bmp')
    imwrite(imread(imagepath),fullfile(fpath,strcat(fname,'.bmp')))
%     error('Image must be .bmp formatted')
end
if nargin<2
    style = 'fit';
end
% bash = 'reg add "HKEY_CURRENT_USER\Control Panel\Desktop" /v Wallpaper /t REG_SZ /d ';
% system(sprintf('%s%s%s',bash,imagepath,'/f'))
% system('RUNDLL32.EXE user32.dll,UpdatePerUserSystemParameters')
bash = {'echo off',...
    'REG ADD "HKCU\Control Panel\Desktop" /v Wallpaper /t REG_SZ /d "" /f',...
    strrep('REG ADD "HKCU\control panel\desktop" /v wallpaper /t REG_SZ /d "FPATH" /f',...
    'FPATH',imagepath),...
    'REG ADD "HKCU\control panel\desktop" /v WallpaperStyle /t REG_SZ /d 10 /f',...
    'RUNDLL32.EXE user32.dll,UpdatePerUserSystemParameters ',...
    'exit'}';
% switch style
%     case 'fit'
%     case 'fill'
%     case 
for i=1:numel(bash);system(bash{i});end
for i=1:numel('The operation completed successfully.\r\r')*3-3
fprintf('\b')
end
% end

