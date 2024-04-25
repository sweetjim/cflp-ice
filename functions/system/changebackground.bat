echo off
REG ADD "HKCU\control panel\desktop" /v wallpaper /t REG_SZ /d "" /f
REG ADD "HKCU\control panel\desktop" /v wallpaper /t REG_SZ /d "C:\Users\James\Desktop\PhD\MATLAB Drive\backlight.bmp" /f
REG ADD "HKCU\control panel\desktop" /v WallpaperStyle /t REG_SZ /d 2 /f
RUNDLL32.EXE user32.dll,UpdatePerUserSystemParameters