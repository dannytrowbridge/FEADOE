@echo off

set pdrive=%~d0
set pdrive=C:\DAN

set CALCULIX_CODE=dccx.exe 
set CALCULIX_ROOT=%pdrive%\DEV\pywork\FEADOE\calculix
set PATH=%CALCULIX_ROOT%;%PATH%

"%CALCULIX_ROOT%\timethis.bat" "%CALCULIX_ROOT%\%CALCULIX_CODE%" %1 %2 %3 %4 %5 %6
