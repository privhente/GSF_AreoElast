REM Batch file to run all DeSiO-Examples.
REM Results of calculations are stored in output_w.log
@echo off
set mypath="%cd%"
IF EXIST "output_w.log" (
	del "output_w.log" 
)
echo "Start running DeSiO-Aero - examples" %mypath% >> %mypath%\output_w.log
for /f "delims=" %%A in ('dir %mypath% /a:d /b') do (
	cd "%mypath%\%%A"
	for /f "delims=" %%B in ('dir %mypath%\%%A /a:d /b') do (
		cd %mypath%\%%A\%%B\DeSiO_InOut\
		del "*.bat"
		del "*.log"
		copy %mypath%\DeSiO.bat %mypath%\%%A\%%B\DeSiO_InOut
		echo %mypath%\%%A\%%B\DeSiO_InOut
		DeSiO.bat
		set a = 0
		for /f "delims== tokens=1,2" %%i in (check.log) do (
			echo %%A - %%B - %%i >> %mypath%\output_w.log
		)
	)
	echo( >> %mypath%\output_w.log
)
cd %mypath%
rem notepad output_w.log