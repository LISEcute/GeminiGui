@echo off

set "APPDIR=C:\Users\taras\Documents\Codex\Gemini\_install"
set "APP=Gemini.exe"
set "WDEP=C:\Qt\6.7.0\mingw_64\bin\windeployqt.exe"

if not exist "%APPDIR%\%APP%" (
  echo ERROR: Cannot find "%APPDIR%\%APP%"
  echo Check APPDIR / APP path.
  pause
  exit /b 1
)

if not exist "%WDEP%" (
  echo ERROR: Cannot find "%WDEP%"
  echo Check Qt path / version / MinGW kit folder.
  pause
  exit /b 1
)

cd /d "%APPDIR%" || (
  echo ERROR: cd failed to "%APPDIR%"
  pause
  exit /b 1
)

echo Deploying Qt/MinGW runtime for "%APPDIR%\%APP%"...
"%WDEP%" --release --compiler-runtime --dir . "%APP%"

if errorlevel 1 (
  echo.
  echo ERROR: windeployqt failed.
  pause
  exit /b 1
)

echo.
echo DONE: Gemini++ MinGW deployment completed.
pause
