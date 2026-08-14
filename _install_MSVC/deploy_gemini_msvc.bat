@echo off

set "APPDIR=C:\Users\taras\Documents\Codex\Gemini\_install_MSVC"
set "WDEP=C:\Qt\6.7.0\msvc2019_64\bin\windeployqt.exe"

if not exist "%APPDIR%\Gemini.exe" (
  echo ERROR: Cannot find "%APPDIR%\Gemini.exe"
  echo Check APPDIR path.
  pause
  exit /b 1
)

if not exist "%WDEP%" (
  echo ERROR: Cannot find "%WDEP%"
  echo Check Qt path / version / kit folder.
  pause
  exit /b 1
)

cd /d "%APPDIR%" || (
  echo ERROR: cd failed to "%APPDIR%"
  pause
  exit /b 1
)

"%WDEP%" --release --compiler-runtime --dir . "Gemini.exe"

if errorlevel 1 (
  echo.
  echo ERROR: windeployqt failed.
  pause
  exit /b 1
)

echo.
echo DONE: Gemini++ MSVC deployment completed.
pause
