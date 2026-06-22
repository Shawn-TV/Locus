@echo off
setlocal
cd /d "%~dp0"

title Locus Review Editor
echo Starting Locus Review Editor...
echo.

set "BUNDLED_PY=%~dp0runtime\windows\python\python.exe"

if exist "%BUNDLED_PY%" (
  echo Using bundled Windows runtime. No Python installation is required.
  echo.
  "%BUNDLED_PY%" tools\launch_review_editor.py
  goto done
)

echo Bundled Windows runtime was not found.
echo This package may be the source-code version, not the Windows portable version.
echo.

set "PYTHON_CMD="
where py >nul 2>nul
if %ERRORLEVEL%==0 set "PYTHON_CMD=py -3"

if not defined PYTHON_CMD (
  where python >nul 2>nul
  if %ERRORLEVEL%==0 set "PYTHON_CMD=python"
)

if not defined PYTHON_CMD (
  echo Python was not found.
  echo Ask for the Windows portable zip: Locus-review-tool-Windows-portable.zip
  echo That version opens directly without installing Python.
  echo.
  pause
  exit /b 1
)

if not exist ".venv\Scripts\python.exe" (
  echo Creating local Python environment...
  %PYTHON_CMD% -m venv .venv
  if errorlevel 1 goto failed
)

if not exist ".venv\.locus_review_deps_ok" (
  echo Installing required packages. This may take a few minutes the first time...
  ".venv\Scripts\python.exe" -m pip install --upgrade pip
  if errorlevel 1 goto failed
  ".venv\Scripts\python.exe" -m pip install -r requirements-review.txt
  if errorlevel 1 goto failed
  echo ok > ".venv\.locus_review_deps_ok"
)

".venv\Scripts\python.exe" tools\launch_review_editor.py
goto done

:failed
echo.
echo Locus failed to start. Please keep this window open and check the message above.
echo For non-technical Windows use, use Locus-review-tool-Windows-portable.zip.
echo.

:done
pause
