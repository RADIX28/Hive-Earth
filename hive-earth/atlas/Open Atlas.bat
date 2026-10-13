@echo off
REM Opens the Atlas through a small local web server.
REM The map tiles (OpenStreetMap) block pages opened straight from a file.
cd /d "%~dp0"
set PY=python
if exist "%USERPROFILE%\anaconda3\python.exe" set PY="%USERPROFILE%\anaconda3\python.exe"
start "" http://localhost:8765/index.html
%PY% -m http.server 8765
