@echo off
REM Opens Helix Lab (the DNA game) through a small local web server.
REM Opening index.html directly also works; the server just keeps
REM live database lookups and saved progress behaving the same in every browser.
cd /d "%~dp0"
set PY=python
if exist "%USERPROFILE%\anaconda3\python.exe" set PY="%USERPROFILE%\anaconda3\python.exe"
start "" http://localhost:8766/index.html
%PY% -m http.server 8766
