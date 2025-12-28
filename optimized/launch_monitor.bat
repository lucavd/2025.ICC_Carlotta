@echo off
title Monitor Progresso ICC
cd /d %~dp0
powershell -ExecutionPolicy Bypass -File "%~dp0monitor_progress.ps1"
pause
