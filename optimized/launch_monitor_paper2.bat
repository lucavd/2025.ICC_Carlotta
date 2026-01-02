@echo off
title Monitor Paper2 Power
cd /d %~dp0
powershell -ExecutionPolicy Bypass -File monitor_paper2_power.ps1
pause
