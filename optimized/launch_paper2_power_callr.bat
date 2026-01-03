@echo off
title Paper2 Power Callr - Timeout 180s
cd /d %~dp0
echo === PAPER 2 POWER CALLR - TIMEOUT REALE ===
echo Timeout 180s per scenario, 10 worker paralleli
echo Individual + Hospital
echo.
Rscript run_paper2_power_callr.R
pause
