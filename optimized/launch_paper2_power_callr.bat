@echo off
title Paper2 Power Callr - Timeout 60s
cd /d %~dp0
echo === PAPER 2 POWER CALLR - TIMEOUT REALE ===
echo Timeout 60s per scenario, 20 worker paralleli
echo Individual + Hospital
echo.
Rscript run_paper2_power_callr.R
pause
