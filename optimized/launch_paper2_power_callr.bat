@echo off
title Paper2 Power Callr - Pool Dinamico
cd /d %~dp0
echo === PAPER 2 POWER CALLR - POOL DINAMICO ===
echo Max 20 worker, timeout 10 min per processo
echo Individual + Hospital
echo.
Rscript run_paper2_power_callr.R
pause
