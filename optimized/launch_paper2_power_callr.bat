@echo off
title Paper2 Power Callr - Pool Dinamico
cd /d %~dp0
echo === PAPER 2 POWER CALLR - SEQUENZIALE ===
echo 1 scenario alla volta, timeout 15 min
echo Individual + Hospital
echo.
Rscript run_paper2_power_callr.R
pause
