@echo off
title Paper2 Sample Size - callr
cd /d "%~dp0"
echo === PAPER 2 SAMPLE SIZE ===
echo Pool: 10 worker, timeout 1h per scenario
echo.
Rscript run_paper2_ss_callr.R
echo.
echo === TERMINATO ===
pause
