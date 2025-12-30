@echo off
title Individual Callr - Timeout 30s
cd /d %~dp0
echo === INDIVIDUAL CALLR - TIMEOUT REALE ===
echo Timeout 30s per replica
echo Repliche lente vengono killate e riprova
echo.
Rscript run_individual_callr.R
pause
