@echo off
title Hospital Callr - Timeout 30s
cd /d %~dp0
echo === HOSPITAL CALLR - TIMEOUT REALE ===
echo Timeout 30s per replica
echo Repliche lente vengono killate e riprova
echo.
Rscript run_hospital_callr.R
pause
