@echo off
echo === Lancio 9 worker INDIVIDUAL ===
cd /d %~dp0

start "Individual Worker 1" cmd /k "Rscript run_parallel_individual.R 1 9"
timeout /t 2 >nul
start "Individual Worker 2" cmd /k "Rscript run_parallel_individual.R 2 9"
timeout /t 2 >nul
start "Individual Worker 3" cmd /k "Rscript run_parallel_individual.R 3 9"
timeout /t 2 >nul
start "Individual Worker 4" cmd /k "Rscript run_parallel_individual.R 4 9"
timeout /t 2 >nul
start "Individual Worker 5" cmd /k "Rscript run_parallel_individual.R 5 9"
timeout /t 2 >nul
start "Individual Worker 6" cmd /k "Rscript run_parallel_individual.R 6 9"
timeout /t 2 >nul
start "Individual Worker 7" cmd /k "Rscript run_parallel_individual.R 7 9"
timeout /t 2 >nul
start "Individual Worker 8" cmd /k "Rscript run_parallel_individual.R 8 9"
timeout /t 2 >nul
start "Individual Worker 9" cmd /k "Rscript run_parallel_individual.R 9 9"

echo === Tutti i 9 worker INDIVIDUAL avviati! ===
