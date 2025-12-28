@echo off
echo === Lancio 13 worker INDIVIDUAL ===
cd /d %~dp0

start "Individual Worker 1" cmd /k "Rscript run_parallel_individual.R 1 13"
timeout /t 2 >nul
start "Individual Worker 2" cmd /k "Rscript run_parallel_individual.R 2 13"
timeout /t 2 >nul
start "Individual Worker 3" cmd /k "Rscript run_parallel_individual.R 3 13"
timeout /t 2 >nul
start "Individual Worker 4" cmd /k "Rscript run_parallel_individual.R 4 13"
timeout /t 2 >nul
start "Individual Worker 5" cmd /k "Rscript run_parallel_individual.R 5 13"
timeout /t 2 >nul
start "Individual Worker 6" cmd /k "Rscript run_parallel_individual.R 6 13"
timeout /t 2 >nul
start "Individual Worker 7" cmd /k "Rscript run_parallel_individual.R 7 13"
timeout /t 2 >nul
start "Individual Worker 8" cmd /k "Rscript run_parallel_individual.R 8 13"
timeout /t 2 >nul
start "Individual Worker 9" cmd /k "Rscript run_parallel_individual.R 9 13"
timeout /t 2 >nul
start "Individual Worker 10" cmd /k "Rscript run_parallel_individual.R 10 13"
timeout /t 2 >nul
start "Individual Worker 11" cmd /k "Rscript run_parallel_individual.R 11 13"
timeout /t 2 >nul
start "Individual Worker 12" cmd /k "Rscript run_parallel_individual.R 12 13"
timeout /t 2 >nul
start "Individual Worker 13" cmd /k "Rscript run_parallel_individual.R 13 13"

echo === Tutti i 13 worker INDIVIDUAL avviati! ===
