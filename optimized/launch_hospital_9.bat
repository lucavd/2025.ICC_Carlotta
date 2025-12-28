@echo off
echo === Lancio 9 worker HOSPITAL ===
cd /d %~dp0

start "Hospital Worker 1" cmd /k "Rscript run_parallel_hospital.R 1 9"
timeout /t 2 >nul
start "Hospital Worker 2" cmd /k "Rscript run_parallel_hospital.R 2 9"
timeout /t 2 >nul
start "Hospital Worker 3" cmd /k "Rscript run_parallel_hospital.R 3 9"
timeout /t 2 >nul
start "Hospital Worker 4" cmd /k "Rscript run_parallel_hospital.R 4 9"
timeout /t 2 >nul
start "Hospital Worker 5" cmd /k "Rscript run_parallel_hospital.R 5 9"
timeout /t 2 >nul
start "Hospital Worker 6" cmd /k "Rscript run_parallel_hospital.R 6 9"
timeout /t 2 >nul
start "Hospital Worker 7" cmd /k "Rscript run_parallel_hospital.R 7 9"
timeout /t 2 >nul
start "Hospital Worker 8" cmd /k "Rscript run_parallel_hospital.R 8 9"
timeout /t 2 >nul
start "Hospital Worker 9" cmd /k "Rscript run_parallel_hospital.R 9 9"

echo === Tutti i 9 worker HOSPITAL avviati! ===
