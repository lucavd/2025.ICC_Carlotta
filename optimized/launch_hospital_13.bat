@echo off
echo === Lancio 13 worker HOSPITAL ===
cd /d %~dp0

start "Hospital Worker 1" cmd /k "Rscript run_parallel_hospital.R 1 13"
timeout /t 2 >nul
start "Hospital Worker 2" cmd /k "Rscript run_parallel_hospital.R 2 13"
timeout /t 2 >nul
start "Hospital Worker 3" cmd /k "Rscript run_parallel_hospital.R 3 13"
timeout /t 2 >nul
start "Hospital Worker 4" cmd /k "Rscript run_parallel_hospital.R 4 13"
timeout /t 2 >nul
start "Hospital Worker 5" cmd /k "Rscript run_parallel_hospital.R 5 13"
timeout /t 2 >nul
start "Hospital Worker 6" cmd /k "Rscript run_parallel_hospital.R 6 13"
timeout /t 2 >nul
start "Hospital Worker 7" cmd /k "Rscript run_parallel_hospital.R 7 13"
timeout /t 2 >nul
start "Hospital Worker 8" cmd /k "Rscript run_parallel_hospital.R 8 13"
timeout /t 2 >nul
start "Hospital Worker 9" cmd /k "Rscript run_parallel_hospital.R 9 13"
timeout /t 2 >nul
start "Hospital Worker 10" cmd /k "Rscript run_parallel_hospital.R 10 13"
timeout /t 2 >nul
start "Hospital Worker 11" cmd /k "Rscript run_parallel_hospital.R 11 13"
timeout /t 2 >nul
start "Hospital Worker 12" cmd /k "Rscript run_parallel_hospital.R 12 13"
timeout /t 2 >nul
start "Hospital Worker 13" cmd /k "Rscript run_parallel_hospital.R 13 13"

echo === Tutti i 13 worker HOSPITAL avviati! ===
