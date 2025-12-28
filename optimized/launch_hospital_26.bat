@echo off
echo === Lancio 26 worker HOSPITAL ===
cd /d %~dp0

start "Hospital Worker 1" cmd /k "Rscript run_parallel_hospital.R 1 26"
timeout /t 2 >nul
start "Hospital Worker 2" cmd /k "Rscript run_parallel_hospital.R 2 26"
timeout /t 2 >nul
start "Hospital Worker 3" cmd /k "Rscript run_parallel_hospital.R 3 26"
timeout /t 2 >nul
start "Hospital Worker 4" cmd /k "Rscript run_parallel_hospital.R 4 26"
timeout /t 2 >nul
start "Hospital Worker 5" cmd /k "Rscript run_parallel_hospital.R 5 26"
timeout /t 2 >nul
start "Hospital Worker 6" cmd /k "Rscript run_parallel_hospital.R 6 26"
timeout /t 2 >nul
start "Hospital Worker 7" cmd /k "Rscript run_parallel_hospital.R 7 26"
timeout /t 2 >nul
start "Hospital Worker 8" cmd /k "Rscript run_parallel_hospital.R 8 26"
timeout /t 2 >nul
start "Hospital Worker 9" cmd /k "Rscript run_parallel_hospital.R 9 26"
timeout /t 2 >nul
start "Hospital Worker 10" cmd /k "Rscript run_parallel_hospital.R 10 26"
timeout /t 2 >nul
start "Hospital Worker 11" cmd /k "Rscript run_parallel_hospital.R 11 26"
timeout /t 2 >nul
start "Hospital Worker 12" cmd /k "Rscript run_parallel_hospital.R 12 26"
timeout /t 2 >nul
start "Hospital Worker 13" cmd /k "Rscript run_parallel_hospital.R 13 26"
timeout /t 2 >nul
start "Hospital Worker 14" cmd /k "Rscript run_parallel_hospital.R 14 26"
timeout /t 2 >nul
start "Hospital Worker 15" cmd /k "Rscript run_parallel_hospital.R 15 26"
timeout /t 2 >nul
start "Hospital Worker 16" cmd /k "Rscript run_parallel_hospital.R 16 26"
timeout /t 2 >nul
start "Hospital Worker 17" cmd /k "Rscript run_parallel_hospital.R 17 26"
timeout /t 2 >nul
start "Hospital Worker 18" cmd /k "Rscript run_parallel_hospital.R 18 26"
timeout /t 2 >nul
start "Hospital Worker 19" cmd /k "Rscript run_parallel_hospital.R 19 26"
timeout /t 2 >nul
start "Hospital Worker 20" cmd /k "Rscript run_parallel_hospital.R 20 26"
timeout /t 2 >nul
start "Hospital Worker 21" cmd /k "Rscript run_parallel_hospital.R 21 26"
timeout /t 2 >nul
start "Hospital Worker 22" cmd /k "Rscript run_parallel_hospital.R 22 26"
timeout /t 2 >nul
start "Hospital Worker 23" cmd /k "Rscript run_parallel_hospital.R 23 26"
timeout /t 2 >nul
start "Hospital Worker 24" cmd /k "Rscript run_parallel_hospital.R 24 26"
timeout /t 2 >nul
start "Hospital Worker 25" cmd /k "Rscript run_parallel_hospital.R 25 26"
timeout /t 2 >nul
start "Hospital Worker 26" cmd /k "Rscript run_parallel_hospital.R 26 26"

echo === Tutti i 26 worker HOSPITAL avviati! ===
