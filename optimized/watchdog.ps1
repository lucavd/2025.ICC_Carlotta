# =============================================================================
# watchdog.ps1 - Monitora e riavvia worker bloccati
# =============================================================================
# USO: .\watchdog.ps1
# Controlla ogni 5 minuti se i worker stanno producendo output.
# Se bloccati >20 minuti senza file nuovi → kill e riavvia
# =============================================================================

$scriptPath = $PSScriptRoot
if (-not $scriptPath) { $scriptPath = Get-Location }

$checkIntervalSec = 300  # 5 minuti
$maxIdleMin = 30         # 30 minuti senza output = bloccato
$cooldownMin = 35        # Dopo riavvio, aspetta 35 min prima di ricontrollare

# Traccia ultimo riavvio per evitare loop
$lastRestartInd = $null
$lastRestartHosp = $null

Write-Host "=== WATCHDOG ATTIVO ===" -ForegroundColor Green
Write-Host "Check interval: $checkIntervalSec sec"
Write-Host "Max idle time: $maxIdleMin min"
Write-Host "Cooldown dopo riavvio: $cooldownMin min"
Write-Host "Monitoraggio: individual & hospital"
Write-Host ""

while ($true) {
    $now = Get-Date
    
    $indBlocked = $false
    $hospBlocked = $false
    
    # Controlla individual
    $indFiles = Get-ChildItem "$scriptPath\results_optimized\paper1_icc_individual\*" -File -ErrorAction SilentlyContinue
    if ($indFiles.Count -gt 0) {
        $lastInd = ($indFiles | Sort-Object LastWriteTime -Descending | Select-Object -First 1).LastWriteTime
        $idleMin = [math]::Round(($now - $lastInd).TotalMinutes, 1)
        
        Write-Host "[$($now.ToString('HH:mm:ss'))] Individual: ultimo file $idleMin min fa"
        
        if ($idleMin -gt $maxIdleMin) {
            # Controlla cooldown: non riavviare se riavviato di recente
            if ($lastRestartInd -and ($now - $lastRestartInd).TotalMinutes -lt $cooldownMin) {
                $minLeft = [math]::Round($cooldownMin - ($now - $lastRestartInd).TotalMinutes, 0)
                Write-Host "  Individual idle ma in cooldown (ancora $minLeft min)" -ForegroundColor DarkYellow
            } else {
                $indBlocked = $true
                Write-Host "  ALERT: Individual bloccato!" -ForegroundColor Red
            }
        }
    }
    
    # Controlla hospital
    $hospFiles = Get-ChildItem "$scriptPath\results_optimized\paper1_icc_hospital\*" -File -ErrorAction SilentlyContinue
    if ($hospFiles.Count -gt 0) {
        $lastHosp = ($hospFiles | Sort-Object LastWriteTime -Descending | Select-Object -First 1).LastWriteTime
        $idleMinH = [math]::Round(($now - $lastHosp).TotalMinutes, 1)
        
        Write-Host "[$($now.ToString('HH:mm:ss'))] Hospital: ultimo file $idleMinH min fa"
        
        if ($idleMinH -gt $maxIdleMin) {
            # Controlla cooldown: non riavviare se riavviato di recente
            if ($lastRestartHosp -and ($now - $lastRestartHosp).TotalMinutes -lt $cooldownMin) {
                $minLeft = [math]::Round($cooldownMin - ($now - $lastRestartHosp).TotalMinutes, 0)
                Write-Host "  Hospital idle ma in cooldown (ancora $minLeft min)" -ForegroundColor DarkYellow
            } else {
                $hospBlocked = $true
                Write-Host "  ALERT: Hospital bloccato!" -ForegroundColor Red
            }
        }
    }
    
    # Riavvia INDIVIDUAL se bloccato
    if ($indBlocked) {
        Write-Host "  Riavvio Individual..." -ForegroundColor Yellow
        
        # PRIMA chiudi processi R Individual (usando command line)
        $indRProcs = Get-WmiObject Win32_Process -Filter "Name='Rscript.exe'" -ErrorAction SilentlyContinue | 
            Where-Object { $_.CommandLine -like '*run_parallel_individual*' }
        $countIndR = ($indRProcs | Measure-Object).Count
        if ($countIndR -gt 0) {
            $indRProcs | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
            Write-Host "  Terminati $countIndR processi R Individual" -ForegroundColor Gray
        }
        
        # Chiudi anche le finestre cmd che li hanno lanciati
        $indCmdProcs = Get-WmiObject Win32_Process -Filter "Name='cmd.exe'" -ErrorAction SilentlyContinue | 
            Where-Object { $_.CommandLine -like '*launch_individual*' -or $_.CommandLine -like '*run_parallel_individual*' }
        $countIndCmd = ($indCmdProcs | Measure-Object).Count
        if ($countIndCmd -gt 0) {
            $indCmdProcs | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
            Write-Host "  Chiuse $countIndCmd finestre cmd Individual" -ForegroundColor Gray
        }
        
        # Attendi che i processi si chiudano
        Start-Sleep -Seconds 3
        
        # Riavvia Individual
        Start-Process cmd -ArgumentList "/c", "cd /d `"$scriptPath`" && launch_individual_13.bat" -WindowStyle Minimized
        $lastRestartInd = Get-Date
        Write-Host "  Individual riavviato! Cooldown $cooldownMin min" -ForegroundColor Green
        Start-Sleep -Seconds 5
    }
    
    # Riavvia HOSPITAL se bloccato (separatamente)
    if ($hospBlocked) {
        Write-Host "  Riavvio Hospital..." -ForegroundColor Yellow
        
        # PRIMA chiudi processi R Hospital (usando command line)
        $hospRProcs = Get-WmiObject Win32_Process -Filter "Name='Rscript.exe'" -ErrorAction SilentlyContinue | 
            Where-Object { $_.CommandLine -like '*run_parallel_hospital*' }
        $countHospR = ($hospRProcs | Measure-Object).Count
        if ($countHospR -gt 0) {
            $hospRProcs | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
            Write-Host "  Terminati $countHospR processi R Hospital" -ForegroundColor Gray
        }
        
        # Chiudi anche le finestre cmd che li hanno lanciati
        $hospCmdProcs = Get-WmiObject Win32_Process -Filter "Name='cmd.exe'" -ErrorAction SilentlyContinue | 
            Where-Object { $_.CommandLine -like '*launch_hospital*' -or $_.CommandLine -like '*run_parallel_hospital*' }
        $countHospCmd = ($hospCmdProcs | Measure-Object).Count
        if ($countHospCmd -gt 0) {
            $hospCmdProcs | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
            Write-Host "  Chiuse $countHospCmd finestre cmd Hospital" -ForegroundColor Gray
        }
        
        # Attendi che i processi si chiudano
        Start-Sleep -Seconds 3
        
        # Riavvia Hospital
        Start-Process cmd -ArgumentList "/c", "cd /d `"$scriptPath`" && launch_hospital_26.bat" -WindowStyle Minimized
        $lastRestartHosp = Get-Date
        Write-Host "  Hospital riavviato! Cooldown $cooldownMin min" -ForegroundColor Green
    }
    
    # Conta processi attivi
    $rProcs = (Get-Process -Name R*, Rscript* -ErrorAction SilentlyContinue | 
               Where-Object {$_.Name -match '^R$|^Rscript'}).Count
    Write-Host "[$($now.ToString('HH:mm:ss'))] Processi R attivi: $rProcs"
    Write-Host ""
    
    # Attendi prossimo check
    Start-Sleep -Seconds $checkIntervalSec
}
