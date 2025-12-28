# =============================================================================
# watchdog_hospital.ps1 - Monitora SOLO Hospital e riavvia se bloccato
# =============================================================================

$scriptPath = $PSScriptRoot
if (-not $scriptPath) { $scriptPath = Get-Location }

$hospPath = "$scriptPath\results_optimized\paper1_icc_hospital"

$checkIntervalSec = 300  # 5 minuti
$maxIdleMin = 30         # 30 minuti senza output = bloccato
$cooldownMin = 35        # Dopo riavvio, aspetta 35 min prima di ricontrollare

# Traccia ultimo riavvio per evitare loop
$lastRestartHosp = $null

Write-Host "=== WATCHDOG HOSPITAL ATTIVO ===" -ForegroundColor Green
Write-Host "Check interval: $checkIntervalSec sec"
Write-Host "Max idle time: $maxIdleMin min"
Write-Host "Cooldown dopo riavvio: $cooldownMin min"
Write-Host "Monitoraggio: SOLO hospital"
Write-Host ""

while ($true) {
    $now = Get-Date
    $hospBlocked = $false
    
    # Controlla HOSPITAL
    $lastFileHosp = Get-ChildItem $hospPath -File -ErrorAction SilentlyContinue | 
                    Sort-Object LastWriteTime -Descending | 
                    Select-Object -First 1
    
    if ($lastFileHosp) {
        $idleMinH = [math]::Round(($now - $lastFileHosp.LastWriteTime).TotalMinutes, 0)
        Write-Host "[$($now.ToString('HH:mm:ss'))] Hospital ultimo file: $idleMinH min fa" -ForegroundColor Magenta
        
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
    } else {
        Write-Host "[$($now.ToString('HH:mm:ss'))] Hospital: nessun file trovato" -ForegroundColor DarkGray
    }
    
    # Riavvia HOSPITAL se bloccato
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
    
    # Attendi prima del prossimo check
    Start-Sleep -Seconds $checkIntervalSec
}
