# =============================================================================
# monitor_paper2_power.ps1 - Monitora avanzamento Paper 2 Power
# =============================================================================
# 480 scenari base + scenari DE per design (Individual + Hospital)
# =============================================================================

$scriptPath = $PSScriptRoot
if (-not $scriptPath) { $scriptPath = Get-Location }

$totalScenariosBase = 480
$checkIntervalSec = 60  # 1 minuto

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "          MONITOR PAPER 2 POWER                             " -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Obiettivo: $totalScenariosBase scenari base + DE per design" -ForegroundColor Gray
Write-Host "Aggiornamento ogni $($checkIntervalSec) secondi" -ForegroundColor Gray
Write-Host ""

while ($true) {
    $now = Get-Date
    $last1h = $now.AddHours(-1)
    
    # Conta file Individual
    $indBase = (Get-ChildItem "$scriptPath\results_optimized\paper2_power_individual\scenario_[0-9]*.rds" -File -ErrorAction SilentlyContinue).Count
    $indDE = (Get-ChildItem "$scriptPath\results_optimized\paper2_power_individual\scenario_DE_*.rds" -File -ErrorAction SilentlyContinue).Count
    $indTotal = $indBase + $indDE
    $indPercentBase = [math]::Round(($indBase / $totalScenariosBase) * 100, 1)
    
    # File ultima ora (Individual)
    $indAllFiles = Get-ChildItem "$scriptPath\results_optimized\paper2_power_individual\*.rds" -File -ErrorAction SilentlyContinue
    $indLast1h = ($indAllFiles | Where-Object { $_.LastWriteTime -ge $last1h }).Count
    $indRate = $indLast1h
    
    # Conta file Hospital
    $hospBase = (Get-ChildItem "$scriptPath\results_optimized\paper2_power_hospital\scenario_[0-9]*.rds" -File -ErrorAction SilentlyContinue).Count
    $hospDE = (Get-ChildItem "$scriptPath\results_optimized\paper2_power_hospital\scenario_DE_*.rds" -File -ErrorAction SilentlyContinue).Count
    $hospTotal = $hospBase + $hospDE
    $hospPercentBase = [math]::Round(($hospBase / $totalScenariosBase) * 100, 1)
    
    # File ultima ora (Hospital)
    $hospAllFiles = Get-ChildItem "$scriptPath\results_optimized\paper2_power_hospital\*.rds" -File -ErrorAction SilentlyContinue
    $hospLast1h = ($hospAllFiles | Where-Object { $_.LastWriteTime -ge $last1h }).Count
    $hospRate = $hospLast1h
    
    # Velocita combinata
    $scenariosPerHour = $indRate + $hospRate
    
    # Stima tempo rimanente
    $indRemaining = $totalScenariosBase - $indBase
    $hospRemaining = $totalScenariosBase - $hospBase
    
    if ($indRate -gt 0) {
        $indHoursLeft = [math]::Round($indRemaining / $indRate, 1)
    } else { $indHoursLeft = "N/A" }
    
    if ($hospRate -gt 0) {
        $hospHoursLeft = [math]::Round($hospRemaining / $hospRate, 1)
    } else { $hospHoursLeft = "N/A" }
    
    # Barra progresso ASCII
    function Get-ProgressBar {
        param([double]$percent)
        $barLength = 30
        $filled = [math]::Floor(($percent / 100) * $barLength)
        $empty = $barLength - $filled
        $bar = ("#" * $filled) + ("-" * $empty)
        return "[$bar]"
    }
    
    # Output
    Clear-Host
    Write-Host "============================================================" -ForegroundColor Cyan
    Write-Host "          MONITOR PAPER 2 POWER - $(Get-Date -Format 'HH:mm:ss')           " -ForegroundColor Cyan
    Write-Host "============================================================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "[SPEED] Velocita (1h): $scenariosPerHour scenari/ora" -ForegroundColor Yellow
    Write-Host ""
    
    # Individual
    Write-Host "-------------------- INDIVIDUAL ----------------------------" -ForegroundColor Green
    Write-Host "  Base: $indBase / $totalScenariosBase" -ForegroundColor White
    Write-Host "  DE:   $indDE completati" -ForegroundColor White
    $indBar = Get-ProgressBar $indPercentBase
    Write-Host "  $indBar $indPercentBase%" -ForegroundColor Green
    Write-Host "  ETA base: $indHoursLeft ore" -ForegroundColor DarkGreen
    Write-Host ""
    
    # Hospital
    Write-Host "-------------------- HOSPITAL ------------------------------" -ForegroundColor Magenta
    Write-Host "  Base: $hospBase / $totalScenariosBase" -ForegroundColor White
    Write-Host "  DE:   $hospDE completati" -ForegroundColor White
    $hospBar = Get-ProgressBar $hospPercentBase
    Write-Host "  $hospBar $hospPercentBase%" -ForegroundColor Magenta
    Write-Host "  ETA base: $hospHoursLeft ore" -ForegroundColor DarkMagenta
    Write-Host ""
    
    # Totale
    $totalBase = $indBase + $hospBase
    $totalTarget = $totalScenariosBase * 2
    $totalPercent = [math]::Round(($totalBase / $totalTarget) * 100, 1)
    $totalDE = $indDE + $hospDE
    Write-Host "-------------------- TOTALE --------------------------------" -ForegroundColor Cyan
    Write-Host "  Base: $totalBase / $totalTarget | DE: $totalDE" -ForegroundColor White
    $totalBar = Get-ProgressBar $totalPercent
    Write-Host "  $totalBar $totalPercent%" -ForegroundColor Cyan
    Write-Host "============================================================" -ForegroundColor Cyan
    Write-Host ""
    
    Write-Host "Prossimo aggiornamento: $($now.AddSeconds($checkIntervalSec).ToString('HH:mm:ss'))" -ForegroundColor DarkGray
    
    Start-Sleep -Seconds $checkIntervalSec
}
