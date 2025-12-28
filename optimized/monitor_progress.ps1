# =============================================================================
# monitor_progress.ps1 - Monitora avanzamento simulazioni
# =============================================================================
# Totale atteso per design: 432 scenari x (5 chunk + 1 summary) = 2592 file
# =============================================================================

$scriptPath = $PSScriptRoot
if (-not $scriptPath) { $scriptPath = Get-Location }

$totalScenarios = 432
$chunksPerScenario = 5
$totalChunks = $totalScenarios * $chunksPerScenario
$totalSummaries = $totalScenarios
$totalFilesPerDesign = $totalChunks + $totalSummaries

$checkIntervalSec = 300  # 5 minuti

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "          MONITOR AVANZAMENTO SIMULAZIONI ICC               " -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Obiettivo per design: $totalScenarios scenari x 6 file = $totalFilesPerDesign file" -ForegroundColor Gray
Write-Host "Aggiornamento ogni $($checkIntervalSec/60) minuti" -ForegroundColor Gray
Write-Host ""

while ($true) {
    $now = Get-Date
    $last24h = $now.AddHours(-24)
    
    # Conta file Individual
    $indAllSummaries = Get-ChildItem "$scriptPath\results_optimized\paper1_icc_individual\*summary*.rds" -File -ErrorAction SilentlyContinue
    $indChunks = (Get-ChildItem "$scriptPath\results_optimized\paper1_icc_individual\*chunk*.rds" -File -ErrorAction SilentlyContinue).Count
    $indSummaries = $indAllSummaries.Count
    $indTotal = $indChunks + $indSummaries
    $indPercent = [math]::Round(($indTotal / $totalFilesPerDesign) * 100, 1)
    $indScenComplete = $indSummaries
    
    # Scenari completati ultime 24h (Individual)
    $indLast24h = ($indAllSummaries | Where-Object { $_.LastWriteTime -ge $last24h }).Count
    $indRate = [math]::Round($indLast24h / 24, 2)
    
    # Conta file Hospital
    $hospAllSummaries = Get-ChildItem "$scriptPath\results_optimized\paper1_icc_hospital\*summary*.rds" -File -ErrorAction SilentlyContinue
    $hospChunks = (Get-ChildItem "$scriptPath\results_optimized\paper1_icc_hospital\*chunk*.rds" -File -ErrorAction SilentlyContinue).Count
    $hospSummaries = $hospAllSummaries.Count
    $hospTotal = $hospChunks + $hospSummaries
    $hospPercent = [math]::Round(($hospTotal / $totalFilesPerDesign) * 100, 1)
    $hospScenComplete = $hospSummaries
    
    # Scenari completati ultime 24h (Hospital)
    $hospLast24h = ($hospAllSummaries | Where-Object { $_.LastWriteTime -ge $last24h }).Count
    $hospRate = [math]::Round($hospLast24h / 24, 2)
    
    # Velocita combinata (ultime 24h)
    $scenariosPerHour = [math]::Round(($indLast24h + $hospLast24h) / 24, 2)
    
    # Stima tempo rimanente per ognuno
    $indRemaining = $totalScenarios - $indScenComplete
    $hospRemaining = $totalScenarios - $hospScenComplete
    
    if ($indRate -gt 0) {
        $indHoursLeft = [math]::Round($indRemaining / $indRate, 1)
        $indDaysLeft = [math]::Round($indHoursLeft / 24, 1)
    } else { $indHoursLeft = 0; $indDaysLeft = 0 }
    
    if ($hospRate -gt 0) {
        $hospHoursLeft = [math]::Round($hospRemaining / $hospRate, 1)
        $hospDaysLeft = [math]::Round($hospHoursLeft / 24, 1)
    } else { $hospHoursLeft = 0; $hospDaysLeft = 0 }
    
    # Stima globale (prende il max)
    $daysRemaining = [math]::Max($indDaysLeft, $hospDaysLeft)
    $hoursRemaining = [math]::Max($indHoursLeft, $hospHoursLeft)
    
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
    Write-Host "          MONITOR AVANZAMENTO SIMULAZIONI ICC               " -ForegroundColor Cyan
    Write-Host "============================================================" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "[SPEED] Velocita (24h): $scenariosPerHour scenari/ora (Ind: $indRate, Hosp: $hospRate)" -ForegroundColor Yellow
    if ($daysRemaining -gt 0) {
        Write-Host "[STIMA] Fine tutto: ~$daysRemaining giorni" -ForegroundColor Yellow
        Write-Host "        Individual: ~$indDaysLeft giorni | Hospital: ~$hospDaysLeft giorni" -ForegroundColor DarkYellow
    }
    Write-Host ""
    
    # Individual
    Write-Host "-------------------- INDIVIDUAL ----------------------------" -ForegroundColor Green
    Write-Host "  Scenari: $indScenComplete / $totalScenarios completati" -ForegroundColor White
    Write-Host "  Files:   $indTotal / $totalFilesPerDesign ($indChunks chunk + $indSummaries summary)" -ForegroundColor White
    $indBar = Get-ProgressBar $indPercent
    Write-Host "  $indBar $indPercent pct" -ForegroundColor Green
    Write-Host ""
    
    # Hospital
    Write-Host "-------------------- HOSPITAL ------------------------------" -ForegroundColor Magenta
    Write-Host "  Scenari: $hospScenComplete / $totalScenarios completati" -ForegroundColor White
    Write-Host "  Files:   $hospTotal / $totalFilesPerDesign ($hospChunks chunk + $hospSummaries summary)" -ForegroundColor White
    $hospBar = Get-ProgressBar $hospPercent
    Write-Host "  $hospBar $hospPercent pct" -ForegroundColor Magenta
    Write-Host ""
    
    # Totale
    $totalComplete = $indScenComplete + $hospScenComplete
    $totalTarget = $totalScenarios * 2
    $totalPercent = [math]::Round(($totalComplete / $totalTarget) * 100, 1)
    Write-Host "-------------------- TOTALE --------------------------------" -ForegroundColor Cyan
    Write-Host "  Scenari: $totalComplete / $totalTarget completati" -ForegroundColor White
    $totalBar = Get-ProgressBar $totalPercent
    Write-Host "  $totalBar $totalPercent pct" -ForegroundColor Cyan
    Write-Host "============================================================" -ForegroundColor Cyan
    Write-Host ""
    
    Write-Host "Prossimo aggiornamento: $($now.AddSeconds($checkIntervalSec).ToString('HH:mm:ss'))" -ForegroundColor DarkGray
    Write-Host ""
    
    # Attendi 5 minuti
    Start-Sleep -Seconds $checkIntervalSec
}
