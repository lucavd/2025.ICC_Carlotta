# =============================================================================
# launch_parallel.ps1 - Lancia N worker in parallelo
# =============================================================================
# USO: .\launch_parallel.ps1 -Design individual -Workers 4
#      .\launch_parallel.ps1 -Design hospital -Workers 4
# =============================================================================

param(
    [Parameter(Mandatory=$true)]
    [ValidateSet("individual", "hospital")]
    [string]$Design,
    
    [Parameter(Mandatory=$true)]
    [int]$Workers
)

$scriptPath = $PSScriptRoot
if (-not $scriptPath) { $scriptPath = Get-Location }

Write-Host "=== Lancio $Workers worker per $Design ===" -ForegroundColor Green

for ($i = 1; $i -le $Workers; $i++) {
    $script = "run_parallel_$Design.R"
    $cmd = "cd '$scriptPath'; Rscript $script $i $Workers"
    
    Write-Host "Avvio Worker $i/$Workers..." -ForegroundColor Yellow
    Start-Process powershell -ArgumentList "-NoExit", "-Command", $cmd
    
    Start-Sleep -Seconds 2  # Piccola pausa tra avvii
}

Write-Host "`n=== Tutti i worker avviati! ===" -ForegroundColor Green
Write-Host "Controlla le finestre PowerShell aperte."
