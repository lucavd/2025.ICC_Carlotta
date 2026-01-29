#!/bin/bash
# Paper2 Sample Size - callr (Linux version)

cd "$(dirname "$0")"

echo "=== PAPER 2 SAMPLE SIZE ==="
echo "Pool: 110 worker, timeout 1h per scenario"
echo ""

Rscript run_paper2_ss_callr.R

echo ""
echo "=== TERMINATO ==="
