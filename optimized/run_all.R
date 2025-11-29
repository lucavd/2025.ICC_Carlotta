# =============================================================================
# run_all.R - Script master per eseguire tutte le simulazioni
# =============================================================================
# Esegui questo script per lanciare l'intera pipeline.
# Puoi anche eseguire i singoli script separatamente.
#
# USO:
#   Rscript run_all.R          # esegue tutto
#   Rscript run_all.R paper1   # solo Paper 1 (ICC estimation)
#   Rscript run_all.R paper2   # solo Paper 2 (power + sample size)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
run_what <- if (length(args) > 0) args[1] else "all"

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║     ICC Simulation Study - Optimized Pipeline                 ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Carica configurazione
source("01_config.R")

# Mostra sommario riduzione carico
print_reduction_summary()

cat(sprintf("\nWorkers disponibili: %d\n", N_CORES))
cat(sprintf("Run mode: %s\n\n", run_what))

start_time <- Sys.time()

# --- Paper 1 ---
if (run_what %in% c("all", "paper1")) {
  cat("\n========================================\n")
  cat("         PAPER 1 - ICC ESTIMATION       \n")
  cat("========================================\n")
  
  source("03_paper1_icc.R")
}

# --- Paper 2 ---
if (run_what %in% c("all", "paper2")) {
  cat("\n========================================\n")
  cat("         PAPER 2 - POWER                \n")
  cat("========================================\n")
  
  source("04_paper2_power.R")
  
  cat("\n========================================\n")
  cat("         PAPER 2 - SAMPLE SIZE          \n")
  cat("========================================\n")
  
  source("05_paper2_sample_size.R")
}

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "hours")

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat(sprintf("║  COMPLETATO in %.2f ore                                      ║\n", as.numeric(elapsed)))
cat("╚═══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Chiudi workers
plan(sequential)

cat("Output files:\n")
cat("  - results_optimized/paper1_ICC_ALL_RESULTS.rds\n")
cat("  - results_optimized/paper2_POWER_ALL_RESULTS.rds\n")
cat("  - results_optimized/paper2_SAMPLE_SIZE_ALL_RESULTS.rds\n")
