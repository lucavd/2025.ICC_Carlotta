# =============================================================================
# 03_paper1_icc.R - Paper 1: ICC Estimation Methods Comparison
# =============================================================================
# Esegue simulazioni Monte Carlo per confrontare 10 metodi di stima ICC
# su coorti survival clusterizzate (ospedali)
#
# Output: per ogni scenario, riassunto delle stime ICC per metodo
#
# USO:
#   1. Modifica parametri in 01_config.R se necessario
#   2. Esegui questo script
#   3. I risultati sono salvati in results_optimized/paper1_icc_*
# =============================================================================

# --- Setup ---
source("01_config.R")
source("02_functions.R")

# Crea directory
create_dirs()

# Inizializza parallelizzazione
setup_parallel(N_CORES)

# =============================================================================
# FUNZIONE PER PROCESSARE UN SINGOLO SCENARIO
# =============================================================================

process_icc_scenario <- function(scenario_id, scen, n_rep, chunk_size, 
                                  dir_out, simula_fun) {
  
  # Check se già completato
  summary_file <- file.path(dir_out, sprintf("scenario_%03d_summary.rds", scenario_id))
  if (file.exists(summary_file)) {
    message(sprintf("Scenario %d già completato, skip.", scenario_id))
    return(readRDS(summary_file))
  }
  
  message(sprintf(">>> Scenario %d | icc=%.2f, n=%d, hosp=%d, beta=%.1f",
                  scenario_id, scen$icc, scen$sample_size, scen$num_hosp, scen$beta))
  
  # Funzione per singola replica
  single_rep <- function() {
    cohort <- simula_fun(
      num_hosp = scen$num_hosp,
      sample_size = scen$sample_size,
      balancing_mode = scen$balancing_mode,
      icc = scen$icc,
      pop_treat_effect = scen$beta,
      lambda = scen$lambda,
      gammas = 1,
      cens = scen$cens
    )
    
    stima_icc <- icc_estimation(cohort)
    stima_icc$prop_cens <- mean(cohort$status == 0)
    stima_icc
  }
  
  # Esegui con chunking
  rep_results <- run_with_chunks(
    scenario_id = sprintf("%03d", scenario_id),
    n_rep = n_rep,
    chunk_size = chunk_size,
    dir_out = dir_out,
    rep_fun = single_rep
  )
  
  if (nrow(rep_results) == 0) {
    message(sprintf("Scenario %d fallito - nessuna replica valida.", scenario_id))
    return(NULL)
  }
  
  # Riassunto per metodo
  summary_icc <- rep_results %>%
    group_by(Method) %>%
    summarise(
      Mean_ICC = mean(ICC, na.rm = TRUE),
      SD_ICC = sd(ICC, na.rm = TRUE),
      Median_ICC = median(ICC, na.rm = TRUE),
      Lower_95CI = quantile(ICC, probs = 0.025, na.rm = TRUE),
      Upper_95CI = quantile(ICC, probs = 0.975, na.rm = TRUE),
      n_valid = sum(!is.na(ICC)),
      prop_cens = mean(prop_cens, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      scenario_id = scenario_id,
      beta = scen$beta,
      lambda = scen$lambda,
      sample_size = scen$sample_size,
      num_hosp = scen$num_hosp,
      num_pat_group_mean = scen$sample_size / scen$num_hosp,
      icc_input = scen$icc,
      cens = scen$cens,
      balancing_mode = scen$balancing_mode
    )
  
  # Salva summary
  saveRDS(summary_icc, summary_file)
  message(sprintf("Scenario %d salvato.", scenario_id))
  
  return(summary_icc)
}


# =============================================================================
# ESECUZIONE PAPER 1 - DISEGNO INDIVIDUAL
# =============================================================================

run_paper1_individual <- function() {
  
  cat("\n========== PAPER 1 - DESIGN INDIVIDUAL ==========\n")
  
  scenarios <- make_paper1_scenarios()
  cat(sprintf("Scenari da processare: %d\n", nrow(scenarios)))
  cat(sprintf("Repliche per scenario: %d\n", PAPER1_N_REP))
  cat(sprintf("Chunk size: %d\n", PAPER1_CHUNK_SIZE))
  cat(sprintf("Totale simulazioni: %s\n\n", 
              format(nrow(scenarios) * PAPER1_N_REP, big.mark = ",")))
  
  # Esegui in parallelo con progress bar
  with_progress({
    p <- progressor(steps = nrow(scenarios))
    
    results <- future_map(1:nrow(scenarios), function(i) {
      res <- process_icc_scenario(
        scenario_id = i,
        scen = scenarios[i, ],
        n_rep = PAPER1_N_REP,
        chunk_size = PAPER1_CHUNK_SIZE,
        dir_out = DIR_PAPER1_IND,
        simula_fun = simulate_survival_cohort_individual
      )
      p()
      res
    }, .options = furrr_options(seed = TRUE))
  })
  
  # Combina tutti i risultati
  final_results <- bind_rows(results)
  
  # Salva risultato finale
  final_file <- file.path(DIR_PAPER1_IND, "ICC_results_individual_FINAL.rds")
  saveRDS(final_results, final_file)
  cat(sprintf("\nRisultati finali salvati in: %s\n", final_file))
  
  return(final_results)
}


# =============================================================================
# ESECUZIONE PAPER 1 - DISEGNO HOSPITAL
# =============================================================================

run_paper1_hospital <- function() {
  
  cat("\n========== PAPER 1 - DESIGN HOSPITAL ==========\n")
  
  scenarios <- make_paper1_scenarios()
  cat(sprintf("Scenari da processare: %d\n", nrow(scenarios)))
  cat(sprintf("Repliche per scenario: %d\n", PAPER1_N_REP))
  cat(sprintf("Chunk size: %d\n", PAPER1_CHUNK_SIZE))
  cat(sprintf("Totale simulazioni: %s\n\n", 
              format(nrow(scenarios) * PAPER1_N_REP, big.mark = ",")))
  
  with_progress({
    p <- progressor(steps = nrow(scenarios))
    
    results <- future_map(1:nrow(scenarios), function(i) {
      res <- process_icc_scenario(
        scenario_id = i,
        scen = scenarios[i, ],
        n_rep = PAPER1_N_REP,
        chunk_size = PAPER1_CHUNK_SIZE,
        dir_out = DIR_PAPER1_HOSP,
        simula_fun = simulate_survival_cohort_hospital
      )
      p()
      res
    }, .options = furrr_options(seed = TRUE))
  })
  
  final_results <- bind_rows(results)
  
  final_file <- file.path(DIR_PAPER1_HOSP, "ICC_results_hospital_FINAL.rds")
  saveRDS(final_results, final_file)
  cat(sprintf("\nRisultati finali salvati in: %s\n", final_file))
  
  return(final_results)
}


# =============================================================================
# MAIN - Esegui entrambi i disegni
# =============================================================================

if (interactive()) {
  cat("\n")
  cat("Per eseguire Paper 1:\n")
  cat("  results_ind  <- run_paper1_individual()\n")
  cat("  results_hosp <- run_paper1_hospital()\n")
  cat("\n")
  cat("Oppure esegui tutto con:\n")
  cat("  source('03_paper1_icc.R')\n")
  cat("  # poi decommentare le righe sotto\n")
} else {
  # Esecuzione batch (quando si fa source() da terminale)
  results_individual <- run_paper1_individual()
  results_hospital   <- run_paper1_hospital()
  
  # Combina i due disegni
  results_combined <- bind_rows(
    results_individual %>% mutate(design = "individual"),
    results_hospital %>% mutate(design = "hospital")
  )
  
  saveRDS(results_combined, file.path(DIR_BASE, "paper1_ICC_ALL_RESULTS.rds"))
  cat("\n\nTutti i risultati Paper 1 salvati in: results_optimized/paper1_ICC_ALL_RESULTS.rds\n")
}

# Chiudi workers
# plan(sequential)
