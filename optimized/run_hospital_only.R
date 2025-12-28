# =============================================================================
# run_hospital_only.R - Esegue solo Paper 1 Hospital design
# =============================================================================
# Usa questo script in un secondo terminale mentre individual gira nel primo.
# Le directory di output sono separate, nessun conflitto.
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("     ICC Simulation - Paper 1 HOSPITAL ONLY                     \n")
cat("================================================================\n")
cat("\n")

# Carica configurazione e funzioni
source("01_config.R")
source("02_functions.R")

# Crea directory
create_dirs()

# Funzione per processare scenario (copiata da 03_paper1_icc.R)
process_icc_scenario <- function(scenario_id, scen, n_rep, chunk_size, 
                                  dir_out, simula_fun) {
  
  summary_file <- file.path(dir_out, sprintf("scenario_%03d_summary.rds", scenario_id))
  if (file.exists(summary_file)) {
    message(sprintf("Scenario %d già completato, skip.", scenario_id))
    return(readRDS(summary_file))
  }
  
  message(sprintf(">>> Scenario %d | icc=%.2f, n=%d, hosp=%d, beta=%.1f",
                  scenario_id, scen$icc, scen$sample_size, scen$num_hosp, scen$beta))
  
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
  
  saveRDS(summary_icc, summary_file)
  message(sprintf("Scenario %d salvato.", scenario_id))
  
  return(summary_icc)
}

start_time <- Sys.time()

cat("\n========== PAPER 1 - DESIGN HOSPITAL ==========\n")

scenarios <- make_paper1_scenarios()
cat(sprintf("Scenari da processare: %d\n", nrow(scenarios)))
cat(sprintf("Repliche per scenario: %d\n", PAPER1_N_REP))
cat(sprintf("Chunk size: %d\n", PAPER1_CHUNK_SIZE))
cat(sprintf("Totale simulazioni: %s\n\n", 
            format(nrow(scenarios) * PAPER1_N_REP, big.mark = ",")))

# Loop sequenziale (mostra output)
results <- vector("list", nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  cat(sprintf("\n[%d/%d] ", i, nrow(scenarios)))
  results[[i]] <- process_icc_scenario(
    scenario_id = i,
    scen = scenarios[i, ],
    n_rep = PAPER1_N_REP,
    chunk_size = PAPER1_CHUNK_SIZE,
    dir_out = DIR_PAPER1_HOSP,
    simula_fun = simulate_survival_cohort_hospital
  )
}

final_results <- bind_rows(results)

final_file <- file.path(DIR_PAPER1_HOSP, "ICC_results_hospital_FINAL.rds")
saveRDS(final_results, final_file)
cat(sprintf("\nRisultati finali salvati in: %s\n", final_file))

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "hours")

cat("\n")
cat("================================================================\n")
cat(sprintf("  HOSPITAL COMPLETATO in %.2f ore\n", as.numeric(elapsed)))
cat("================================================================\n")
