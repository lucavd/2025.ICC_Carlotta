# =============================================================================
# 05_paper2_sample_size.R - Paper 2: Sample Size for Target Power
# =============================================================================
# Trova la sample size minima per ottenere power >= 0.8
# Usa binary search invece di incremento +1 (molto più veloce)
#
# Output: per ogni scenario, n minimo per power 80%
# =============================================================================

source("01_config.R")
source("02_functions.R")

create_dirs()
setup_parallel(N_CORES)

# =============================================================================
# FUNZIONE PER PROCESSARE UN SINGOLO SCENARIO SAMPLE SIZE
# =============================================================================

process_ss_scenario <- function(scenario_id, scen, nsim, max_n, 
                                 dir_out, simula_fun) {
  
  filename <- file.path(dir_out, sprintf("scenario_%03d.rds", scenario_id))
  
  if (file.exists(filename)) {
    return(readRDS(filename))
  }
  
  message(sprintf(">>> Sample size scenario %d | icc=%.2f, hosp=%d, effect=%.1f",
                  scenario_id, scen$icc, scen$num_hosp, scen$pop_treat_effect))
  
  res <- tryCatch({
    sample_size_icc_binary(
      icc = scen$icc,
      num_hosp = scen$num_hosp,
      pop_treat_effect = scen$pop_treat_effect,
      nsim = nsim,
      min_n = 2,
      max_n = max_n,
      target_power = 0.8,
      lambda = scen$lambda,
      cens = scen$cens,
      balancing_mode = scen$balancing_mode,
      simula_fun = simula_fun
    )
  }, error = function(e) {
    message(sprintf("Errore scenario %d: %s", scenario_id, e$message))
    NULL
  })
  
  if (is.null(res)) {
    res <- tibble(
      nsim = nsim,
      icc = scen$icc,
      num_hosp = scen$num_hosp,
      lambda = scen$lambda,
      cens = scen$cens,
      pop_treat_effect = scen$pop_treat_effect,
      balancing_mode = scen$balancing_mode,
      num_pat_group = NA_integer_,
      sample_size = NA_integer_,
      power = NA_real_
    )
  }
  
  res <- res %>%
    mutate(scenario_id = scenario_id)
  
  saveRDS(res, filename)
  return(res)
}


# =============================================================================
# RUN PAPER 2 SAMPLE SIZE - DESIGN INDIVIDUAL
# =============================================================================

run_paper2_ss_individual <- function() {
  
  cat("\n========== PAPER 2 SAMPLE SIZE - DESIGN INDIVIDUAL ==========\n")
  
  scenarios <- make_paper2_ss_scenarios()
  cat(sprintf("Scenari: %d\n", nrow(scenarios)))
  cat(sprintf("nsim per step: %d\n", PAPER2_SS_NSIM))
  cat(sprintf("Max n per gruppo: %d\n", PAPER2_SS_MAX_PAT))
  
  with_progress({
    p <- progressor(steps = nrow(scenarios))
    
    results <- future_map(1:nrow(scenarios), function(i) {
      res <- process_ss_scenario(
        scenario_id = i,
        scen = scenarios[i, ],
        nsim = PAPER2_SS_NSIM,
        max_n = PAPER2_SS_MAX_PAT,
        dir_out = DIR_PAPER2_SS,
        simula_fun = simulate_survival_cohort_individual
      )
      p()
      res
    }, .options = furrr_options(seed = TRUE))
  })
  
  results_df <- bind_rows(results) %>%
    mutate(design = "individual")
  
  final_file <- file.path(DIR_PAPER2_SS, "sample_size_individual_FINAL.rds")
  saveRDS(results_df, final_file)
  cat(sprintf("\nSalvato: %s\n", final_file))
  
  return(results_df)
}


# =============================================================================
# RUN PAPER 2 SAMPLE SIZE - DESIGN HOSPITAL
# =============================================================================

run_paper2_ss_hospital <- function() {
  
  cat("\n========== PAPER 2 SAMPLE SIZE - DESIGN HOSPITAL ==========\n")
  
  scenarios <- make_paper2_ss_scenarios()
  cat(sprintf("Scenari: %d\n", nrow(scenarios)))
  cat(sprintf("nsim per step: %d\n", PAPER2_SS_NSIM))
  cat(sprintf("Max n per gruppo: %d\n", PAPER2_SS_MAX_PAT))
  
  # Directory separata per hospital
  dir_ss_hosp <- file.path(DIR_BASE, "paper2_sample_size_hospital")
  dir.create(dir_ss_hosp, recursive = TRUE, showWarnings = FALSE)
  
  with_progress({
    p <- progressor(steps = nrow(scenarios))
    
    results <- future_map(1:nrow(scenarios), function(i) {
      res <- process_ss_scenario(
        scenario_id = i,
        scen = scenarios[i, ],
        nsim = PAPER2_SS_NSIM,
        max_n = PAPER2_SS_MAX_PAT,
        dir_out = dir_ss_hosp,
        simula_fun = simulate_survival_cohort_hospital
      )
      p()
      res
    }, .options = furrr_options(seed = TRUE))
  })
  
  results_df <- bind_rows(results) %>%
    mutate(design = "hospital")
  
  final_file <- file.path(dir_ss_hosp, "sample_size_hospital_FINAL.rds")
  saveRDS(results_df, final_file)
  cat(sprintf("\nSalvato: %s\n", final_file))
  
  return(results_df)
}


# =============================================================================
# MAIN
# =============================================================================

if (interactive()) {
  cat("\n")
  cat("Per eseguire Paper 2 Sample Size:\n")
  cat("  ss_ind  <- run_paper2_ss_individual()\n")
  cat("  ss_hosp <- run_paper2_ss_hospital()\n")
  cat("\n")
} else {
  ss_individual <- run_paper2_ss_individual()
  ss_hospital   <- run_paper2_ss_hospital()
  
  ss_combined <- bind_rows(ss_individual, ss_hospital)
  saveRDS(ss_combined, file.path(DIR_BASE, "paper2_SAMPLE_SIZE_ALL_RESULTS.rds"))
  cat("\n\nTutti i risultati Paper 2 Sample Size salvati.\n")
}
