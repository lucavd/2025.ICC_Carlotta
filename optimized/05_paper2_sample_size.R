# --- Per DESIGN INDIVIDUAL (Cerca N pazienti) ---
process_ss_scenario_ind <- function(scenario_id, scen, nsim, max_n, dir_out) {
  filename <- file.path(dir_out, sprintf("scenario_ind_%03d.rds", scenario_id))
  if (file.exists(filename)) return(readRDS(filename))
  
  message(sprintf(">>> Individual SS Scenario %d | ICC=%.2f, Hosp=%d", scenario_id, scen$icc, scen$num_hosp))
  
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
      simula_fun = simulate_survival_cohort_individual
    )
  }, error = function(e) { message("Errore: ", e$message); return(NULL) })
  
  if (!is.null(res)) {
    res$scenario_id <- scenario_id
    saveRDS(res, filename)
  }
  return(res)
}

# --- Per DESIGN HOSPITAL (Cerca N ospedali) ---
process_ss_scenario_hosp <- function(scenario_id, scen, nsim, max_hosp, dir_out) {
  filename <- file.path(dir_out, sprintf("scenario_hosp_%03d.rds", scenario_id))
  if (file.exists(filename)) return(readRDS(filename))
  
  message(sprintf(">>> Hospital SS Scenario %d | ICC=%.2f, Pat/Hosp=%d", scenario_id, scen$icc, scen$n_per_hosp))
  
  res <- tryCatch({
    sample_hosp_size_binary(
      icc = scen$icc,
      num_hosp = scen$num_hosp, 
      pop_treat_effect = scen$pop_treat_effect,
      nsim = nsim,
      min_hosp = 4,
      max_hosp = max_hosp,         
      target_power = 0.8,
      lambda = scen$lambda,
      cens = scen$cens,
      balancing_mode = scen$balancing_mode
    )
  }, error = function(e) { message("Errore: ", e$message); return(NULL) })
  
  if (!is.null(res)) {
    res$scenario_id <- scenario_id
    saveRDS(res, filename)
  }
  return(res)
}


# RUN INDIVIDUAL
run_paper2_ss_individual <- function() {
  cat("\n========== DESIGN INDIVIDUAL (Ricerca Pazienti) ==========\n")
  scenarios <- make_scenarios_individual() # Carica griglia IND
  
  results <- future_map(1:nrow(scenarios), function(i) {
    process_ss_scenario_ind(i, scenarios[i,], PAPER2_SS_NSIM, PAPER2_SS_MAX_PAT, DIR_PAPER2_SS)
  }, .options = furrr_options(seed = TRUE))
  
  bind_rows(results) %>% mutate(design = "individual")
}

# RUN HOSPITAL
run_paper2_ss_hospital <- function() {
  cat("\n========== DESIGN HOSPITAL (Ricerca Ospedali) ==========\n")
  scenarios <- make_scenarios_hospital()   # Carica griglia HOSP
  
  # Directory dedicata
  dir_ss_hosp <- file.path(DIR_BASE, "paper2_sample_size_hospital")
  dir.create(dir_ss_hosp, recursive = TRUE, showWarnings = FALSE)
  
  results <- future_map(1:nrow(scenarios), function(i) {
    process_ss_scenario_hosp(i, scenarios[i,], PAPER2_SS_NSIM, 500, dir_ss_hosp) # max_hosp = 500
  }, .options = furrr_options(seed = TRUE))
  
  bind_rows(results) %>% mutate(design = "hospital")
}



# =============================================================================
# MAIN
# =============================================================================
# Esecuzione (sia da RStudio che da terminale)
ss_individual <- run_paper2_ss_individual()
ss_hospital   <- run_paper2_ss_hospital()

ss_combined <- bind_rows(ss_individual, ss_hospital)
saveRDS(ss_combined, file.path(DIR_BASE, "paper2_SAMPLE_SIZE_ALL_RESULTS.rds"))
cat("\n\nTutti i risultati Paper 2 Sample Size salvati.\n")
