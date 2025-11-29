# =============================================================================
# 04_paper2_power.R - Paper 2: ICC Impact on Statistical Power
# =============================================================================
# Valuta come ICC influenza la potenza statistica e calcola Design Effect (DE)
#
# Workflow:
#   1. Per ogni scenario, stima power via Monte Carlo
#   2. Per scenari con ICC > 0, calcola DE e sample_size corretto
#   3. Ri-esegui con sample_size_DE per verificare che power ≈ power(ICC=0)
#
# Output: tabella con power, DE, power_DE, confronto formule Schoenfeld/Xie
# =============================================================================

source("01_config.R")
source("02_functions.R")

create_dirs()
setup_parallel(N_CORES)

# =============================================================================
# FUNZIONE PER PROCESSARE UN SINGOLO SCENARIO POWER
# =============================================================================

process_power_scenario <- function(scenario_id, scen, nsim, dir_out, simula_fun) {
  
  filename <- file.path(dir_out, sprintf("scenario_%03d.rds", scenario_id))
  
  # Check cache
  if (file.exists(filename)) {
    return(readRDS(filename))
  }
  
  message(sprintf(">>> Power scenario %d | icc=%.2f, n=%d, hosp=%d, effect=%.1f",
                  scenario_id, scen$icc, scen$sample_size, scen$num_hosp, 
                  scen$pop_treat_effect))
  
  res <- tryCatch({
    surv_power_function(
      simula_coorte_fun = simula_fun,
      simula_args = list(
        num_hosp = scen$num_hosp,
        sample_size = scen$sample_size,
        balancing_mode = scen$balancing_mode,
        pop_treat_effect = scen$pop_treat_effect,
        lambda = scen$lambda,
        icc = scen$icc,
        cens = scen$cens
      ),
      nsim = nsim
    )
  }, error = function(e) NULL)
  
  if (is.null(res)) {
    res <- tibble(
      nsim = nsim,
      icc = scen$icc,
      num_hosp = scen$num_hosp,
      num_pat_group_mean = scen$sample_size / scen$num_hosp,
      sample_size = scen$sample_size,
      power = NA_real_,
      prop_cens = NA_real_
    )
  }
  
  # Aggiungi info scenario
  res <- res %>%
    mutate(
      scenario_id = scenario_id,
      pop_treat_effect = scen$pop_treat_effect,
      lambda = scen$lambda,
      cens = scen$cens,
      balancing_mode = scen$balancing_mode,
      # Design Effect
      DE = 1 + (sample_size / num_hosp - 1) * icc,
      sample_size_DE = ceiling(sample_size * DE)
    )
  
  saveRDS(res, filename)
  return(res)
}


# =============================================================================
# FUNZIONE PER SCENARIO DE (verifica power con sample size corretto)
# =============================================================================

process_power_DE_scenario <- function(scenario_id, scen, nsim, dir_out, simula_fun) {
  
  filename <- file.path(dir_out, sprintf("scenario_DE_%03d.rds", scenario_id))
  
  if (file.exists(filename)) {
    return(readRDS(filename))
  }
  
  # Skip se ICC = 0 (DE = 1, niente da verificare)
  if (scen$icc == 0) {
    return(NULL)
  }
  
  message(sprintf(">>> Power DE scenario %d | icc=%.2f, n_DE=%d", 
                  scenario_id, scen$icc, scen$sample_size_DE))
  
  res <- tryCatch({
    surv_power_function(
      simula_coorte_fun = simula_fun,
      simula_args = list(
        num_hosp = scen$num_hosp,
        sample_size = scen$sample_size_DE,  # <- sample size corretto per DE
        balancing_mode = scen$balancing_mode,
        pop_treat_effect = scen$pop_treat_effect,
        lambda = scen$lambda,
        icc = scen$icc,
        cens = scen$cens
      ),
      nsim = nsim
    )
  }, error = function(e) NULL)
  
  if (!is.null(res)) {
    res <- res %>%
      mutate(
        scenario_id = scenario_id,
        sample_size_original = scen$sample_size
      ) %>%
      rename(power_DE = power, prop_cens_DE = prop_cens)
    
    saveRDS(res, filename)
  }
  
  return(res)
}


# =============================================================================
# RUN PAPER 2 POWER - DESIGN INDIVIDUAL
# =============================================================================

run_paper2_power_individual <- function() {
  
  cat("\n========== PAPER 2 POWER - DESIGN INDIVIDUAL ==========\n")
  
  scenarios <- make_paper2_power_scenarios()
  cat(sprintf("Scenari base: %d\n", nrow(scenarios)))
  cat(sprintf("nsim per scenario: %d\n", PAPER2_POWER_NSIM))
  
  # --- STEP 1: Power base ---
  cat("\n--- Step 1: Calcolo power base ---\n")
  
  with_progress({
    p <- progressor(steps = nrow(scenarios))
    
    results_base <- future_map(1:nrow(scenarios), function(i) {
      res <- process_power_scenario(
        scenario_id = i,
        scen = scenarios[i, ],
        nsim = PAPER2_POWER_NSIM,
        dir_out = DIR_PAPER2_POWER_IND,
        simula_fun = simulate_survival_cohort_individual
      )
      p()
      res
    }, .options = furrr_options(seed = TRUE))
  })
  
  results_df <- bind_rows(results_base)
  
  # --- STEP 2: Power con DE (solo ICC > 0) ---
  cat("\n--- Step 2: Calcolo power con sample size DE ---\n")
  
  # Prepara scenari DE
  scenarios_DE <- results_df %>%
    filter(icc > 0, !is.na(sample_size_DE))
  
  cat(sprintf("Scenari DE da processare: %d\n", nrow(scenarios_DE)))
  
  with_progress({
    p <- progressor(steps = nrow(scenarios_DE))
    
    results_DE <- future_map(1:nrow(scenarios_DE), function(i) {
      res <- process_power_DE_scenario(
        scenario_id = scenarios_DE$scenario_id[i],
        scen = scenarios_DE[i, ],
        nsim = PAPER2_POWER_NSIM,
        dir_out = DIR_PAPER2_POWER_IND,
        simula_fun = simulate_survival_cohort_individual
      )
      p()
      res
    }, .options = furrr_options(seed = TRUE))
  })
  
  results_DE_df <- bind_rows(results_DE)
  
  # --- STEP 3: Combina e aggiungi formule teoriche ---
  cat("\n--- Step 3: Combinazione risultati ---\n")
  
  if (nrow(results_DE_df) > 0) {
    results_final <- results_df %>%
      left_join(
        results_DE_df %>% select(scenario_id, power_DE, prop_cens_DE),
        by = "scenario_id"
      )
  } else {
    results_final <- results_df %>%
      mutate(power_DE = NA_real_, prop_cens_DE = NA_real_)
  }
  
  # Aggiungi calcoli Schoenfeld/Freedman e Xie
  results_final <- results_final %>%
    mutate(
      alpha = 0.05,
      z_alpha = qnorm(1 - alpha / 2),
      p = 0.5,
      HR = exp(pop_treat_effect),
      z_beta = qnorm(pmax(power, 0.01, na.rm = TRUE)),  # evita -Inf
      E_schoenfeld = ((z_alpha + z_beta)^2) / (p * (1 - p) * (log(HR))^2),
      prop_eventi = 1 - prop_cens,
      N_schoenfeld = ceiling(E_schoenfeld / pmax(prop_eventi, 0.01)),
      N_xie = ceiling(N_schoenfeld * DE),
      design = "individual"
    ) %>%
    select(-alpha, -z_alpha, -p, -HR, -z_beta, -E_schoenfeld, -prop_eventi)
  
  # Salva
  final_file <- file.path(DIR_PAPER2_POWER_IND, "power_results_individual_FINAL.rds")
  saveRDS(results_final, final_file)
  cat(sprintf("\nSalvato: %s\n", final_file))
  
  return(results_final)
}


# =============================================================================
# RUN PAPER 2 POWER - DESIGN HOSPITAL
# =============================================================================

run_paper2_power_hospital <- function() {
  
  cat("\n========== PAPER 2 POWER - DESIGN HOSPITAL ==========\n")
  
  scenarios <- make_paper2_power_scenarios()
  cat(sprintf("Scenari base: %d\n", nrow(scenarios)))
  cat(sprintf("nsim per scenario: %d\n", PAPER2_POWER_NSIM))
  
  # --- STEP 1: Power base ---
  cat("\n--- Step 1: Calcolo power base ---\n")
  
  with_progress({
    p <- progressor(steps = nrow(scenarios))
    
    results_base <- future_map(1:nrow(scenarios), function(i) {
      res <- process_power_scenario(
        scenario_id = i,
        scen = scenarios[i, ],
        nsim = PAPER2_POWER_NSIM,
        dir_out = DIR_PAPER2_POWER_HOSP,
        simula_fun = simulate_survival_cohort_hospital
      )
      p()
      res
    }, .options = furrr_options(seed = TRUE))
  })
  
  results_df <- bind_rows(results_base)
  
  # --- STEP 2: Power con DE ---
  cat("\n--- Step 2: Calcolo power con sample size DE ---\n")
  
  scenarios_DE <- results_df %>%
    filter(icc > 0, !is.na(sample_size_DE))
  
  cat(sprintf("Scenari DE da processare: %d\n", nrow(scenarios_DE)))
  
  with_progress({
    p <- progressor(steps = nrow(scenarios_DE))
    
    results_DE <- future_map(1:nrow(scenarios_DE), function(i) {
      res <- process_power_DE_scenario(
        scenario_id = scenarios_DE$scenario_id[i],
        scen = scenarios_DE[i, ],
        nsim = PAPER2_POWER_NSIM,
        dir_out = DIR_PAPER2_POWER_HOSP,
        simula_fun = simulate_survival_cohort_hospital
      )
      p()
      res
    }, .options = furrr_options(seed = TRUE))
  })
  
  results_DE_df <- bind_rows(results_DE)
  
  # --- STEP 3: Combina ---
  cat("\n--- Step 3: Combinazione risultati ---\n")
  
  if (nrow(results_DE_df) > 0) {
    results_final <- results_df %>%
      left_join(
        results_DE_df %>% select(scenario_id, power_DE, prop_cens_DE),
        by = "scenario_id"
      )
  } else {
    results_final <- results_df %>%
      mutate(power_DE = NA_real_, prop_cens_DE = NA_real_)
  }
  
  results_final <- results_final %>%
    mutate(
      alpha = 0.05,
      z_alpha = qnorm(1 - alpha / 2),
      p = 0.5,
      HR = exp(pop_treat_effect),
      z_beta = qnorm(pmax(power, 0.01, na.rm = TRUE)),
      E_schoenfeld = ((z_alpha + z_beta)^2) / (p * (1 - p) * (log(HR))^2),
      prop_eventi = 1 - prop_cens,
      N_schoenfeld = ceiling(E_schoenfeld / pmax(prop_eventi, 0.01)),
      N_xie = ceiling(N_schoenfeld * DE),
      design = "hospital"
    ) %>%
    select(-alpha, -z_alpha, -p, -HR, -z_beta, -E_schoenfeld, -prop_eventi)
  
  final_file <- file.path(DIR_PAPER2_POWER_HOSP, "power_results_hospital_FINAL.rds")
  saveRDS(results_final, final_file)
  cat(sprintf("\nSalvato: %s\n", final_file))
  
  return(results_final)
}


# =============================================================================
# MAIN
# =============================================================================

if (interactive()) {
  cat("\n")
  cat("Per eseguire Paper 2 Power:\n")
  cat("  power_ind  <- run_paper2_power_individual()\n")
  cat("  power_hosp <- run_paper2_power_hospital()\n")
  cat("\n")
} else {
  power_individual <- run_paper2_power_individual()
  power_hospital   <- run_paper2_power_hospital()
  
  power_combined <- bind_rows(power_individual, power_hospital)
  saveRDS(power_combined, file.path(DIR_BASE, "paper2_POWER_ALL_RESULTS.rds"))
  cat("\n\nTutti i risultati Paper 2 Power salvati.\n")
}
