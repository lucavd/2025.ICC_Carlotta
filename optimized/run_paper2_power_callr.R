# =============================================================================
# run_paper2_power_callr.R - Paper 2 Power con timeout REALE via callr
# 20 scenari paralleli, timeout 180s per scenario
# =============================================================================

cat("\n=== PAPER 2 POWER CON TIMEOUT CALLR ===\n")
cat("Timeout per scenario: 180 secondi\n")
cat("20 scenari in parallelo\n\n")

library(callr)
library(dplyr)
source("01_config.R")
source("02_functions.R")

# Parametri
TIMEOUT_SEC <- 180
NSIM <- PAPER2_POWER_NSIM  # 1000
N_PARALLEL <- 20

# Directory output
dir_out_ind <- DIR_PAPER2_POWER_IND
dir_out_hosp <- DIR_PAPER2_POWER_HOSP
dir.create(dir_out_ind, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out_hosp, recursive = TRUE, showWarnings = FALSE)

# Genera scenari
scenarios <- make_paper2_power_scenarios()
n_scenarios <- nrow(scenarios)

cat(sprintf("Totale scenari: %d\n", n_scenarios))
cat(sprintf("nsim per scenario: %d\n", NSIM))
cat(sprintf("Timeout per scenario: %d secondi\n", TIMEOUT_SEC))
cat(sprintf("Scenari paralleli: %d\n\n", N_PARALLEL))

# =============================================================================
# FUNZIONE PER ESEGUIRE UN SINGOLO SCENARIO IN PROCESSO SEPARATO
# =============================================================================

# Lancia processo in background (NON aspetta)
launch_scenario <- function(scenario_id, scen, nsim, dir_out, simula_fun_name) {
  
  filename <- file.path(dir_out, sprintf("scenario_%03d.rds", scenario_id))
  if (file.exists(filename)) return(NULL)  # già completato
  
  r_bg(
    function(scenario_id, scen, nsim, dir_out, simula_fun_name) {
      suppressPackageStartupMessages({
        source("01_config.R")
        source("02_functions.R")
      })
      
      simula_fun <- get(simula_fun_name)
      
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
          cv = NA_real_,
          prop_cens = NA_real_
        )
      }
      
      res <- res %>%
        mutate(
          scenario_id = scenario_id,
          pop_treat_effect = scen$pop_treat_effect,
          lambda = scen$lambda,
          cens = scen$cens,
          balancing_mode = scen$balancing_mode,
          DE = 1 + (sample_size / num_hosp - 1) * icc,
          DE_CV = 1 + ((cv^2 + 1) * (sample_size / num_hosp) - 1) * icc,
          sample_size_DE = ceiling(sample_size * DE),
          sample_size_DE_CV = ceiling(sample_size * DE_CV)
        )
      
      filename <- file.path(dir_out, sprintf("scenario_%03d.rds", scenario_id))
      saveRDS(res, filename)
      
      return(res)
    },
    args = list(
      scenario_id = scenario_id, 
      scen = scen, 
      nsim = nsim, 
      dir_out = dir_out, 
      simula_fun_name = simula_fun_name
    ),
    package = TRUE
  )
}

# Aspetta lista di processi con timeout globale
wait_processes <- function(procs, timeout_sec) {
  start <- Sys.time()
  
  repeat {
    # Conta quanti ancora vivi
    alive <- sum(sapply(procs, function(p) !is.null(p) && p$is_alive()))
    
    if (alive == 0) break
    
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    if (elapsed > timeout_sec) {
      # Kill tutti quelli ancora vivi
      for (p in procs) {
        if (!is.null(p) && p$is_alive()) p$kill()
      }
      break
    }
    
    Sys.sleep(0.5)
  }
  
  # Conta successi/fallimenti
  ok <- 0
  fail <- 0
  for (p in procs) {
    if (is.null(p)) {
      ok <- ok + 1  # era cached
    } else {
      res <- tryCatch(p$get_result(), error = function(e) NULL)
      if (!is.null(res)) ok <- ok + 1 else fail <- fail + 1
    }
  }
  
  list(ok = ok, fail = fail)
}

# =============================================================================
# FUNZIONE PER ESEGUIRE SCENARIO DE
# =============================================================================

# Lancia processo DE in background (NON aspetta)
launch_scenario_DE <- function(scenario_id, scen, nsim, dir_out, simula_fun_name) {
  
  filename <- file.path(dir_out, sprintf("scenario_DE_%03d.rds", scenario_id))
  if (file.exists(filename)) return(NULL)  # già completato
  
  # Skip se ICC = 0
  if (scen$icc == 0) return(NULL)
  
  # Determina sample size da usare
  sample_size_used <- if (scen$balancing_mode == 1) {
    scen$sample_size_DE
  } else {
    scen$sample_size_DE_CV
  }
  
  # Skip se troppo grande o NA
  if (is.na(sample_size_used) || sample_size_used > PAPER2_POWER_MAX_SS_DE) {
    return(NULL)
  }
  
  r_bg(
    function(scenario_id, scen, nsim, dir_out, simula_fun_name, sample_size_used) {
      suppressPackageStartupMessages({
        source("01_config.R")
        source("02_functions.R")
      })
      
      simula_fun <- get(simula_fun_name)
      
      res <- tryCatch({
        surv_power_function(
          simula_coorte_fun = simula_fun,
          simula_args = list(
            num_hosp = scen$num_hosp,
            sample_size = sample_size_used,
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
        
        filename <- file.path(dir_out, sprintf("scenario_DE_%03d.rds", scenario_id))
        saveRDS(res, filename)
      }
      
      return(res)
    },
    args = list(
      scenario_id = scenario_id,
      scen = scen,
      nsim = nsim,
      dir_out = dir_out,
      simula_fun_name = simula_fun_name,
      sample_size_used = sample_size_used
    ),
    package = TRUE
  )
}

# =============================================================================
# PROCESSAMENTO BATCH PARALLELO
# =============================================================================

process_scenarios_batch <- function(scenarios, dir_out, simula_fun_name, phase_name) {
  
  cat(sprintf("\n=== %s ===\n", phase_name))
  n_total <- nrow(scenarios)
  
  # Trova scenari già completati
  existing <- list.files(dir_out, pattern = "^scenario_[0-9]+\\.rds$")
  completed_ids <- as.integer(gsub("scenario_([0-9]+)\\.rds", "\\1", existing))
  
  to_process <- setdiff(1:n_total, completed_ids)
  cat(sprintf("Scenari da processare: %d/%d (già completati: %d)\n", 
              length(to_process), n_total, length(completed_ids)))
  
  if (length(to_process) == 0) {
    cat("Tutti gli scenari già completati.\n")
    return(NULL)
  }
  
  start_total <- Sys.time()
  processed <- 0
  timeout_count <- 0
  
  # Processa in batch di N_PARALLEL
  while (length(to_process) > 0) {
    batch_ids <- head(to_process, N_PARALLEL)
    to_process <- setdiff(to_process, batch_ids)
    
    cat(sprintf("\nBatch di %d scenari (rimanenti: %d)... ", length(batch_ids), length(to_process)))
    
    # Lancia TUTTI i processi in parallelo (senza aspettare)
    procs <- lapply(batch_ids, function(id) {
      launch_scenario(
        scenario_id = id,
        scen = scenarios[id, ],
        nsim = NSIM,
        dir_out = dir_out,
        simula_fun_name = simula_fun_name
      )
    })
    
    # Aspetta TUTTI con timeout globale
    results <- wait_processes(procs, TIMEOUT_SEC)
    batch_ok <- results$ok
    batch_timeout <- results$fail
    
    processed <- processed + batch_ok
    timeout_count <- timeout_count + batch_timeout
    
    cat(sprintf("%d OK, %d timeout\n", batch_ok, batch_timeout))
    
    # Progresso
    completed_now <- length(list.files(dir_out, pattern = "^scenario_[0-9]+\\.rds$"))
    elapsed <- as.numeric(difftime(Sys.time(), start_total, units = "mins"))
    if (completed_now > 0) {
      rate <- elapsed / completed_now
      remaining <- n_total - completed_now
      eta <- remaining * rate
      cat(sprintf("  Progresso: %d/%d | Elapsed: %.1f min | ETA: %.1f min\n", 
                  completed_now, n_total, elapsed, eta))
    }
  }
  
  cat(sprintf("\n%s completato: %d processati, %d timeout\n", phase_name, processed, timeout_count))
}

process_scenarios_DE_batch <- function(base_results, dir_out, simula_fun_name, phase_name) {
  
  cat(sprintf("\n=== %s ===\n", phase_name))
  
  # Filtra solo ICC > 0
  scenarios_DE <- base_results %>% filter(icc > 0, !is.na(sample_size_DE))
  n_total <- nrow(scenarios_DE)
  
  cat(sprintf("Scenari DE da processare: %d\n", n_total))
  
  if (n_total == 0) return(NULL)
  
  # Trova già completati
  existing <- list.files(dir_out, pattern = "^scenario_DE_[0-9]+\\.rds$")
  completed_ids <- as.integer(gsub("scenario_DE_([0-9]+)\\.rds", "\\1", existing))
  
  to_process_ids <- setdiff(scenarios_DE$scenario_id, completed_ids)
  cat(sprintf("Da processare: %d (già completati: %d)\n", 
              length(to_process_ids), length(completed_ids)))
  
  if (length(to_process_ids) == 0) {
    cat("Tutti gli scenari DE già completati.\n")
    return(NULL)
  }
  
  start_total <- Sys.time()
  timeout_count <- 0
  
  while (length(to_process_ids) > 0) {
    batch_ids <- head(to_process_ids, N_PARALLEL)
    to_process_ids <- setdiff(to_process_ids, batch_ids)
    
    cat(sprintf("\nBatch DE di %d scenari... ", length(batch_ids)))
    
    # Lancia TUTTI i processi DE in parallelo
    procs <- lapply(batch_ids, function(id) {
      scen_row <- scenarios_DE %>% filter(scenario_id == id)
      if (nrow(scen_row) == 0) return(NULL)
      
      launch_scenario_DE(
        scenario_id = id,
        scen = scen_row,
        nsim = NSIM,
        dir_out = dir_out,
        simula_fun_name = simula_fun_name
      )
    })
    
    # Aspetta TUTTI con timeout globale (DE può essere più lento)
    results <- wait_processes(procs, TIMEOUT_SEC * 2)
    batch_ok <- results$ok
    batch_timeout <- results$fail
    timeout_count <- timeout_count + batch_timeout
    
    cat(sprintf("%d OK, %d timeout/skip\n", batch_ok, batch_timeout))
  }
  
  cat(sprintf("\n%s completato.\n", phase_name))
}

# =============================================================================
# MAIN - INDIVIDUAL
# =============================================================================

cat("\n========== PAPER 2 POWER - INDIVIDUAL ==========\n")
start_total <- Sys.time()

process_scenarios_batch(scenarios, dir_out_ind, "simulate_survival_cohort_individual", "STEP 1: Power base Individual")

# Carica risultati base
base_files <- list.files(dir_out_ind, pattern = "^scenario_[0-9]+\\.rds$", full.names = TRUE)
if (length(base_files) > 0) {
  results_base_ind <- bind_rows(lapply(base_files, readRDS))
  
  # Step 2: DE
  process_scenarios_DE_batch(results_base_ind, dir_out_ind, "simulate_survival_cohort_individual", "STEP 2: Power DE Individual")
  
  # Combina risultati
  de_files <- list.files(dir_out_ind, pattern = "^scenario_DE_[0-9]+\\.rds$", full.names = TRUE)
  if (length(de_files) > 0) {
    results_DE_ind <- bind_rows(lapply(de_files, readRDS))
    results_final_ind <- results_base_ind %>%
      left_join(results_DE_ind %>% dplyr::select(scenario_id, power_DE, prop_cens_DE), by = "scenario_id")
  } else {
    results_final_ind <- results_base_ind %>% mutate(power_DE = NA_real_, prop_cens_DE = NA_real_)
  }
  
  results_final_ind <- results_final_ind %>% mutate(design = "individual")
  saveRDS(results_final_ind, file.path(dir_out_ind, "power_results_individual_FINAL.rds"))
  cat("\nSalvato: power_results_individual_FINAL.rds\n")
}

# =============================================================================
# MAIN - HOSPITAL
# =============================================================================

cat("\n========== PAPER 2 POWER - HOSPITAL ==========\n")

process_scenarios_batch(scenarios, dir_out_hosp, "simulate_survival_cohort_hospital", "STEP 1: Power base Hospital")

# Carica risultati base
base_files <- list.files(dir_out_hosp, pattern = "^scenario_[0-9]+\\.rds$", full.names = TRUE)
if (length(base_files) > 0) {
  results_base_hosp <- bind_rows(lapply(base_files, readRDS))
  
  # Step 2: DE
  process_scenarios_DE_batch(results_base_hosp, dir_out_hosp, "simulate_survival_cohort_hospital", "STEP 2: Power DE Hospital")
  
  # Combina risultati
  de_files <- list.files(dir_out_hosp, pattern = "^scenario_DE_[0-9]+\\.rds$", full.names = TRUE)
  if (length(de_files) > 0) {
    results_DE_hosp <- bind_rows(lapply(de_files, readRDS))
    results_final_hosp <- results_base_hosp %>%
      left_join(results_DE_hosp %>% dplyr::select(scenario_id, power_DE, prop_cens_DE), by = "scenario_id")
  } else {
    results_final_hosp <- results_base_hosp %>% mutate(power_DE = NA_real_, prop_cens_DE = NA_real_)
  }
  
  results_final_hosp <- results_final_hosp %>% mutate(design = "hospital")
  saveRDS(results_final_hosp, file.path(dir_out_hosp, "power_results_hospital_FINAL.rds"))
  cat("\nSalvato: power_results_hospital_FINAL.rds\n")
}

# =============================================================================
# COMBINA TUTTO
# =============================================================================

cat("\n=== COMBINAZIONE FINALE ===\n")

if (exists("results_final_ind") && exists("results_final_hosp")) {
  power_combined <- bind_rows(results_final_ind, results_final_hosp)
  saveRDS(power_combined, file.path(DIR_BASE, "paper2_POWER_ALL_RESULTS.rds"))
  cat("Salvato: paper2_POWER_ALL_RESULTS.rds\n")
}

elapsed_total <- as.numeric(difftime(Sys.time(), start_total, units = "hours"))
cat(sprintf("\n=== COMPLETATO in %.1f ore ===\n", elapsed_total))
