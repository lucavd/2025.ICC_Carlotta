# =============================================================================
# run_paper2_ss_callr.R - Paper 2 Sample Size con timeout REALE via callr
# Pool dinamico: 10 worker, timeout configurabile per processo
# =============================================================================

cat("\n=== PAPER 2 SAMPLE SIZE CON TIMEOUT CALLR ===\n")
cat("Pool dinamico: 10 worker\n")
cat("Timeout per processo: da config\n\n")

library(callr)
library(dplyr)
source("01_config.R")
source("02_functions.R")

# Parametri
TIMEOUT_SEC <- min(PAPER2_SS_TIMEOUT, 3600)  # max 1 ora per scenario
NSIM <- PAPER2_SS_NSIM
MAX_PAT <- PAPER2_SS_MAX_PAT
MAX_WORKERS <- 10

# Directory output
dir_out_ind <- DIR_PAPER2_SS
dir_out_hosp <- file.path(DIR_BASE, "paper2_sample_size_hospital")
dir.create(dir_out_ind, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out_hosp, recursive = TRUE, showWarnings = FALSE)

# Genera scenari
scenarios <- make_paper2_ss_scenarios()
n_scenarios <- nrow(scenarios)

cat(sprintf("Totale scenari: %d\n", n_scenarios))
cat(sprintf("nsim per step binary search: %d\n", NSIM))
cat(sprintf("Max pazienti per gruppo: %d\n", MAX_PAT))
cat(sprintf("Timeout per processo: %d secondi\n", TIMEOUT_SEC))
cat(sprintf("Max worker paralleli: %d\n\n", MAX_WORKERS))

# =============================================================================
# FUNZIONE PER ESEGUIRE UN SINGOLO SCENARIO SAMPLE SIZE
# =============================================================================

launch_ss_scenario <- function(scenario_id, scen, nsim, max_n, dir_out, simula_fun_name) {
  
  filename <- file.path(dir_out, sprintf("scenario_%03d.rds", scenario_id))
  if (file.exists(filename)) return(NULL)  # già completato
  
  r_bg(
    function(scenario_id, scen, nsim, max_n, dir_out, simula_fun_name) {
      suppressPackageStartupMessages({
        source("01_config.R")
        source("02_functions.R")
      })
      
      simula_fun <- get(simula_fun_name)
      
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
      }, error = function(e) NULL)
      
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
          power = NA_real_,
          cv = NA_real_
        )
      }
      
      res <- res %>%
        mutate(scenario_id = scenario_id)
      
      filename <- file.path(dir_out, sprintf("scenario_%03d.rds", scenario_id))
      saveRDS(res, filename)
      
      return(res)
    },
    args = list(
      scenario_id = scenario_id, 
      scen = scen, 
      nsim = nsim, 
      max_n = max_n,
      dir_out = dir_out, 
      simula_fun_name = simula_fun_name
    ),
    package = TRUE
  )
}

# =============================================================================
# POOL DINAMICO
# =============================================================================

run_pool_ss <- function(scenario_ids, scenarios, dir_out, simula_fun_name, max_workers, timeout_sec) {
  
  queue <- scenario_ids
  active <- list()
  completed <- 0
  failed <- 0
  
  start_total <- Sys.time()
  
  while (length(queue) > 0 || length(active) > 0) {
    
    # 1. Lancia nuovi processi se c'è spazio
    while (length(active) < max_workers && length(queue) > 0) {
      id <- queue[1]
      queue <- queue[-1]
      
      proc <- launch_ss_scenario(id, scenarios[id,], NSIM, MAX_PAT, dir_out, simula_fun_name)
      
      if (!is.null(proc)) {
        active[[length(active) + 1]] <- list(
          proc = proc,
          id = id,
          start = Sys.time()
        )
      } else {
        completed <- completed + 1  # Era già completato (cached)
      }
    }
    
    if (length(active) == 0) break
    
    # 2. Controlla processi attivi
    still_active <- list()
    for (w in active) {
      elapsed <- as.numeric(difftime(Sys.time(), w$start, units = "secs"))
      
      if (!w$proc$is_alive()) {
        res <- tryCatch(w$proc$get_result(), error = function(e) NULL)
        if (!is.null(res)) {
          completed <- completed + 1
        } else {
          failed <- failed + 1
        }
      } else if (elapsed > timeout_sec) {
        w$proc$kill()
        failed <- failed + 1
        cat(sprintf(" [T%d]", w$id))
      } else {
        still_active[[length(still_active) + 1]] <- w
      }
    }
    active <- still_active
    
    # 3. Mostra progresso
    total_done <- completed + failed
    total <- length(scenario_ids)
    elapsed_min <- as.numeric(difftime(Sys.time(), start_total, units = "mins"))
    cat(sprintf("\r  [%.1fm] %d/%d (ok:%d fail:%d) active:%d queue:%d    ",
                elapsed_min, total_done, total, completed, failed, 
                length(active), length(queue)))
    
    Sys.sleep(1)
  }
  
  cat("\n")
  list(completed = completed, failed = failed)
}

# =============================================================================
# PROCESSAMENTO
# =============================================================================

process_ss_scenarios <- function(scenarios, dir_out, simula_fun_name, phase_name) {
  
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
  
  results <- run_pool_ss(to_process, scenarios, dir_out, simula_fun_name, MAX_WORKERS, TIMEOUT_SEC)
  
  cat(sprintf("\n%s completato: %d OK, %d timeout/fail\n", phase_name, results$completed, results$failed))
}

# =============================================================================
# MAIN - INDIVIDUAL
# =============================================================================

cat("\n========== PAPER 2 SAMPLE SIZE - INDIVIDUAL ==========\n")
start_total <- Sys.time()

process_ss_scenarios(scenarios, dir_out_ind, "simulate_survival_cohort_individual", "Sample Size Individual")

# Carica e combina risultati
ind_files <- list.files(dir_out_ind, pattern = "^scenario_[0-9]+\\.rds$", full.names = TRUE)
if (length(ind_files) > 0) {
  results_ind <- bind_rows(lapply(ind_files, readRDS))
  results_ind <- results_ind %>% mutate(design = "individual")
  saveRDS(results_ind, file.path(dir_out_ind, "sample_size_individual_FINAL.rds"))
  cat("\nSalvato: sample_size_individual_FINAL.rds\n")
}

# =============================================================================
# MAIN - HOSPITAL
# =============================================================================

cat("\n========== PAPER 2 SAMPLE SIZE - HOSPITAL ==========\n")

process_ss_scenarios(scenarios, dir_out_hosp, "simulate_survival_cohort_hospital", "Sample Size Hospital")

# Carica e combina risultati
hosp_files <- list.files(dir_out_hosp, pattern = "^scenario_[0-9]+\\.rds$", full.names = TRUE)
if (length(hosp_files) > 0) {
  results_hosp <- bind_rows(lapply(hosp_files, readRDS))
  results_hosp <- results_hosp %>% mutate(design = "hospital")
  saveRDS(results_hosp, file.path(dir_out_hosp, "sample_size_hospital_FINAL.rds"))
  cat("\nSalvato: sample_size_hospital_FINAL.rds\n")
}

# =============================================================================
# COMBINA TUTTO
# =============================================================================

cat("\n=== COMBINAZIONE FINALE ===\n")

if (exists("results_ind") && exists("results_hosp")) {
  ss_combined <- bind_rows(results_ind, results_hosp)
  saveRDS(ss_combined, file.path(DIR_BASE, "paper2_SAMPLE_SIZE_ALL_RESULTS.rds"))
  cat("Salvato: paper2_SAMPLE_SIZE_ALL_RESULTS.rds\n")
}

elapsed_total <- as.numeric(difftime(Sys.time(), start_total, units = "hours"))
cat(sprintf("\n=== COMPLETATO in %.1f ore ===\n", elapsed_total))
