# =============================================================================
# run_paper2_power_callr.R - Paper 2 Power con timeout REALE via callr
# Pool dinamico: 10 worker, timeout 15 min per processo
# =============================================================================

cat("\n=== PAPER 2 POWER CON TIMEOUT CALLR ===\n")
cat("Pool dinamico: 10 worker\n")
cat("Timeout per processo: 15 minuti\n\n")

library(callr)
library(dplyr)
source("01_config.R")
source("02_functions.R")

# Parametri
TIMEOUT_SEC <- 900  # 15 minuti per processo
NSIM <- PAPER2_POWER_NSIM  # 1000
MAX_WORKERS <- 10

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
cat(sprintf("Timeout per processo: %d secondi\n", TIMEOUT_SEC))
cat(sprintf("Max worker paralleli: %d\n\n", MAX_WORKERS))

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

# Pool dinamico: gestisce max_workers processi, lancia nuovo appena uno finisce
run_pool <- function(scenario_ids, scenarios, dir_out, simula_fun_name, max_workers, timeout_sec) {
  
  # Stato del pool
  queue <- scenario_ids  # Coda scenari da processare
  active <- list()       # Lista processi attivi: list(proc=, id=, start=)
  completed <- 0
  failed <- 0
  
  start_total <- Sys.time()
  
  while (length(queue) > 0 || length(active) > 0) {
    
    # 1. Lancia nuovi processi se c'è spazio
    while (length(active) < max_workers && length(queue) > 0) {
      id <- queue[1]
      queue <- queue[-1]
      
      proc <- launch_scenario(id, scenarios[id,], NSIM, dir_out, simula_fun_name)
      
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
        # Processo terminato
        res <- tryCatch(w$proc$get_result(), error = function(e) NULL)
        if (!is.null(res)) {
          completed <- completed + 1
        } else {
          failed <- failed + 1
        }
      } else if (elapsed > timeout_sec) {
        # Timeout
        w$proc$kill()
        failed <- failed + 1
        cat(sprintf(" [T%d]", w$id))
      } else {
        # Ancora attivo
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
    
    Sys.sleep(0.5)
  }
  
  cat("\n")
  list(completed = completed, failed = failed)
}

# Pool dinamico per scenari DE
run_pool_DE <- function(scenario_ids, scenarios_DE, dir_out, simula_fun_name, max_workers, timeout_sec) {
  
  queue <- scenario_ids
  active <- list()
  completed <- 0
  failed <- 0
  
  start_total <- Sys.time()
  
  while (length(queue) > 0 || length(active) > 0) {
    
    # 1. Lancia nuovi processi
    while (length(active) < max_workers && length(queue) > 0) {
      id <- queue[1]
      queue <- queue[-1]
      
      scen_row <- scenarios_DE %>% filter(scenario_id == id)
      if (nrow(scen_row) == 0) {
        completed <- completed + 1
        next
      }
      
      proc <- launch_scenario_DE(id, scen_row, NSIM, dir_out, simula_fun_name)
      
      if (!is.null(proc)) {
        active[[length(active) + 1]] <- list(proc = proc, id = id, start = Sys.time())
      } else {
        completed <- completed + 1
      }
    }
    
    if (length(active) == 0) break
    
    # 2. Controlla processi
    still_active <- list()
    for (w in active) {
      elapsed <- as.numeric(difftime(Sys.time(), w$start, units = "secs"))
      
      if (!w$proc$is_alive()) {
        res <- tryCatch(w$proc$get_result(), error = function(e) {
          cat(sprintf(" [E%d: %s]", w$id, conditionMessage(e)))
          NULL
        })
        if (!is.null(res)) completed <- completed + 1 else failed <- failed + 1
      } else if (elapsed > timeout_sec) {
        w$proc$kill()
        failed <- failed + 1
        cat(sprintf(" [T%d]", w$id))
      } else {
        still_active[[length(still_active) + 1]] <- w
      }
    }
    active <- still_active
    
    # 3. Progresso
    total_done <- completed + failed
    total <- length(scenario_ids)
    elapsed_min <- as.numeric(difftime(Sys.time(), start_total, units = "mins"))
    cat(sprintf("\r  [%.1fm] %d/%d (ok:%d fail:%d) active:%d    ",
                elapsed_min, total_done, total, completed, failed, length(active)))
    
    Sys.sleep(0.5)
  }
  
  cat("\n")
  list(completed = completed, failed = failed)
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

process_scenarios_pool <- function(scenarios, dir_out, simula_fun_name, phase_name) {
  
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
  
  # Usa pool dinamico
  results <- run_pool(to_process, scenarios, dir_out, simula_fun_name, MAX_WORKERS, TIMEOUT_SEC)
  
  cat(sprintf("\n%s completato: %d OK, %d timeout\n", phase_name, results$completed, results$failed))
}

process_scenarios_DE_pool <- function(base_results, dir_out, simula_fun_name, phase_name) {
  
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
  
  # Pool dinamico per DE (timeout più lungo)
  results <- run_pool_DE(to_process_ids, scenarios_DE, dir_out, simula_fun_name, MAX_WORKERS, TIMEOUT_SEC * 2)
  
  cat(sprintf("\n%s completato: %d OK, %d timeout\n", phase_name, results$completed, results$failed))
}

# =============================================================================
# MAIN - INDIVIDUAL
# =============================================================================

cat("\n========== PAPER 2 POWER - INDIVIDUAL ==========\n")
start_total <- Sys.time()

process_scenarios_pool(scenarios, dir_out_ind, "simulate_survival_cohort_individual", "STEP 1: Power base Individual")

# Carica risultati base
base_files <- list.files(dir_out_ind, pattern = "^scenario_[0-9]+\\.rds$", full.names = TRUE)
if (length(base_files) > 0) {
  results_base_ind <- bind_rows(lapply(base_files, readRDS))
  
  # Step 2: DE
  process_scenarios_DE_pool(results_base_ind, dir_out_ind, "simulate_survival_cohort_individual", "STEP 2: Power DE Individual")
  
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

process_scenarios_pool(scenarios, dir_out_hosp, "simulate_survival_cohort_hospital", "STEP 1: Power base Hospital")

# Carica risultati base
base_files <- list.files(dir_out_hosp, pattern = "^scenario_[0-9]+\\.rds$", full.names = TRUE)
if (length(base_files) > 0) {
  results_base_hosp <- bind_rows(lapply(base_files, readRDS))
  
  # Step 2: DE
  process_scenarios_DE_pool(results_base_hosp, dir_out_hosp, "simulate_survival_cohort_hospital", "STEP 2: Power DE Hospital")
  
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
