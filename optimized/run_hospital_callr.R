# =============================================================================
# run_hospital_callr.R - Esegue scenari Hospital con timeout REALE
# Una replica alla volta in processo separato, kill dopo 30s se non finisce
# =============================================================================

cat("\n=== HOSPITAL SEQUENTIAL CON TIMEOUT CALLR ===\n")
cat("Timeout per replica: 30 secondi\n")
cat("Se una replica non finisce in tempo, viene killata e si riprova\n\n")

library(callr)
source("01_config.R")
source("02_functions.R")

# Parametri
TIMEOUT_SEC <- 30
N_REP <- PAPER1_N_REP  # 1000
MAX_ATTEMPTS_MULTIPLIER <- 3  # Max 3x tentativi per evitare loop infiniti
N_PARALLEL <- 20  # Worker paralleli

# Directory output
dir_out <- file.path(DIR_BASE, "paper1_icc_hospital")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# Genera scenari
scenarios <- make_paper1_scenarios()
n_scenarios <- nrow(scenarios)

cat(sprintf("Totale scenari: %d\n", n_scenarios))
cat(sprintf("Repliche per scenario: %d\n", N_REP))
cat(sprintf("Timeout per replica: %d secondi\n", TIMEOUT_SEC))
cat(sprintf("Worker paralleli: %d\n\n", N_PARALLEL))

# Funzione per eseguire UNA replica in processo separato con timeout
run_single_replica <- function(scen, timeout_sec, simula_fun_name) {
  proc <- r_bg(
    function(scen, simula_fun_name) {
      suppressPackageStartupMessages({
        source("01_config.R")
        source("02_functions.R")
      })
      
      simula_fun <- get(simula_fun_name)
      
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
    },
    args = list(scen = scen, simula_fun_name = simula_fun_name),
    package = TRUE
  )
  
  start <- Sys.time()
  while (proc$is_alive()) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    if (elapsed > timeout_sec) {
      proc$kill()
      return(list(success = FALSE, reason = "TIMEOUT", time = elapsed))
    }
    Sys.sleep(0.2)
  }
  
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  
  tryCatch({
    result <- proc$get_result()
    list(success = TRUE, result = result, time = elapsed)
  }, error = function(e) {
    list(success = FALSE, reason = conditionMessage(e), time = elapsed)
  })
}

# Funzione per processare uno scenario
process_scenario_callr <- function(scenario_id, scen, n_rep, dir_out, timeout_sec) {
  
  summary_file <- file.path(dir_out, sprintf("scenario_%03d_summary.rds", scenario_id))
  
  # Skip se già completato
  if (file.exists(summary_file)) {
    cat(sprintf("Scenario %d già completato, skip.\n", scenario_id))
    return(NULL)
  }
  
  cat(sprintf("\n>>> Scenario %d | icc=%.2f, n=%d, hosp=%d, beta=%.1f, cens=%.1f, bal=%d\n",
              scenario_id, scen$icc, scen$sample_size, scen$num_hosp, 
              scen$beta, scen$cens, scen$balancing_mode))
  
  valid_results <- list()
  attempt <- 0
  timeout_count <- 0
  max_attempts <- n_rep * MAX_ATTEMPTS_MULTIPLIER
  
  start_scenario <- Sys.time()
  
  while (length(valid_results) < n_rep && attempt < max_attempts) {
    # Quante repliche servono ancora?
    needed <- min(N_PARALLEL, n_rep - length(valid_results))
    
    cat(sprintf("  Batch %d repliche (valide: %d/%d)... ", needed, length(valid_results), n_rep))
    
    # Lancia N_PARALLEL processi in parallelo
    procs <- lapply(1:needed, function(i) {
      r_bg(
        function(scen) {
          suppressPackageStartupMessages({
            source("01_config.R")
            source("02_functions.R")
          })
          cohort <- simulate_survival_cohort_hospital(
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
        },
        args = list(scen = scen),
        package = TRUE
      )
    })
    
    # Aspetta tutti con timeout
    start_batch <- Sys.time()
    batch_ok <- 0
    batch_timeout <- 0
    
    for (i in seq_along(procs)) {
      proc <- procs[[i]]
      
      # Aspetta questo processo con timeout
      while (proc$is_alive()) {
        elapsed <- as.numeric(difftime(Sys.time(), start_batch, units = "secs"))
        if (elapsed > timeout_sec) {
          proc$kill()
          break
        }
        Sys.sleep(0.1)
      }
      
      # Raccogli risultato
      res <- tryCatch({
        if (proc$is_alive()) {
          proc$kill()
          NULL
        } else {
          proc$get_result()
        }
      }, error = function(e) NULL)
      
      if (!is.null(res)) {
        valid_results[[length(valid_results) + 1]] <- res
        batch_ok <- batch_ok + 1
      } else {
        batch_timeout <- batch_timeout + 1
        timeout_count <- timeout_count + 1
      }
    }
    
    attempt <- attempt + needed
    cat(sprintf("%d OK, %d timeout\n", batch_ok, batch_timeout))
  }
  
  elapsed_scenario <- as.numeric(difftime(Sys.time(), start_scenario, units = "mins"))
  
  if (length(valid_results) < n_rep) {
    cat(sprintf("  ATTENZIONE: Solo %d/%d repliche valide (troppi timeout)\n", 
                length(valid_results), n_rep))
  }
  
  # Combina risultati
  if (length(valid_results) == 0) {
    cat(sprintf("  ERRORE: Scenario %d fallito completamente\n", scenario_id))
    return(NULL)
  }
  
  rep_results <- bind_rows(valid_results)
  
  # Calcola statistiche
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
      balancing_mode = scen$balancing_mode,
      total_attempts = attempt,
      timeout_count = timeout_count,
      elapsed_min = elapsed_scenario
    )
  
  saveRDS(summary_icc, summary_file)
  
  cat(sprintf("  Completato: %d/%d valide | %d timeout | %.1f min\n", 
              length(valid_results), n_rep, timeout_count, elapsed_scenario))
  
  return(summary_icc)
}

# =============================================================================
# ESECUZIONE PRINCIPALE
# =============================================================================

cat("Inizio elaborazione scenari...\n")
start_total <- Sys.time()

for (i in 1:n_scenarios) {
  scen <- scenarios[i, ]
  
  process_scenario_callr(
    scenario_id = i,
    scen = scen,
    n_rep = N_REP,
    dir_out = dir_out,
    timeout_sec = TIMEOUT_SEC
  )
  
  # Stima tempo rimanente ogni 10 scenari
  if (i %% 10 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_total, units = "hours"))
    completed <- length(list.files(dir_out, pattern = "_summary.rds"))
    remaining <- n_scenarios - completed
    if (completed > 0) {
      rate <- elapsed / completed
      eta <- remaining * rate
      cat(sprintf("\n=== PROGRESSO: %d/%d completati | ETA: %.1f ore ===\n\n", 
                  completed, n_scenarios, eta))
    }
  }
}

elapsed_total <- as.numeric(difftime(Sys.time(), start_total, units = "hours"))
cat(sprintf("\n=== COMPLETATO in %.1f ore ===\n", elapsed_total))
