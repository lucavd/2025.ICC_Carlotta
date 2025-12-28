# =============================================================================
# run_parallel_hospital.R - Esegue un subset di scenari hospital
# =============================================================================
# USO: Rscript run_parallel_hospital.R <worker_id> <total_workers>
# Esempio: 4 terminali -> Rscript run_parallel_hospital.R 1 4
#                         Rscript run_parallel_hospital.R 2 4
#                         Rscript run_parallel_hospital.R 3 4
#                         Rscript run_parallel_hospital.R 4 4
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("USO: Rscript run_parallel_hospital.R <worker_id> <total_workers>")
}

worker_id <- as.integer(args[1])
total_workers <- as.integer(args[2])

cat(sprintf("\n=== WORKER %d/%d - HOSPITAL ===\n\n", worker_id, total_workers))

source("01_config.R")
source("02_functions.R")

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
    tryCatch({
      R.utils::withTimeout({
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
      }, timeout = 300)  # 5 minuti
    }, TimeoutException = function(e) {
      message("Replica timeout dopo 5 minuti - skip")
      NULL
    })
  }
  
  rep_results <- run_with_chunks(
    scenario_id = sprintf("%03d", scenario_id),
    n_rep = n_rep,
    chunk_size = chunk_size,
    dir_out = dir_out,
    rep_fun = single_rep
  )
  
  if (nrow(rep_results) == 0) {
    message(sprintf("Scenario %d fallito.", scenario_id))
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

# Directory output
dir_out <- file.path(DIR_BASE, "paper1_icc_hospital")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# Genera scenari
scenarios <- make_paper1_scenarios()
n_scenarios <- nrow(scenarios)

# Dividi scenari tra workers
my_scenarios <- seq(worker_id, n_scenarios, by = total_workers)

cat(sprintf("Totale scenari: %d\n", n_scenarios))
cat(sprintf("Scenari per questo worker: %d\n", length(my_scenarios)))
cat(sprintf("Scenari assegnati: %s\n\n", paste(head(my_scenarios, 10), collapse=", ")))

# Esegui scenari assegnati
for (i in my_scenarios) {
  scen <- scenarios[i, ]
  
  result <- process_icc_scenario(
    scenario_id = i,
    scen = scen,
    n_rep = PAPER1_N_REP,
    chunk_size = PAPER1_CHUNK_SIZE,
    dir_out = dir_out,
    simula_fun = simulate_survival_cohort_hospital
  )
}

cat(sprintf("\n=== WORKER %d COMPLETATO ===\n", worker_id))
