source("/home/user/2025.ICC_carlotta/optimized/01_config.R")
source("/home/user/2025.ICC_carlotta/optimized/02_functions.R")

library(future.apply)

# ===============================
# INDIVIDUAL - 241 scenari mancanti (con salvataggio parziale)
# ===============================

scenarios <- read_csv("scenari_individual_mancanti.csv")
dir_out <- "individual_mancanti_partial"
dir.create(dir_out, showWarnings = FALSE)

required_cols <- c(
  "scenario_id", "num_hosp", "sample_size", "balancing_mode",
  "icc", "pop_treat_effect", "lambda", "cens")
stopifnot(all(required_cols %in% names(scenarios)))

cat(sprintf("Scenari individual: %d\n", nrow(scenarios)))

run_one_scenario <- function(row, nsim = 1000, dir_out) {
  
  filename <- file.path(dir_out, sprintf("scenario_%03d.rds", row$scenario_id))
  
  if (file.exists(filename)) {
    return(readRDS(filename))
  }
  
  cat(sprintf(">>> Scenario %d | hosp=%d, ss=%d, icc=%.2f\n", 
              row$scenario_id, row$num_hosp, row$sample_size, row$icc))
  
  res <- surv_power_function(
    simula_coorte_fun = simulate_survival_cohort_individual,
    simula_args = list(
      num_hosp = row$num_hosp,
      sample_size = row$sample_size,
      balancing_mode = row$balancing_mode,
      icc = row$icc,
      pop_treat_effect = row$pop_treat_effect,
      lambda = row$lambda,
      gammas = 1,
      cens = row$cens
    ),
    nsim = nsim
  ) %>%
    mutate(
      scenario_id = row$scenario_id,
      balancing_mode = row$balancing_mode,
      pop_treat_effect = row$pop_treat_effect,
      lambda = row$lambda,
      cens = row$cens
    )
  
  saveRDS(res, filename)
  return(res)
}

set.seed(123)

plan(multisession, workers = N_CORES)

cat("Avvio simulazioni individual...\n")

results <- future_lapply(
  split(scenarios, scenarios$scenario_id),
  function(s) run_one_scenario(s[1, ], nsim = 1000, dir_out = dir_out),
  future.seed = TRUE
) %>%
  bind_rows()

plan(sequential)

saveRDS(results, "surv_power_individual_mancanti.rds")
cat("Salvato: surv_power_individual_mancanti.rds\n")
cat(sprintf("Parziali in: %s/\n", dir_out))
