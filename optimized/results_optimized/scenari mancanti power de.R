library(survival)
library(tidyverse)
library(parallel)
library(future)
library(future.apply)

##INDIVIDUAL
# ===============================
# 1) CARICA SCENARI
# ===============================

scenarios <- read_csv("scenari_individual_mancanti.csv")

required_cols <- c(
  "scenario_id", "num_hosp", "sample_size", "balancing_mode",
  "icc", "pop_treat_effect", "lambda", "cens")
stopifnot(all(required_cols %in% names(scenarios)))

# ===============================
# 2) FUNZIONE: UN SINGOLO SCENARIO
# ===============================

run_one_scenario <- function(row, nsim = 1000) {
  
  surv_power_function(
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
}

# ===============================
# 3) ESECUZIONE DI TUTTI GLI SCENARI
# ===============================

set.seed(123)

plan(multisession, workers = max(1, availableCores() - 1))

results <- future_lapply(
  split(scenarios, scenarios$scenario_id),
  function(s) run_one_scenario(s[1, ], nsim = 1000),
  future.seed = TRUE) %>%
  bind_rows()

plan(sequential)

# ===============================
# 4) SALVATAGGIO
# ===============================
saveRDS(results, "surv_power_all_scenarios.rds")


##INDIVIDUAL
# ===============================
# 1) CARICA SCENARI
# ===============================

scenarios_hosp <- read_csv("scenari_hospital_mancanti.csv")

required_cols <- c(
  "scenario_id", "num_hosp", "sample_size", "balancing_mode",
  "icc", "pop_treat_effect", "lambda", "cens")
stopifnot(all(required_cols %in% names(scenarios_hosp)))

# ===============================
# 2) FUNZIONE: UN SINGOLO SCENARIO
# ===============================

run_one_scenario <- function(row, nsim = 1000) {
  
  surv_power_function(
    simula_coorte_fun = simulate_survival_cohort_hospital(),
    simula_args = list(
      num_hosp = row$num_hosp,
      sample_size = row$sample_size,
      balancing_mode = row$balancing_mode,
      icc = row$icc,
      pop_treat_effect = row$pop_treat_effect,
      lambda = row$lambda,
      gammas = 1,
      cens = row$cens  ),
    nsim = nsim  ) %>%
    mutate(
      scenario_id = row$scenario_id,
      balancing_mode = row$balancing_mode,
      pop_treat_effect = row$pop_treat_effect,
      lambda = row$lambda,
      cens = row$cens )}

# ===============================
# 3) ESECUZIONE DI TUTTI GLI SCENARI
# ===============================

set.seed(123)

plan(multisession, workers = max(1, availableCores() - 1))

results_hosp <- future_lapply(
  split(scenarios_hosp, scenarios_hosp$scenario_id),
  function(s) run_one_scenario(s[1, ], nsim = 500),
  future.seed = TRUE ) %>%
  bind_rows()

plan(sequential)

# ===============================
# 4) SALVATAGGIO
# ===============================
saveRDS(results_hosp, "surv_power_all_scenarios_hosp.rds")
