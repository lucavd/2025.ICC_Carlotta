library(tidyverse)
# library(ggplot2)
# library(patchwork)
# library(flextable)
library(parallel)

##### Simulate scenarios here; tables and graphs to be added after datasets are created

#source("funzioni paperSeptember2025.R")
# source("funzioni paper - Ott 2025.R")
source("Paper Functions – Nov 2025.R")


# DB RESULTS -----------------------------------------------------

# Final scenario parameters

betas <- c(- 0.2, - 0.5) ## corrispondono a HR: 0.8 medio, 0.6 forte
lambdas <- c(0.115, 0.012) #means median: 6 months, 60 months
#num_pat_groups <- c(10, 20, 40, 60, 100) # removed num_pat_group; using total sample size instead
sample_sizes <- c(100, 500, 1000, 6000)  ### c(50, 100, 150, 200, 300, 500, 600, 900, 1200, 1500, 1800, 2400, 3000, 3600, 6000) ## tutte le combinazioni hosp*pazienti
num_hosps <- c(5, 15, 60) ## tolgo 30
iccs <- c(0.01, 0.06, 0.1, 0.3) ### farei solo 3 icc basso medio alto - ricontrolla da letteratura
balancing_modes <- c(1, 2) ##new add !
cens_value <- c( 0.5, 5) ##tolgo 2 

###2 times this evaluation

# ## test:
# betas <- c(-0.5, 0.5)
# lambdas <- c(0.11)
# # num_pat_groups <- c(10, 100)
# sample_sizes <- c(100, 1000)
# num_hosps <- c(5)
# iccs <- c(0.01, 0.05, 0.25)
# cens_value <- 2 
# balancing_modes <- 2


# scenario grid
### def 720 scenarios
scenarios <- expand.grid(
  beta = betas,
  lambda = lambdas,
  #num_pat_group = num_pat_groups,
  sample_size = sample_sizes,
  num_hosp = num_hosps,
  icc = iccs,
  cens = cens_value,
  balancing_mode = balancing_modes)


# OPZIONE A BOOTSTRAP PER OSPEDALE - CLUSTER-----------------------------------------------------

## I would remove this option A) and keep only Monte Carlo option B (see below) for the simulation, if you agree...


# PURPOSE OF THIS SCRIPT
# This script simulates multiple survival scenarios by varying key parameters
# (number of clusters, ICC, treatment effect, baseline hazard, censoring, etc.).
#
# For each scenario, the following steps are performed:
#   1) simulate a survival cohort using simulate_survival_cohort_individual() OR simulate_survival_cohort_hospital ## this last one, if we decide to keep it
#   2) estimate point ICC with different methods
#   3) compute bootstrap confidence intervals
#
# Since the number of scenarios can be large and each bootstrap procedure is
# computationally intensive, parallel computing is used.
#
# To ensure computational safety and avoid losing progress in case of errors or
# interruptions, each scenario is saved individually 

# # DIRECTORY PARTIAL RESULTS ---
# unlink("results/ICC_BT_results_paper1", recursive = TRUE) ## delete, per testing
# dir.create("results/ICC_BT_results_paper1", showWarnings = FALSE)
# 
# # function to process a single scenario 
# process_scenario <- function(i, scenarios) {
#   scen <- scenarios[i, ]
#   filename <- paste0("results/ICC_BT_results_paper1/scenario_", i, ".rds")
#   
#   if(file.exists(filename)) {
#     cat("Scenario", i, "already exists. Skipping.\n")
#     return(NULL)  }
#   
#   cat("Processing scenario:", i, "of", nrow(scenarios), "\n")
#   
#   cohort <- tryCatch({
#     simulate_survival_cohort_individual( ## change if including second design type: simulate_survival_cohort_hospital
#       num_hosp = scen$num_hosp,
#       sample_size = scen$sample_size,
#       icc = scen$icc,
#       pop_treat_effect = scen$beta,
#       lambda = scen$lambda,
#       gammas = 1,
#       cens = scen$cens,
#       balancing_mode = scen$balancing_mode)
#   }, error = function(e) {
#     cat("ERROR in scenario", i, "- skipping cohort.\n")
#     return(NULL)   })
#   
#   if(is.null(cohort)) return(NULL)
#   
#   boot_res <- tryCatch({
#     bootstrap_icc(cohort, B = 10)  # B = 2000 per versione definitiva
#   }, error = function(e) {
#     cat("ERROR in bootstrap scenario", i, "- skipping.\n")
#     return(NULL) })
#   
#   if(is.null(boot_res)) return(NULL)
#   
#   boot_summary <- boot_res$summary %>%
#     mutate(
#       beta = scen$beta,
#       lambda = scen$lambda,
#       sample_size = scen$sample_size,
#       num_pat_group_mean = scen$sample_size/scen$num_hosp,
#       num_hosp = scen$num_hosp,
#       icc_input = scen$icc,
#       cens = scen$cens,
#       prop_cens = mean(cohort$status == 0, rm.na = TRUE),
#       balancing_mode <- scen$balancing_mode)
#   
#   saveRDS(boot_summary, file = filename)
#   cat("Saved scenario", i, "\n")
#   
#   return(NULL) }
# 
# num_cores <- detectCores() - 1
# cl <- makeCluster(num_cores)
# 
# clusterExport(cl, varlist = c(
#   "scenarios",
#   "simulate_survival_cohort_individual", ## change if including second design type: simulate_survival_cohort_hospital
#   "bootstrap_icc",
#   "icc_estimation" ))
# 
# clusterEvalQ(cl, {
#   library(dplyr)
#   library(simsurv)
#   library(survival)
#   library(lme4)
#   library(coxme) })
# 
# # parallelizzation
# parLapply(cl, 1:nrow(scenarios), process_scenario, scenarios = scenarios)
# 
# stopCluster(cl)
# 
# # all results
# files <- list.files("results/ICC_BT_results_paper1", pattern = "*.rds", full.names = TRUE)
# results_list <- lapply(files, readRDS)
# final_icc_table_BT <- bind_rows(results_list)
# 
# 
# 
# 


# OPZIONE B – MONTE CARLO REPLICATES (NO BOOTSTRAP) -----------------------------------------------
### I'll keep only this option for the simulation. Do you agree?


# PURPOSE OF THIS SCRIPT
# This script performs a Monte Carlo simulation study across multiple survival scenarios

# For each scenario, the following steps are performed:
#   1) simulate a survival cohort using `simulate_survival_cohort_individual()`OR simulate_survival_cohort_hospital ## this last one, if we decide to keep it
#   2) estimate the ICC using multiple statistical methods 
#   3) repeat the process n_rep times (Monte Carlo replications)
#   4) compute the mean ICC, standard deviation, and 95% confidence interval 
#
# Unlike OPTION A this version does NOT use bootstrap resampling. 
# Each repetition corresponds to a new independently simulated cohort.


n_rep <- 2000  
dir_out <-"results/ICC_MC_results_paper1_prova3"

unlink(dir_out, recursive = TRUE) ##delete, for test
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# total number of scenarios
n_scen <- nrow(scenarios)

process_scenario <- function(i, scenarios, n_rep, dir_out) {
  scen <- scenarios[i, ]
  file_out <- file.path(dir_out, paste0("scenario_", i, ".rds"))
  
  if (file.exists(file_out)) {
    cat("Scenario", i, "alredy existing. Skipping.\n")
    return(NULL)  }
  
  cat(">>> Scenario", i, "su", nrow(scenarios), " | beta =", scen$beta,
      "icc =", scen$icc, "num_pat_group =", scen$num_pat_group, "\n")
  
  # --- loop repliche ---
  rep_results <- map_dfr(1:n_rep, function(r) {
    cohort <- tryCatch({
      simulate_survival_cohort_individual( ## change if including second design type: simulate_survival_cohort_hospital
        num_hosp = scen$num_hosp,
        sample_size = scen$sample_size,
        balancing_mode = scen$balancing_mode,   # sbilanciato, puoi parametrizzarlo
        icc = scen$icc,
        pop_treat_effect = scen$beta,
        lambda = scen$lambda,
        gammas = 1,
        cens = scen$cens      )
    }, error = function(e) NULL)
    
    if (is.null(cohort)) return(NULL)
    
    # ICC estimation
    stima_icc <- tryCatch({
      icc_estimation(cohort)
    }, error = function(e) NULL)
    
    if (is.null(stima_icc)) return(NULL)
    
    stima_icc$prop_cens <- mean(cohort$status == 0)
    stima_icc$rep <- r
    return(stima_icc)  })
  
  if (nrow(rep_results) == 0) {
    cat("Scenario", i, "failed – no valid replication.\n")
    return(NULL)  }
  
  # summaRy per scenario
  summary_icc <- rep_results %>%
    group_by(Method) %>%
    summarise(
      Mean_ICC = mean(ICC, na.rm = TRUE),
      SD_ICC = sd(ICC, na.rm = TRUE),
      Lower_95CI = Mean_ICC - 1.96 * SD_ICC / sqrt(n_rep),    ## QUESTION: Mean & SD vs. Quantile? Likely similar due to many repetitions. Which do you prefer?
      Upper_95CI = Mean_ICC + 1.96 * SD_ICC / sqrt(n_rep),
      Lower_95CI_Median = quantile(ICC, probs = 0.025, na.rm = TRUE),
      Upper_95CI_Median = quantile(ICC, probs = 0.975, na.rm = TRUE),
      Median_ICC = quantile(ICC, probs = 0.5, na.rm = TRUE),
      prop_cens = mean(prop_cens, na.rm = TRUE),
      .groups = "drop"   ) %>%
    mutate(
      beta = scen$beta,
      lambda = scen$lambda,
      sample_size = scen$sample_size,
      num_pat_group = scen$sample_size/scen$num_hosp,
      num_hosp = scen$num_hosp,
      icc_input = scen$icc,
      cens = scen$cens,
      balancing_mode = scen$balancing_mode)
  
  # --- save results ---
  saveRDS(summary_icc, file = file_out)
  cat("Scenario", i, "salvato in", file_out, "\n")
  return(NULL) }


# PARALLEL CLUSTER

num_cores <- max(1, detectCores() - 1)
cl <- makeCluster(num_cores)
cat("Usando", num_cores, "core.\n")

clusterExport(cl, varlist = c(
  "simulate_survival_cohort_individual", ## change if including second design type: simulate_survival_cohort_hospital
  "icc_estimation",
  "scenarios",
  "n_rep",
  "dir_out",
  "process_scenario" ))

clusterEvalQ(cl, {
  library(dplyr)
  library(purrr)
  library(survival)
  library(simsurv)
  library(lme4)
  library(coxme)})


parLapply(cl, 1:n_scen, process_scenario,
          scenarios = scenarios,
          n_rep = n_rep,
          dir_out = dir_out)

stopCluster(cl)


# results FINAL
files <- list.files(dir_out, pattern = "*.rds", full.names = TRUE)
results_list <- lapply(files, readRDS)
final_icc_table_MC <- bind_rows(results_list)

saveRDS(final_icc_table, file = file.path("results", "ICC_MC_results_paper1"))








# GRAFICI -----------------------------------------------------------------
### prendi spunto da Kalia et al (2016)
