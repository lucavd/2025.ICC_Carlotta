library(simsurv)
library(survival)
library(tidyverse)
library(lme4)
library(ggplot2)
library(progressr)
library(R.utils)
library(parallel)

#set.seed(123)
source("Paper Functions – Nov 2025.R")

# 2 PAPER: Impact of ICC on Power and Sample Size --------

# 2) ICC E SAMPLE SIZE -------------------------------------------------------------
### How different levels of ICC influence statistical power and required sample.
## In this code, a loop is used to estimate the required sample size to achieve a target power 
## of 0.8 across different scenarios.


# Final scenario parameters
pte_values <- c(- 0.2, - 0.5) ## ==  HR: 0.8 medium, 0.6 high effect
lambda_values <- c(0.115, 0.012) #means median: 6 months, 60 months
#num_pat_groups <- c(10, 20, 40, 60, 100)
#sample_sizes <- c(100, 500, 1000, 2000, 6000)  
n_hosp_values <-  c(5, 60) #c(5, 15, 30, 60)
icc_values <- c(0, 0.01, 0.04) #seq(0.00, 0.32, by = 0.04)
cens_value <- c(0.5, 5) #tolto2
balancing_mode <- c(1,2)

# ##FOR TESTING
# icc_values <- c(0.1, 0.3, 0.5)#seq(0.00, 0.2, by = 0.04)
# pte_values <- c(-0.5, 0.5)
# n_hosp_values <- 50
# lambda_values <- c(0.01, 0.1)
# cens_value <- 2 #c(0.5, 5)
# #num_pat_group_values <- 50
# #sample_sizes <- n_hosp_values *3
# balancing_mode <- c(1, 2)  # added the option for unbalanced groups (== 2); starts from an mean n_int (not same n_init for group!!!). 
#                            # IS IT WORTH? IS IT INTERESTING?


scenarios <- expand.grid(
  pop_treat_effect = pte_values,
  lambda = lambda_values,
  icc = icc_values,
  num_hosp = n_hosp_values,
  cens = cens_value,
  balancing_mode = balancing_mode,
  stringsAsFactors = FALSE)


unlink("results/sample_size_results", recursive = TRUE) ## delete for testing
dir.create("results/sample_size_results", showWarnings = FALSE)


n_cores <- detectCores() - 1

cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(dplyr)
  library(simsurv)
  library(survival)
  library(R.utils)})

clusterExport(cl, c("sample_size_icc", "simulate_survival_cohort_individual"))

## parallel scenarios
parLapply(cl, seq_len(nrow(scenarios)), function(i, scenarios) {
  scen <- scenarios[i, ]
 
  filename <- paste0("results/sample_size_results/scenario_",
                     scen$pop_treat_effect, "_",
                     scen$lambda, "_",
                     scen$num_hosp, "_",
                     scen$cens, "_",
                     scen$balancing_mode, "_",
                     scen$icc, ".rds")
  
  if(file.exists(filename)) return(NULL)
  
  res <- tryCatch({
    sample_size_icc(
      icc = scen$icc,
      balancing_mode = scen$balancing_mode,
      num_hosp = scen$num_hosp,
      pop_treat_effect = scen$pop_treat_effect,
      nsim = 2000,  ## 2000 for REAL simulation
      lambda = scen$lambda,
      n_init = 2    )
  }, error = function(e) return(NULL))
  
  if(!is.null(res)) saveRDS(res, filename)
  return(res)
}, scenarios = scenarios)

stopCluster(cl)


## all final results
res_files <- list.files("results/sample_size_results", full.names = TRUE)
risultati_df <- bind_rows(lapply(res_files, readRDS))
risultati_df
saveRDS(risultati_df, file = "results/sample_size_results/sample_size_DEF.rds")





# # Table & Graph ICC & Sample size -----------------------------------------------------
# sample_size_db <- readRDS("risultati/sample_size_icc_ott.rds")

## SALVATI POSSIBILI GRAFICI NEL VECCHIO SCRIPT - RESULT_PAPER 2
# 




