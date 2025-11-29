library(simsurv)
library(survival)
library(tidyverse)
library(lme4)
library(coxme)
library(progressr)
library(R.utils) ## timeout
library(parallel)

#set.seed(123)

#source("funzioni paperSeptember2025.R") 
#source("funzioni paper - Ott 2025.R")
source("Paper Functions – Nov 2025.R")



# 2 PAPER: Impact of ICC on Power and Sample Size --------
### How different levels of ICC influence statistical power and required sample.

# 1) ICC E POWER -------------------------------------------------------------

### Apply surv_power_def to different scenarios
## This script assesses how statistical power changes with varying ICC values using the "surv_power_function", while keeping all other parameters fixed.
##
## Once the power is computed for each scenario at varying ICC values, the power with ICC = 0 is used as the reference.
##
## For ICC > 0, in the tibble of results returned by surv_power_function:
## 1. Compute the Design Effect (DE) and the DE-adjusted sample size.
## 2. Re-run the simulation with the adjusted sample size, keeping other parameters constant.
##
## Finally, compare the resulting power (after DE-based sample size adjustment) with the reference power (ICC = 0)
## By theory, the power obtained using the DE-adjusted sample size should be approximately the same as the power with ICC = 0.
##
## However, preliminary results suggest that the DE-adjusted sample size may overestimate the true sample size required to maintain equivalent power.


# Final scenario parameters
treat_effects <- c(- 0.2, - 0.5) ## HR: 0.8 medium, 0.6 high
lambdas <- c(0.115, 0.012) #means median: 6 months, 60 months
#num_pat_groups <- c(10, 20, 40, 60, 100)
sample_sizes <- c(100, 500, 2000, 6000)  
n_hosp_values <- c(5, 15, 30, 60)
icc_values <- c(0, 0.01, 0.02, 0.03, 0.05, 0.01, 0.4) 
cens_value <- c(0.5, 2, 5)
balancing_mode <- c(1, 2)

# ## test:
# icc_values <- c(0, 0.1) ## qui ne servono tanti - valuta
# lambdas <- 0.1 #c(0.01, 0.11)
# n_hosp_values <- c(10, 15)
# sample_sizes <- c(100, 1000)
# treat_effects <- 0.15
# cens_value <- c(0.5, 5)
# balancing_mode <- c(1, 2)


###### parallel computation + partial saving
unlink("results/power_results_individual", recursive = TRUE) ## delete, for testing
unlink("results/power_results_DE_individual", recursive = TRUE) ## delete, for testing

dir.create("results/power_results_individual", showWarnings = FALSE)
dir.create("results/power_results_DE_individual", showWarnings = FALSE)

scenarios <- expand.grid(
  pop_treat_effect = treat_effects,
  lambda = lambdas,
  icc = icc_values,
  num_hosp = n_hosp_values,
  #num_pat_group = num_pat_group_values,
  sample_size = sample_sizes,
  cens = cens_value,
  balancing= balancing_mode,
  stringsAsFactors = FALSE )

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
clusterExport(cl, c( "scenarios",
                     "simulate_survival_cohort_individual", ## or simulate_survival_cohort_hospital
                    "surv_power_function"))
clusterEvalQ(cl, {
  library(dplyr)
  library(simsurv)
  library(survival) })


# --- Funzione per calcolare scenario singolo ---
process_scenario <- function(i) {
  scen <- scenarios[i, ]
  filename <- paste0("results/power_results_individual/scenario_", i, ".rds")
  
 if (file.exists(filename)) return(readRDS(filename))  
  
  res <- tryCatch({
    surv_power_function(
      simula_coorte_fun = simulate_survival_cohort_individual,
      simula_args = list(
        num_hosp = scen$num_hosp,
        #num_pat_group = scen$num_pat_group,
        sample_size = scen$sample_size,
        balancing_mode= scen$balancing,
        pop_treat_effect = scen$pop_treat_effect,
        lambda = scen$lambda,
        icc = scen$icc,
        cens = scen$cens  ),
      nsim = 2000   )
  }, error = function(e) return(NULL))
  
  if (is.null(res)) {
    res <- cbind(scen, power = NA_real_)
    saveRDS(res, filename)
        return(res)  }   ## even if null, return the row

 res <- res %>%
    mutate(
      pop_treat_effect = scen$pop_treat_effect,
      lambda = scen$lambda,
      cens = scen$cens,
      icc = scen$icc,
      num_hosp = scen$num_hosp,
      num_pat_group_mean = scen$sample_size/scen$num_hosp,
      sample_size = scen$sample_size,
      balancing= scen$balancing,
      DE = 1 + ((scen$sample_size/scen$num_hosp) - 1) * scen$icc,
      sample_size_DE = ceiling(scen$sample_size * DE),
      m_DE = ceiling((scen$sample_size/scen$num_hosp) * DE) )
  
  saveRDS(res, filename)
  return(res) }


results_list <- parLapply(cl, 1:nrow(scenarios), process_scenario)

stopCluster(cl)  

# all results
res_fin <- list.files("results/power_results_individual", full.names = TRUE)
results_df <- bind_rows(lapply(res_fin, readRDS))


# --- DE scenario ---
scenarios_DE <- scenarios %>%
  filter(icc != 0) %>%
  left_join(
    results_df %>%
      select(icc, pop_treat_effect, lambda, num_hosp, cens, balancing, sample_size, sample_size_DE),
    by = c("icc", "pop_treat_effect", "lambda", "num_hosp", "cens", "balancing", "sample_size") )






cl <- makeCluster(n_cores)
clusterExport(cl, c("simulate_survival_cohort_individual", ## or simulate_survival_cohort_hospital
                    "surv_power_function", "scenarios_DE"))

clusterEvalQ(cl, {
  library(dplyr)
  library(simsurv)
  library(survival) })


process_DE <- function(i) {
  scen <- scenarios_DE[i, ]
  filename <- paste0("results/power_results_DE_individual/scenario_", i, ".rds")
  
  if(file.exists(filename)) return(NULL)
  
  res <- tryCatch({
    surv_power_function(
      simula_coorte_fun = simulate_survival_cohort_individual,
      simula_args = list(
        num_hosp = scen$num_hosp,
        sample_size =  scen$sample_size_DE, ##new sample size (correct)
        balancing_mode = scen$balancing,
        #num_pat_group = scen$num_pat_group,
        pop_treat_effect = scen$pop_treat_effect,
        icc = scen$icc,
        lambda = scen$lambda,
        cens = scen$cens ),
      nsim = 2000  )
  }, error = function(e) return(NULL))
  
  if(is.null(res)) return(NULL)
  
  res <- res %>%
    mutate(
      pop_treat_effect = scen$pop_treat_effect,
      lambda = scen$lambda,
      cens = scen$cens,
      num_hosp = scen$num_hosp,
      sample_size = scen$sample_size,
      #m_DE = scen$num_pat_group_values_corrected,
      sample_size_DE = scen$sample_size_DE,
      balancing= scen$balancing    )
  
  saveRDS(res, filename)
  return(res)}

results_DE_list <- parLapply(cl, 1:nrow(scenarios_DE), process_DE)

stopCluster(cl)

# --- Combine first tibble of results with results DE ---
res_fin2 <- list.files("results/power_results_DE_individual", full.names = TRUE)
results_df2 <- bind_rows(lapply(res_fin2, readRDS))


results_df2 <- results_df2 %>%
  select(-num_pat_group_mean, - sample_size_DE) %>%
  rename(power_DE = power)


results_def <- results_df %>%
  select(-nsim)%>%
 left_join(results_df2,
           by = c("pop_treat_effect", "lambda", "num_hosp", "icc", "cens", "balancing", "sample_size"))%>%
  mutate( design = "individual")

results_def <- results_def[, c("nsim", "icc", "lambda", "cens", "prop_cens.x", "prop_cens.y", "pop_treat_effect", "num_hosp", "balancing", "num_pat_group_mean", "sample_size", "power", "DE", "power_DE", "sample_size_DE", "m_DE")]

#saveRDS(results_df, file = "results/risultati_def_power_prova_grafici.rds")
#results_df<- readRDS("results/risultatioggi_prova_grafic.rds")


## Add two sample size calculations:
## - Schoenfeld and Freedman (do not account for ICC)
## - Xie correction (adjusts Schoenfeld/Freedman by design effect DE)
## Schoenfeld == Freedman when p = 0.5 (equal groups) and proportional hazards hold

results_def <- results_def %>%
  mutate(
    alpha = 0.05,
    z_alpha = qnorm(1 - alpha / 2), 
    p = 0.5,  # equal randomization by definition (code)
    prop_cens_mean = rowMeans(select(., prop_cens.x, prop_cens.y), na.rm = TRUE),
    HR = exp(pop_treat_effect),  # hazard ratio from treatment effect
    z_beta = qnorm(power), 
    E_schoenfeld = ((z_alpha + z_beta)^2) / (p * (1 - p) * (log(HR))^2),  ## Schoenfeld/Freedman: number of required events
    prop_eventi = 1 - prop_cens_mean,
    N_schoe_free = ceiling(E_schoenfeld / prop_eventi),
    ## Xie (2003): applies Schoenfeld/Freedman formula and multiplies by DE
    DE = 1 + (num_pat_group_mean - 1) * icc,  ## DE for cluster correction per Xie (2003)
    N_xie_freedman = ceiling(N_schoe_free * DE) ) %>% # corrected total sample size
  select(-z_beta, -HR, -prop_eventi, -E_schoenfeld, -alpha, -z_alpha, -p)


saveRDS(results_def, file = "results/power_results_individual/power_DEF_individual_prova_grafici.rds")

# ###METODO 3 * CV   - dal paper : è però per crt 9.	Practical considerations for sample size calculation for cluster randomized trials Clemence Leyrat 2024 
# ####cv fisso dal mio codice
# 
# de_cv <- 1+(0.29*0.29 + 1)*((mean_num_pat_group - 1)* icc) ### per cluster randomized trial
# n_cv_leyrat <- de_cv * n ###  ... (ma quale n) N START RAGIONARE SE TENERE 
# 



# #### HOSPITAL - different design ----------------------------------------


###### parallel computation + partial saving
# unlink("results/power_results_hospital", recursive = TRUE) ## delete, for testing
# unlink("results/power_results_DE_hospital", recursive = TRUE) ## delete, for testing

dir.create("results/power_results_hospital", showWarnings = FALSE)
dir.create("results/power_results_DE_hospital", showWarnings = FALSE)

scenarios <- expand.grid(
  pop_treat_effect = treat_effects,
  lambda = lambdas,
  icc = icc_values,
  num_hosp = n_hosp_values,
  #num_pat_group = num_pat_group_values,
  sample_size = sample_sizes,
  cens = cens_value,
  balancing= balancing_mode,
  stringsAsFactors = FALSE )

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
clusterExport(cl, c( "scenarios",
                     "simulate_survival_cohort_hospital", 
                     "surv_power_function"))
clusterEvalQ(cl, {
  library(dplyr)
  library(simsurv)
  library(survival) })


# --- Funzione per calcolare scenario singolo ---
process_scenario <- function(i) {
  scen <- scenarios[i, ]
  filename <- paste0("results/power_results_hospital/scenario_", i, ".rds")
  
  if (file.exists(filename)) return(readRDS(filename))  
  
  res <- tryCatch({
    surv_power_function(
      simula_coorte_fun = simulate_survival_cohort_hospital,
      simula_args = list(
        num_hosp = scen$num_hosp,
        #num_pat_group = scen$num_pat_group,
        sample_size = scen$sample_size,
        balancing_mode= scen$balancing,
        pop_treat_effect = scen$pop_treat_effect,
        lambda = scen$lambda,
        icc = scen$icc,
        cens = scen$cens  ),
      nsim = 2000   )
  }, error = function(e) return(NULL))
  
  if (is.null(res)) {
    res <- cbind(scen, power = NA_real_)
    saveRDS(res, filename)
    return(res)  }   ## even if null, return the row
  
  res <- res %>%
    mutate(
      pop_treat_effect = scen$pop_treat_effect,
      lambda = scen$lambda,
      cens = scen$cens,
      icc = scen$icc,
      num_hosp = scen$num_hosp,
      num_pat_group_mean = scen$sample_size/scen$num_hosp,
      sample_size = scen$sample_size,
      balancing= scen$balancing,
      DE = 1 + ((scen$sample_size/scen$num_hosp) - 1) * scen$icc,
      sample_size_DE = ceiling(scen$sample_size * DE),
      m_DE = ceiling((scen$sample_size/scen$num_hosp) * DE) )
  
  saveRDS(res, filename)
  return(res) }


results_list <- parLapply(cl, 1:nrow(scenarios), process_scenario)

stopCluster(cl)  

# all results
res_fin <- list.files("results/power_results_hospital", full.names = TRUE)
results_df <- bind_rows(lapply(res_fin, readRDS))


# --- DE scenario ---
scenarios_DE <- scenarios %>%
  filter(icc != 0) %>%
  left_join(
    results_df %>%
      select(icc, pop_treat_effect, lambda, num_hosp, cens, balancing, sample_size, sample_size_DE),
    by = c("icc", "pop_treat_effect","lambda", "num_hosp", "cens", "balancing", "sample_size") )



cl <- makeCluster(n_cores)
clusterExport(cl, c("simulate_survival_cohort_hospital", ## or simulate_survival_cohort_hospital
                    "surv_power_function", "scenarios_DE"))

clusterEvalQ(cl, {
  library(dplyr)
  library(simsurv)
  library(survival) })


process_DE <- function(i) {
  scen <- scenarios_DE[i, ]
  filename <- paste0("results/power_results_DE_hospital/scenario_", i, ".rds")
  
  if(file.exists(filename)) return(NULL)
  
  res <- tryCatch({
    surv_power_function(
      simula_coorte_fun = simulate_survival_cohort_hospital,
      simula_args = list(
        num_hosp = scen$num_hosp,
        sample_size =  scen$sample_size_DE, ##new sample size (correct)
        balancing_mode = scen$balancing,
        #num_pat_group = scen$num_pat_group,
        pop_treat_effect = scen$pop_treat_effect,
        icc = scen$icc,
        lambda = scen$lambda,
        cens = scen$cens ),
      nsim = 2000  )
  }, error = function(e) return(NULL))
  
  if(is.null(res)) return(NULL)
  
  res <- res %>%
    mutate(
      pop_treat_effect = scen$pop_treat_effect,
      lambda = scen$lambda,
      cens = scen$cens,
      num_hosp = scen$num_hosp,
      sample_size = scen$sample_size,
      #m_DE = scen$num_pat_group_values_corrected,
      sample_size_DE = scen$sample_size_DE,
      balancing= scen$balancing    )
  
  saveRDS(res, filename)
  return(res)}

results_DE_list <- parLapply(cl, 1:nrow(scenarios_DE), process_DE)

stopCluster(cl)

# --- Combine first tibble of results with results DE ---
res_fin2 <- list.files("results/power_results_DE_hospital", full.names = TRUE)
results_df2 <- bind_rows(lapply(res_fin2, readRDS))

results_df2 <- results_df2 %>%
  select(-num_pat_group_mean, - sample_size_DE) %>%
  rename(power_DE = power)


results_def <- results_df %>%
  select(-nsim)%>%
  left_join(results_df2,
            by = c("pop_treat_effect", "lambda", "num_hosp", "icc", "cens", "balancing", "sample_size"))%>%
  mutate( design = "hospital")
  

results_def <- results_def[, c("nsim", "icc", "lambda", "cens", "prop_cens.x", "prop_cens.y", "pop_treat_effect", "num_hosp", "balancing", "num_pat_group_mean", "sample_size", "power", "DE", "power_DE", "sample_size_DE", "m_DE")]

#saveRDS(results_df, file = "results/risultati_def_power.rds")
#results_df<- readRDS("results/risultatioggi.rds")


## Add two sample size calculations:
## - Schoenfeld and Freedman (do not account for ICC)
## - Xie correction (adjusts Schoenfeld/Freedman by design effect DE)
## Schoenfeld == Freedman when p = 0.5 (equal groups) and proportional hazards hold

results_def <- results_def %>%
  mutate(
    alpha = 0.05,
    z_alpha = qnorm(1 - alpha / 2), 
    p = 0.5,  # equal randomization by definition (code)
    prop_cens_mean = rowMeans(select(., prop_cens.x, prop_cens.y), na.rm = TRUE),
    HR = exp(pop_treat_effect),  # hazard ratio from treatment effect
    z_beta = qnorm(power), 
    E_schoenfeld = ((z_alpha + z_beta)^2) / (p * (1 - p) * (log(HR))^2),  ## Schoenfeld/Freedman: number of required events
    prop_eventi = 1 - prop_cens_mean,
    N_schoe_free = ceiling(E_schoenfeld / prop_eventi),
    ## Xie (2003): applies Schoenfeld/Freedman formula and multiplies by DE
    DE = 1 + (num_pat_group_mean - 1) * icc,  ## DE for cluster correction per Xie (2003)
    N_xie_freedman = ceiling(N_schoe_free * DE) ) %>% # corrected total sample size
  select(-z_beta, -HR, -prop_eventi, -E_schoenfeld, -alpha, -z_alpha, -p)


saveRDS(results_def, file = "results/power_results_hospital/power_DEF_hospital.rds")

# ###METODO 3 * CV   - dal paper : è però per crt 9.	Practical considerations for sample size calculation for cluster randomized trials Clemence Leyrat 2024 
# ####cv fisso dal mio codice
# 
# de_cv <- 1+(0.29*0.29 + 1)*((mean_num_pat_group - 1)* icc) ### per cluster randomized trial
# n_cv_leyrat <- de_cv * n ###  ... (ma quale n) N START RAGIONARE SE TENERE 
# 







# Table & GraPHS: ICC & Power -----------------------------------------------------

## SALVATI POSSIBILI GRAFICI NEL VECCHIO SCRIPT - RESULT_PAPER 2
# 
# 
