
### FUNCTIONS IN THIS FILE:
#
# 1) Simulating data 
# 2) Assess ICC
# 3) Bootstrap
# 4) Power
# 5) Sample Size



## TO DO / Notes:
##
## 1. Review the individual formulas used for ICC estimation (function icc_estimation)
##    – ensure they match the methods reported in the paper/draft.
##    – particular attention to Oliveira (2016), and optionally consider including Jung (2007) if you think it is appropriate.

## 2. Check the validity of the simulation code:
##    – does it correctly simulate the ICC target? (target_ICC) - Are we correctly simulating data with the intended target ICC?

## 3. Check on the general structure and logic of the code

## 4. Proposed (NEW) addition: simulate scenarios where treatment is assigned at the cluster level instead of individually, as an additional set of scenarios To explore. 
##    – does this approach make sense and provide added value?


library(simsurv)
library(survival)
library(tidyverse)
library(lme4)
library(coxme)
library(R.utils) ## timeout


# 1 A) Simulate survival cohort with NO interaction (hospital × treat) --------
#
# Input: ICC
# 
# This function simulates a cohort of survival data with patients clustered 
# within hospitals. The user provides:
#   - num_hosp          : number of hospitals (clusters)
#   - num_pat_group     : number of patients per hospital
#   - icc               : intra-class correlation to define cluster effect
#   - pop_treat_effect  : fixed treatment effect
#   - lambdas           : represents how frequently events occur
#   - gammas            : shape parameter (fixed to 1 for exponential baseline)
#   - maxt              : fixed 2*median
#
# Output: a data frame containing:
#   - id                : patient identifier
#   - hospital          : hospital/cluster identifier
#   - treat             : treatment assignment (0/1)
#   - eventtime         : observed survival/censoring time
#   - status            : event indicator (1=event, 0=censored)
#   - sigma_hosp        : standard deviation of hospital random effect 
#   - cluster_eff       : hospital-specific random effect (frailty)
#
# Notes:
#   - Data simulation inspired by: Brilleman SL. "Simulating Survival Data Using the simsurv R Package"
#     we use a weibull distribution (with gamma == 1, exponential) and a gaussian fraility term to model the cluster effect
#   - The hazard function assumes NO interaction between hospital and treatment:
#                  h_ij(t) = h0(t) * exp(cluster_eff_j + treat_ij * beta)
#     Cluster_eff is distributed as N(0, sigma_hosp), which corresponds to the between-hospital variance.
#     ICC: The within-cluster variance (residual variance) for the model is π^2/(6) - approximation
#     the treatment effect (beta) is the same across all hospitals.
#   - Sigma_hosp is derived from ICC using assumption of log-normal frailty distribution
#   - al momento non considerare: "funzioni altre" allow  1) specifying sigma_hosp for simulating  
#                             2) a model with treatment × hospital interaction: h_ij(t) = h0(t) * exp(treat_ij * (beta + u_j))    interaction model
#   - The maximum follow-up time (maxt) is chosen empirically (2× median expect) according to Klein & Moeschberger ("Survival Analysis: Techniques for Censored and Truncated Data" - 2003)

## NEW: vary administrative censoring; scenarios: cens = 0.5 → high (~60-70%), cens = 2 → medium (~20-30%), cens = 5 → low (<10%)
## NEW 2!! Added the option to create unbalanced groups: balancing_mode = 1: balanced groups, = 2`: unbalanced groups with SD = 0.5 and CV (coefficient of variation of group sizes) ≈ 0.29


## Attention: In this formula, patients within each group receive treatment with a 50/50 split

simulate_survival_cohort_individual <- function(num_hosp, 
                                                        sample_size,
                                                        balancing_mode = 1,  # 1 = balanced, 2 = unbalanced
                                                        icc, 
                                                        pop_treat_effect, 
                                                        lambda, 
                                                        gammas = 1,
                                                        cens = 2) {

  if (balancing_mode == 1) {
  
    num_pat_group <- rep(floor(sample_size / num_hosp), num_hosp)   # balanced
    resto <- sample_size - sum(num_pat_group)
    if (resto > 0) num_pat_group[1:resto] <- num_pat_group[1:resto] + 1     } 
  
  if (balancing_mode == 2) {

    mean_size <- sample_size / num_hosp # unbalanced groups while keeping the total sample size constant

    pattern_name <- paste0("fixed_var_factor_", num_hosp) # unbalancing with a fixed coefficient of variation for each number of hospitals (num_hosp)
    
    if (!exists(pattern_name, envir = .GlobalEnv)) {
      set.seed(num_hosp)  # ensures the same pattern for the same number of hospitals
      assign(pattern_name, runif(num_hosp, 0.5, 1.5), envir = .GlobalEnv)    }   ## in this way the variability of the groups is about 30%. Is that okay? Should it be increased?
    
    var_factor <- get(pattern_name, envir = .GlobalEnv)
    
    num_pat_group <- mean_size * var_factor
    num_pat_group <- num_pat_group / sum(num_pat_group) * sample_size
    num_pat_group <- round(num_pat_group)
    
    num_pat_group[num_pat_group <= 0] <- 1  # no hospital with 0 patients 
     
   diff <- sample_size - sum(num_pat_group) # corrects to match exact sample_size 
    if (diff != 0) {
      adjust_index <- sample(1:num_hosp, 1)
      num_pat_group[adjust_index] <- num_pat_group[adjust_index] + diff    }  }
  
  
  tot_patients <- sum(num_pat_group)
  var_resid <- (pi^2) / 6 ### residual variance assumption (see Paper 1)
  sigma2_hosp <- abs((icc / (1 - icc)) * var_resid)
  sigma_hosp <- sqrt(sigma2_hosp)
  
  hospital <- rep(1:num_hosp, times = num_pat_group)
  cov <- data.frame(
    id = 1:tot_patients,
    hospital = hospital,
    treat = rbinom(tot_patients, 1, 0.5)  )
  
  cluster_effect <- rnorm(num_hosp, mean = 0, sd = sigma_hosp)
  cov$cluster_eff <- cluster_effect[cov$hospital]
  
  beta <- c(treat = pop_treat_effect, cluster_eff = 1)
  
  dati_sopravvivenza <- simsurv(
    dist = "weibull",
    lambdas = lambda,
    gammas = gammas,
    x = cov[, c("treat", "cluster_eff")],
    betas = beta,
    maxt = cens * (log(2) / lambda)  )
  
  cohort <- merge(cov, dati_sopravvivenza, by = "id")
  cohort$sigma_hosp <- as.numeric(sigma_hosp[1])
  cohort$icc <- as.numeric(icc[1])
  cohort$cens_prop <- mean(cohort$status == 0)
  
  
  cohort[] <- lapply(cohort, function(x) if (is.numeric(x)) round(x, 3) else x)
  
  return(cohort) }


# cohort <- simulate_survival_cohort_individual(
#   num_hosp = 15,            # number of hospitals
#   sample_size = 100,        # total number of patients
#   balancing_mode = 2,       # 1 = balanced, 2 = unbalanced
#   icc = 0.01,               # intraclass correlation coefficient
#   pop_treat_effect = 0.15,  # treatment effect
#   lambda = 0.1,             # Weibull scale parameter
#   gammas = 1,               # shape parameter
#   cens = 0.5 )              # censoring level
# 
# 
# cohort %>%
#   gtsummary::tbl_summary(by = treat)



## 1 B)
# ATTENTION!! In this formula, each cluster receives only one treatment: 50% of hospitals get treatment = 0, 50% get treatment = 1

simulate_survival_cohort_hospital <- function(num_hosp, 
                                                                            sample_size,
                                                                            balancing_mode = 1,  # 1 = bilanciato, 2 = sbilanciato
                                                                            icc, 
                                                                            pop_treat_effect, 
                                                                            lambda, 
                                                                            gammas = 1,
                                                                            cens = 2) {
  
  if (balancing_mode == 1) { # balanced
    num_pat_group <- rep(floor(sample_size / num_hosp), num_hosp)
    resto <- sample_size - sum(num_pat_group)
    if (resto > 0) num_pat_group[1:resto] <- num_pat_group[1:resto] + 1     } 
  
  if (balancing_mode == 2) {
    mean_size <- sample_size / num_hosp # unbalanced groups while keeping the total sample size constant
    
    pattern_name <- paste0("fixed_var_factor_", num_hosp) # unbalancing with a fixed coefficient of variation for each number of hospitals (num_hosp)
    
    if (!exists(pattern_name, envir = .GlobalEnv)) {
      set.seed(num_hosp)  # garantisce stesso pattern per stesso numero di ospedali
      assign(pattern_name, runif(num_hosp, 0.5, 1.5), envir = .GlobalEnv)    }
    
    var_factor <- get(pattern_name, envir = .GlobalEnv)
    
    num_pat_group <- mean_size * var_factor
    num_pat_group <- num_pat_group / sum(num_pat_group) * sample_size
    num_pat_group <- round(num_pat_group)
    
     num_pat_group[num_pat_group <= 0] <- 1 # no hospital with 0 patients
    
    diff <- sample_size - sum(num_pat_group) # corrects to match exact sample_size
    if (diff != 0) {
      adjust_index <- sample(1:num_hosp, 1)
      num_pat_group[adjust_index] <- num_pat_group[adjust_index] + diff    }  }
  
  
  tot_patients <- sum(num_pat_group)
  var_resid <- (pi^2) / 6
  sigma2_hosp <- abs((icc / (1 - icc)) * var_resid)
  sigma_hosp <- sqrt(sigma2_hosp)
  
  hospital <- rep(1:num_hosp, times = num_pat_group)
  
  hospital_treat <- rbinom(num_hosp, 1, 0.5)   ## section differing from the first function
  treat <- hospital_treat[hospital]
  
  cov <- data.frame(
    id = 1:tot_patients,
    hospital = hospital,
    treat = treat   )
  
  
  cluster_effect <- rnorm(num_hosp, mean = 0, sd = sigma_hosp)
  cov$cluster_eff <- cluster_effect[cov$hospital]
  
  beta <- c(treat = pop_treat_effect, cluster_eff = 1)
  
  dati_sopravvivenza <- simsurv(
    dist = "weibull",
    lambdas = lambda,
    gammas = gammas,
    x = cov[, c("treat", "cluster_eff")],
    betas = beta,
    maxt = cens * (log(2) / lambda)  )
  
  cohort <- merge(cov, dati_sopravvivenza, by = "id")
  cohort$sigma_hosp <- as.numeric(sigma_hosp[1])
  cohort$icc <- as.numeric(icc[1])
  cohort$cens_prop <- mean(cohort$status == 0)
  
  
  cohort[] <- lapply(cohort, function(x) if (is.numeric(x)) round(x, 3) else x)
  
  return(cohort) }

# 
# cohort <- simulate_survival_cohort_hospital(
#   num_hosp = 15,             # numero ospedali
#   sample_size = 100,        # numero totale pazienti
#   balancing_mode = 2,      # 1 = bilanciato, 2 = sbilanciato
#   icc = 0.0,               # intraclass correlation coefficient
#   pop_treat_effect = 0.15,  # effetto del trattamento
#   lambda = 0.1,             # parametro di scala della Weibull
#   gammas = 1,               # parametro di forma
#   cens = 0.5    )             # censura (max follow-up time)
# 
# cohort%>%
#   gtsummary::tbl_summary(by = treat)
# 



# 2) Estimated ICC ------------------------------------------------------------
# ICC is estimated using methods from the literature, applied to a simulated cohort.
# Each method is described in the table provided in draft folfer - paper 1, with the reference paper listed  next to the code. All papers are available in the shared Zotero library.


icc_estimation <- function(data) {
  
  # Method 1 – Weibull - from:  Williams (1995) Crowder & Hand (1990)
  data$logT <- log(data$eventtime) # log-tempi di sopravvivenza
  fit <- lmer(logT ~ 1 +  (1 | hospital), data = data) # fit modello misto: random intercept per cluster
  var_cluster <- as.numeric(VarCorr(fit)$hospital[1])# varianza del cluster
 
  var_resid_theoretical <- pi^2 / (6 * 1) # gamma = 1   # varianza residua
  icc_weibull <- var_cluster / (var_cluster + var_resid_theoretical)
  
  
  # Method 2 – CoxME Gaussian frailty - from: (McCune, 2025) 
  mod_coxme <- tryCatch({
    coxme(Surv(eventtime, status) ~ 1 + (1|hospital), data = data)
  }, error = function(e) NULL)
  
  icc_coxme <- if (!is.null(mod_coxme)) {
    var_coxme <- VarCorr(mod_coxme)$hospital[[1]]
    var_coxme / (var_coxme + pi^2 / 6)
  } else { NA }
  
  
  # Method 3 – Cox with gamma frailty (non-parametric ICC) - from: (McCune, 2025) 
  cox_frailty_gamma <- tryCatch({
    coxph(Surv(eventtime, status) ~ 1 + frailty(hospital, distribution = "gamma"), data = data)
  }, error = function(e) NULL)
  
  var_gamma <- if (!is.null(cox_frailty_gamma)) {
    as.numeric(cox_frailty_gamma$history$`frailty(hospital, distribution = "gamma")`[1])
  } else { NA }
  
  if (!is.na(var_gamma)) {
    g <- function(w, k, s, sigma2) { -k * w + exp(w) * s + w^2 / (2 * sigma2) }
    g2 <- function(w, k, s, sigma2) { exp(w) * s + 1 / sigma2 }
    
    Lapl <- Vectorize(function(s, k, sigma2) {
      wTilde <- optimize(f = g, c(-1e10, 1e10), maximum = FALSE,
                         k = k, s = s, sigma2 = sigma2)$minimum
      (-1)^k * exp(-g(w = wTilde, k = k, s = s, sigma2 = sigma2)) /
        sqrt(sigma2 * g2(w = wTilde, k = k, s = s, sigma2 = sigma2))
    }, "s")
    
    fr.lognormal <- function(k, s, sigma2) {
      intTau <- Vectorize(function(x, intTau.sigma2 = sigma2) {
        x * Lapl(s = x, k = 0, sigma2 = intTau.sigma2) *
          Lapl(s = x, k = 2, sigma2 = intTau.sigma2)
      }, "x")
      4 * integrate(f = intTau, lower = 0, upper = Inf,
                    intTau.sigma2 = sigma2)$value - 1 }
    
    icc_gamma_np <- fr.lognormal(k = 1, s = 1, sigma2 = var_gamma)
  } else { icc_gamma_np <- NA }
  
 
  # - Method 4 – Binomial GLMM (cloglog/logit link) on discretized time - from: (McCune, 2025 + Lam & Ip 2003) 
  ### We choose the cloglog link, which is more suitable for survival data
  
  # Time discretized into 4 quantile-based intervals for generalization across datasets
  q <- quantile(data$eventtime, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  q <- unique(q) ##remove any duplicates
  tmin <- min(data$eventtime, na.rm = TRUE)
  q <- q[q > tmin]
  
  # discretized time db
  discretized_db <- survSplit(
    Surv(eventtime, status) ~ hospital,
    data = data,
    cut = q,
    start = "tstart",
    end = "tstop"  )
  
 discretized_db$interval <- case_when(
    discretized_db$tstart == 0                   ~ 1,
    discretized_db$tstart == q[1]               ~ 2,
    length(q) > 1 & discretized_db$tstart == q[2] ~ 3,
    length(q) > 2 & discretized_db$tstart == q[3] ~ 4  )
  
  # GLMM with link cloglog
  mod_glmm <- tryCatch({
    glmer(
      status ~ as.factor(interval) + (1|hospital),
      data = discretized_db,
      family = binomial(link = "cloglog"), 
      nAGQ = 5    ) ## The paper used 7, but 5 gives slightly more variability while being much faster
  }, error = function(e) NULL)
  
  # ICC
  icc_glmm <- if (!is.null(mod_glmm) && !isSingular(mod_glmm)) {
    var_glmm <- mod_glmm@theta^2
    var_glmm / (var_glmm + pi^2 /6)
  } else { NA }
  
  
  
  # Method 5 – CoxME on discretized time - from: McCune, 2025
  mod_cox_exploded <- tryCatch({
    coxme(Surv(tstart, tstop, status) ~ 1 + (1|hospital), data = discretized_db)
  }, error = function(e) NULL)
  
  icc_cox_exploded <- if (!is.null(mod_cox_exploded)) {
    var_exploded <- VarCorr(mod_cox_exploded)$hospital[[1]]
    var_exploded / (var_exploded + pi^2 / 6)  ## approximating the residual standard deviation (sigma)
  } else { NA }
  
  
  # Method 6 – ICC with censoring indicators - from: Kalia, Klar & Donner (2016)
  grand_mean <- mean(data$status, na.rm = TRUE) # 2. overall mean 
  
  cluster_n <- tapply(data$status, data$hospital, length)   # numero di pazienti per cluster
  cluster_mean <- tapply(data$status, data$hospital, mean)  # proporzione di eventi nel cluster
  
  
  r <- length(cluster_n)   # numero cluster
  MSB <- sum(cluster_n * (cluster_mean - grand_mean)^2) / (r - 1) # MSB (Between-cluster Mean Square)
  
  cluster_mean_ind <- cluster_mean[data$hospital]
  MSW <- sum((data$status - cluster_mean_ind)^2) / (nrow(data) - r) # MSW (Within-cluster Mean Square)
  
  m_bar <- mean(cluster_n)
  icc_censoring <- (MSB - MSW) / (MSB + (m_bar - 1) * MSW)

  
  # Method 7 – ICC with observed event times - from: Kalia, Klar & Donner (2016)

  obs_data <- data[data$status == 1, ] # consider only observed events (exclude censored)
  obs_data$hospital <- as.factor(obs_data$hospital)
  obs_data$hospital <- droplevels(obs_data$hospital)
  
  # cluster-level mean times and overall mean
  cluster_means <- tapply(obs_data$eventtime, obs_data$hospital, mean)
  overall_mean <- mean(obs_data$eventtime)
  
  mi <- table(obs_data$hospital)    # events per hospital
  r <- length(mi)                   # number of hospital
  N_star <- sum(mi)                 # total number of events
  
  m_o <- (1 / (r - 1)) * (N_star - sum(mi^2) / N_star)
  
  MSB <- (1 / (r - 1)) * sum(mi * (cluster_means - overall_mean)^2) # between hospital variance (MSB)
  
  
  MSW <- (1 / (N_star - r)) * sum(tapply(obs_data$eventtime, obs_data$hospital,#within hospital variance (MSW)
                                         function(x) sum((x - mean(x))^2)))
  icc_event <- (MSB - MSW) / (MSB + (m_o - 1) * MSW)

  
  
  #### - Method 9 - Martingale residual-based ICC (Cox, non-parametric) - from: Gagnon, 2004
  
  cox_fit <- coxph(Surv(eventtime, status) ~ 1, data = data) # Fit Cox model ignoring clustering to get martingale residuals
  data$mart_resid <- residuals(cox_fit, type = "martingale")
  
  cluster_means <- data %>% # Compute cluster means of residuals
    group_by(hospital) %>%
    summarise(mean_resid = mean(mart_resid, na.rm = TRUE),
              n = n(), .groups = "drop")
  
  grand_mean <- mean(data$mart_resid, na.rm = TRUE)
  MSB <- sum(cluster_means$n * (cluster_means$mean_resid - grand_mean)^2) /
    (nrow(cluster_means) - 1)
  
  data <- left_join(data, cluster_means, by = "hospital")
  MSW <- sum((data$mart_resid - data$mean_resid)^2) /
    (nrow(data) - nrow(cluster_means))
  
  icc_martingale <- (MSB - MSW) / (MSB + (mean(cluster_means$n) - 1) * MSW)
  
  
  ## - Method 9 – Kaplan-Meier ICC - from: Jung, 2007 (Ying, 1994) 
  ## As of today, I have not been able to implement this method.
  
  ## QUESTIONS: Is it worth considering whether to try to add it or whether the methods already implemented are sufficient?
  
  ## Reference: Jung (2007), nonparametric approach using marginal and bivariate survival function estimators,
  ## extending Kaplan–Meier to clustered data. Allows ICC estimation without frailty or parametric assumptions,
  ## but can be unstable with few clusters or heavy censoring.
  
  
  ## - Method 8 PROPORTION OF EVENT - from: Xie & Waksman. (2003)
  p_hat <- mean(data$status, na.rm = TRUE) # overall event proportion
  
  data$dstatus <- data$status - p_hat # deviation from mean for each subject
  
  agg <- by(data$dstatus, data$hospital, function(x) {
    s <- sum(x)
    ssq <- sum(x^2)
    list(K = length(x), pair_sum = (s^2 - ssq))})
  
  agg_list <- lapply(agg, unclass)
  Ks <- vapply(agg_list, `[[`, numeric(1), "K")
  pair_sums <- vapply(agg_list, `[[`, numeric(1), "pair_sum")
  
  numerator <- sum(pair_sums)
  denom_pairs <- sum(Ks * (Ks - 1))
  denominator <- p_hat * (1 - p_hat) * denom_pairs
  
  icc_proportion <- numerator / denominator
  

  ## Method 10: Weibull model with combined frailty (Normal + Gamma) - from Oliveira et al. (2016) 
  ### I AM UNSURE IF I HAVE IMPLEMENTED THE METHODS CORRECTLY!!
  
  icc_weibull_combined <- tryCatch({
    coxme(Surv(eventtime, status) ~ 1 + (1|hospital), data = data)
  }, error = function(e) NULL)
  
  icc_weibull_combined <- if (!is.null(icc_weibull_combined)) {
    D <- VarCorr(icc_weibull_combined)$hospital[[1]]
    D / (D + 1)  #### overdispersion == 1 in caso fraility gamma esponenziale
  } else { NA }
  
  
  icc_weibull <- ifelse(exists("icc_weibull"), icc_weibull, NA)
  icc_coxme <- ifelse(exists("icc_coxme"), icc_coxme, NA)
  icc_gamma_np <- ifelse(exists("icc_gamma_np"), icc_gamma_np, NA)
  icc_glmm <- ifelse(exists("icc_glmm"), icc_glmm, NA)
  icc_cox_exploded <- ifelse(exists("icc_cox_exploded"), icc_cox_exploded, NA)
  icc_censoring <- ifelse(exists("icc_censoring"), icc_censoring, NA)
  icc_event <- ifelse(exists("icc_event"), icc_event, NA)
  icc_martingale <- ifelse(exists("icc_martingale"), icc_martingale, NA)
  icc_proportion <- ifelse(exists("icc_proportion"), icc_proportion, NA)
  icc_weibull_combined <- ifelse(exists("icc_weibull_combined"), icc_weibull_combined, NA)
  
  icc_theoretical <- if("icc" %in% names(data)) {
    unique(data$icc)[1]   
  } else {
    NA   }
  
  # FINAL TIBBLE 
  risultati <- data.frame(
    Method = c(
      "Analytical Log-Weibull",       # 1
      "CoxME Gaussian frailty",       # 2
      "Cox gamma frailty (NP)",       # 3
      "GLMM logit discretized",       # 4
      "CoxME discretized time",       # 5
      "Censoring indicators",         # 6
      "Observed event times",         # 7
      "Martingale",                   # 8
      #"KM" ,                         # 9
      "Weibull combined model",       #10
      "Proportion of event",
      "Start ICC - 0"    ),
    ICC = round(c(
      icc_weibull,
      icc_coxme,
      icc_gamma_np,
      icc_glmm,
      icc_cox_exploded,
      icc_censoring,
      icc_event,
      icc_martingale,
      #icc_km,
      icc_weibull_combined,
      icc_proportion, 
      icc_theoretical   ), 3)  )
return(risultati)}


# 3) Bootstrap ICC ------------------------------------------------------
# This bootstrap samples entire clusters (e.g., hospitals) with replacement.
# For each selected cluster, all subjects belonging to that cluster are included.
# The ICC is calculated on each bootstrap dataset.

### THIS METHOD IS USEFUL TO ASSESS THE ICC AND HIS CONFIDENCE INTERVAL WHEN WE HAVE A SINGLE REAL DATASET.
### WHEN USING SIMULATION, WE INSTEAD APPLY A MONTE CARLO APPROACH (CREATING MULTIPLE DATASETS WITH FIXED PARAMETERS) TO ESTIMATE ICC AND ITS CI.


bootstrap_icc <- function(data, B = 100, cluster_var = "hospital") {
  
  boot_results <- vector("list", B)
  
  # Cluster
  clusters <- unique(data[[cluster_var]])
  
  for (b in 1:B) {
    # --- Resampling for cluster ---
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
    
    # Assign unique cluster IDs to duplicated clusters to preserve all rows
    data_boot <- do.call(rbind, lapply(seq_along(sampled_clusters), function(i) {
      cl <- sampled_clusters[i]
      subset_data <- data[data[[cluster_var]] == cl, ]
      subset_data[[cluster_var]] <- i
      subset_data   }))
    
    data_boot[[cluster_var]] <- factor(data_boot[[cluster_var]])
    
    # --- Asses ICC on dataset bootstrap ---
    res <- tryCatch({
      icc_estimation(data_boot)$ICC
    }, error = function(e) rep(NA, ncol(icc_estimation(data)$ICC)))
    
    boot_results[[b]] <- res   }
  
  # Matrix
  boot_mat <- do.call(rbind, boot_results)
  colnames(boot_mat) <- icc_estimation(data)$Method
  
  # Mean-Median & CI
  ## QUESTION: Mean & SD vs. Quantile? Likely similar due to many repetitions. Which do you prefer?
  
  boot_summary <- data.frame(
    Method = colnames(boot_mat),
    Mean = apply(boot_mat, 2, mean, na.rm = TRUE),
    SD   = apply(boot_mat, 2, sd, na.rm = TRUE),
    median = apply(boot_mat, 2, quantile, probs = 0.5, na.rm = TRUE),
    CI_low  = apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    CI_high = apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE) )
  
  return(list(
    summary = boot_summary,
    replicates = boot_mat  ))}

#  
# #examples
# #set.seed(123)
# cohort <- simulate_survival_cohort_individual(icc= 0.05, pop_treat_effect = 0.8, sample_size = 1000, num_hosp = 10, lambda = 0.01)
# result <- icc_estimation(cohort)
# result
# boot_out <- bootstrap_icc(cohort, B = 10)
# print(boot_out$summary)
# 
# 
# 


# 4) Power & ICC -----------------------------------------------------------

# FUNCTION: surv_power_function
# 
# PURPOSE:  Estimate statistical power for a survival study using Monte Carlo simulations
#
# HOW IT WORKS:
# 1. Generates 'nsim' simulated cohorts using the provided simulation function
# 2. For each cohort, fits Cox proportional hazards models (with clustering at the hospital level (coxph + cluster(hospital))
# 3. Records p-values for the treatment effect and counts how many simulations are statistically significant (p < 0.05)
# 4. Calculates: power: proportion significant in the clustered model
# 5. Returns a tibble with the estimated powers and input parameters

surv_power_function <- function(simula_coorte_fun, simula_args = list(),
                                    nsim = 100) {
  
  p_values_icc <- numeric(nsim)
  significant <- logical(nsim)
  prop_cens <- numeric(nsim)  # % of censoring 
  
  for (i in 1:nsim) {
    cohort <- do.call(simula_coorte_fun, simula_args)# simulate cohorts
    
    if (i == 1) {
      sigma_hosp <- unique(cohort$sigma_hosp)
      icc <- simula_args$icc
      num_hosp <- length(unique(cohort$hospital))
      num_pat_group_mean <- nrow(cohort) / num_hosp
      sample_size <- nrow(cohort)}
    
    prop_cens[i] <- mean(cohort$status == 0, na.rm = TRUE) # % of censoring
    
    mod_icc <- coxph(Surv(eventtime, status) ~ treat + cluster(hospital), data = cohort) ###cox model with cluster corretion
    
    p_val_icc <- summary(mod_icc)$coefficients["treat", "Pr(>|z|)"]   # p-value
    p_values_icc[i] <- p_val_icc
    significant[i] <- p_val_icc < 0.05 }
  
  power <- mean(significant)  # power estimation
  prop_cens_mean <- mean(prop_cens, na.rm = TRUE)
  
  results <- tibble(
    nsim = nsim,
    icc = icc,
    num_hosp = num_hosp,
    num_pat_group_mean = num_pat_group_mean,
    sample_size = sample_size,
    power = power,
    prop_cens = prop_cens_mean)
  
  return(results)}

 
# res <- surv_power_function(
#   simula_coorte_fun = simulate_survival_cohort_individual,  #### or simulate_survival_cohort_hospital
#   simula_args = list( cens = 0.5,
#     num_hosp = 10,
#     sample_size = 1000,
#     icc = 0.0,
#     pop_treat_effect = 0.15,
#     lambda = 0.1,
#     gammas = 1,
#     balancing_mode = 1  ),
#   nsim = 10)
# 
# print(res)





# 5) ICC & Sample Size ---------------------------------------------------

# FUNCTION: sample_size_icc
#
# PURPOSE:  Estimate the minimum number of patients per group required to achieve a target statistical power (≥ 0.8)
#           for a clustered survival study, accounting for intraclass correlation (ICC) among hospitals
#
# HOW IT WORKS:
# 1. Starts with an initial number of patients per group
# 2. Generates a simulated cohort 
# 3. Runs 'nsim' Monte Carlo simulations, fitting a Cox proportional hazards model with clustering at the hospital level, 
#    and calculates the proportion of simulations where the treatment effect is statistically significant (p < 0.05): this is the estimated power!
# 4. Increments the number of patients per group until the target power is reached
# 5. Returns a tibble (...)

sample_size_icc <- function(icc, num_hosp, pop_treat_effect, nsim = 10,
                            max_pat_group = 500, timeout_sec = 300, ## input to avoid infinite loops
                            lambda = 0.1, gammas = 1, n_init, balancing_mode) {
  
  npt <- n_init  # starting number of patients per group
  found_target <- FALSE
  candidate_npt <- NULL
  candidate_result <- NULL
  
  while (npt <= max_pat_group) {
    message("Testing NPT = ", npt)
    
    result <- tryCatch({
      withTimeout({
        significant <- logical(nsim)
        
        for (i in 1:nsim) {
         
          cohort <- simulate_survival_cohort_individual( # simulated cohort #### or simulate_survival_cohort_hospital
            num_hosp = num_hosp,
            sample_size = npt*num_hosp,
            icc = icc,
            pop_treat_effect = pop_treat_effect,
            lambda = lambda,
            gammas = gammas,
            balancing_mode = balancing_mode)
          
          mod_icc <- coxph(Surv(eventtime, status) ~ treat + cluster(hospital), data = cohort) ## cox model accounting for clustering (as power function above)
          p_val <- summary(mod_icc)$coefficients["treat", "Pr(>|z|)"]
          significant[i] <- p_val < 0.05}
        power <- mean(significant)
        
        tibble(
          nsim = nsim,
          icc = icc,
          num_hosp = num_hosp,
          lambda = lambda,
          pop_treat_effect = pop_treat_effect,
          balancing_mode = balancing_mode,
          num_pat_group_mean = npt,
          sample_size = npt * num_hosp,
          power = power  )  }, timeout = timeout_sec, onTimeout = "silent")
    }, error = function(e) { message("Error: ", conditionMessage(e)); NULL })
    
    valid <- !is.null(result) && !is.null(result$power) && !is.na(result$power)
    
    if (valid && result$power >= 0.8) {
      found_target <- TRUE
      candidate_npt <- npt
      candidate_result <- result
      break  }
    
    npt <- npt + 1  }
  
  if (!found_target) {
    return(tibble(
      num_hosp = num_hosp,
      num_pat_group_mean = NA_integer_,
      nsim = nsim,
      power = NA_real_,
      sample_size = NA_real_,
      icc = icc,
      pop_treat_effect = pop_treat_effect,
      lambda = lambda,
      balancing_mode = balancing_mode)) }
  
  candidate_result }


# #example:
# res <- sample_size_icc(
#   icc = 0.1,
#    num_hosp = 10,
#    pop_treat_effect = 0.5,
#    nsim = 10,             # test veloce
#    lambda = 0.1,
#    n_init = 5,
#   balancing_mode = 2)
# 
# print(res)



