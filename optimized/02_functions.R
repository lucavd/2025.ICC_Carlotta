# =============================================================================
# 02_functions.R - Funzioni di simulazione + helper per chunking/checkpoint
# =============================================================================

library(simsurv)
library(survival)
library(tidyverse)
library(lme4)
library(coxme)

# =============================================================================
# HELPER FUNCTIONS PER CHUNKING E CHECKPOINT
# =============================================================================

#' Esegue una funzione con checkpoint per chunk
#' 
#' @param scenario_id ID univoco dello scenario (usato per nome file)
#' @param n_rep Numero totale di repliche
#' @param chunk_size Dimensione di ogni chunk
#' @param dir_out Directory di output
#' @param rep_fun Funzione da eseguire per ogni replica (deve restituire un data.frame/tibble)
#' @param ... Argomenti passati a rep_fun
#' @return tibble con risultati aggregati
run_with_chunks <- function(scenario_id, n_rep, chunk_size, dir_out, rep_fun, ...) {
  
  n_chunks <- ceiling(n_rep / chunk_size)
  all_results <- vector("list", n_chunks)
  
  for (ch in seq_len(n_chunks)) {
    chunk_file <- file.path(dir_out, sprintf("scenario_%s_chunk_%03d.rds", scenario_id, ch))
    
    # Se il chunk esiste già, caricalo
    if (file.exists(chunk_file)) {
      all_results[[ch]] <- readRDS(chunk_file)
      next
    }
    
    # Calcola range di repliche per questo chunk
    start_rep <- (ch - 1) * chunk_size + 1
    end_rep <- min(ch * chunk_size, n_rep)
    n_this_chunk <- end_rep - start_rep + 1
    
    # Esegui repliche
    chunk_results <- map_dfr(seq_len(n_this_chunk), function(r) {
      tryCatch(rep_fun(...), error = function(e) NULL)
    })
    
    # Salva checkpoint
    if (nrow(chunk_results) > 0) {
      saveRDS(chunk_results, chunk_file)
    }
    
    all_results[[ch]] <- chunk_results
  }
  
  bind_rows(all_results)
}

#' Carica tutti i risultati da una directory di chunk
#' 
#' @param dir_out Directory contenente i file .rds
#' @param pattern Pattern regex per filtrare i file (default: tutti .rds)
#' @return tibble con tutti i risultati
load_all_results <- function(dir_out, pattern = "\\.rds$") {
  files <- list.files(dir_out, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(tibble())
  
  results <- lapply(files, function(f) {
    tryCatch(readRDS(f), error = function(e) NULL)
  })
  
  bind_rows(results)
}

#' Combina chunk di uno scenario in un file riassunto
#' 
#' @param scenario_id ID dello scenario
#' @param dir_out Directory di output
#' @param summary_fun Funzione di aggregazione (prende tibble, restituisce tibble)
combine_chunks <- function(scenario_id, dir_out, summary_fun = NULL) {
  pattern <- sprintf("scenario_%s_chunk_", scenario_id)
  files <- list.files(dir_out, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) return(NULL)
  
  all_data <- bind_rows(lapply(files, readRDS))
  
  if (!is.null(summary_fun)) {
    all_data <- summary_fun(all_data)
  }
  
  # Salva riassunto e rimuovi chunk
  summary_file <- file.path(dir_out, sprintf("scenario_%s_summary.rds", scenario_id))
  saveRDS(all_data, summary_file)
  
  # Opzionale: rimuovi chunk dopo combinazione
  # file.remove(files)
  
  return(all_data)
}

#' Verifica se uno scenario è già completo
#' 
#' @param scenario_id ID dello scenario
#' @param n_rep Numero totale di repliche richieste
#' @param chunk_size Dimensione chunk
#' @param dir_out Directory di output
#' @return TRUE se completo, FALSE altrimenti
is_scenario_complete <- function(scenario_id, n_rep, chunk_size, dir_out) {
  # Prima controlla se esiste il file summary
  summary_file <- file.path(dir_out, sprintf("scenario_%s_summary.rds", scenario_id))
  if (file.exists(summary_file)) return(TRUE)
  

  # Altrimenti conta i chunk
  n_chunks_expected <- ceiling(n_rep / chunk_size)
  pattern <- sprintf("scenario_%s_chunk_", scenario_id)
  files <- list.files(dir_out, pattern = pattern)
  
  return(length(files) >= n_chunks_expected)
}


# =============================================================================
# FUNZIONI DI SIMULAZIONE (da Paper Functions – Nov 2025.R)
# =============================================================================

# 1A) Simulate survival cohort - INDIVIDUAL design
# Entrambi i trattamenti presenti in ogni ospedale (50/50 split)
simulate_survival_cohort_individual <- function(num_hosp, 
                                                 sample_size,
                                                 balancing_mode = 1,
                                                 icc, 
                                                 pop_treat_effect, 
                                                 lambda, 
                                                 gammas = 1,
                                                 cens = 2) {
  
  if (balancing_mode == 1) {
    num_pat_group <- rep(floor(sample_size / num_hosp), num_hosp)
    resto <- sample_size - sum(num_pat_group)
    if (resto > 0) num_pat_group[1:resto] <- num_pat_group[1:resto] + 1
  } 
  
  if (balancing_mode == 2) {
    mean_size <- sample_size / num_hosp
    pattern_name <- paste0("fixed_var_factor_", num_hosp)
    
    if (!exists(pattern_name, envir = .GlobalEnv)) {
      set.seed(num_hosp)
      assign(pattern_name, runif(num_hosp, 0.5, 1.5), envir = .GlobalEnv)
    }
    
    var_factor <- get(pattern_name, envir = .GlobalEnv)
    num_pat_group <- mean_size * var_factor
    num_pat_group <- num_pat_group / sum(num_pat_group) * sample_size
    num_pat_group <- round(num_pat_group)
    num_pat_group[num_pat_group <= 0] <- 1
    
    diff <- sample_size - sum(num_pat_group)
    if (diff != 0) {
      adjust_index <- sample(1:num_hosp, 1)
      num_pat_group[adjust_index] <- num_pat_group[adjust_index] + diff
    }
  }
  
  tot_patients <- sum(num_pat_group)
  var_resid <- (pi^2) / 6
  sigma2_hosp <- abs((icc / (1 - icc)) * var_resid)
  sigma_hosp <- sqrt(sigma2_hosp)
  
  hospital <- rep(1:num_hosp, times = num_pat_group)
  cov <- data.frame(
    id = 1:tot_patients,
    hospital = hospital,
    treat = rbinom(tot_patients, 1, 0.5)
  )
  
  cluster_effect <- rnorm(num_hosp, mean = 0, sd = sigma_hosp)
  cov$cluster_eff <- cluster_effect[cov$hospital]
  
  beta <- c(treat = pop_treat_effect, cluster_eff = 1)
  
  dati_sopravvivenza <- simsurv(
    dist = "weibull",
    lambdas = lambda,
    gammas = gammas,
    x = cov[, c("treat", "cluster_eff")],
    betas = beta,
    maxt = cens * (log(2) / lambda)
  )
  
  cohort <- merge(cov, dati_sopravvivenza, by = "id")
  cohort$sigma_hosp <- as.numeric(sigma_hosp[1])
  cohort$icc <- as.numeric(icc[1])
  cohort$cens_prop <- mean(cohort$status == 0)
  
  cohort[] <- lapply(cohort, function(x) if (is.numeric(x)) round(x, 3) else x)
  
  return(cohort)
}


# 1B) Simulate survival cohort - HOSPITAL design
# Ogni ospedale riceve un solo trattamento (cluster-randomized)
simulate_survival_cohort_hospital <- function(num_hosp, 
                                               sample_size,
                                               balancing_mode = 1,
                                               icc, 
                                               pop_treat_effect, 
                                               lambda, 
                                               gammas = 1,
                                               cens = 2) {
  
  if (balancing_mode == 1) {
    num_pat_group <- rep(floor(sample_size / num_hosp), num_hosp)
    resto <- sample_size - sum(num_pat_group)
    if (resto > 0) num_pat_group[1:resto] <- num_pat_group[1:resto] + 1
  } 
  
  if (balancing_mode == 2) {
    mean_size <- sample_size / num_hosp
    pattern_name <- paste0("fixed_var_factor_", num_hosp)
    
    if (!exists(pattern_name, envir = .GlobalEnv)) {
      set.seed(num_hosp)
      assign(pattern_name, runif(num_hosp, 0.5, 1.5), envir = .GlobalEnv)
    }
    
    var_factor <- get(pattern_name, envir = .GlobalEnv)
    num_pat_group <- mean_size * var_factor
    num_pat_group <- num_pat_group / sum(num_pat_group) * sample_size
    num_pat_group <- round(num_pat_group)
    num_pat_group[num_pat_group <= 0] <- 1
    
    diff <- sample_size - sum(num_pat_group)
    if (diff != 0) {
      adjust_index <- sample(1:num_hosp, 1)
      num_pat_group[adjust_index] <- num_pat_group[adjust_index] + diff
    }
  }
  
  tot_patients <- sum(num_pat_group)
  var_resid <- (pi^2) / 6
  sigma2_hosp <- abs((icc / (1 - icc)) * var_resid)
  sigma_hosp <- sqrt(sigma2_hosp)
  
  hospital <- rep(1:num_hosp, times = num_pat_group)
  
  # DIFFERENZA: trattamento assegnato a livello ospedale
  hospital_treat <- rbinom(num_hosp, 1, 0.5)
  treat <- hospital_treat[hospital]
  
  cov <- data.frame(
    id = 1:tot_patients,
    hospital = hospital,
    treat = treat
  )
  
  cluster_effect <- rnorm(num_hosp, mean = 0, sd = sigma_hosp)
  cov$cluster_eff <- cluster_effect[cov$hospital]
  
  beta <- c(treat = pop_treat_effect, cluster_eff = 1)
  
  dati_sopravvivenza <- simsurv(
    dist = "weibull",
    lambdas = lambda,
    gammas = gammas,
    x = cov[, c("treat", "cluster_eff")],
    betas = beta,
    maxt = cens * (log(2) / lambda)
  )
  
  cohort <- merge(cov, dati_sopravvivenza, by = "id")
  cohort$sigma_hosp <- as.numeric(sigma_hosp[1])
  cohort$icc <- as.numeric(icc[1])
  cohort$cens_prop <- mean(cohort$status == 0)
  
  cohort[] <- lapply(cohort, function(x) if (is.numeric(x)) round(x, 3) else x)
  
  return(cohort)
}


# =============================================================================
# 2) ICC ESTIMATION (tutti i 10 metodi)
# =============================================================================

icc_estimation <- function(data) {
  
  # Method 1 – Weibull / Log-normal
  data$logT <- log(data$eventtime)
  fit <- lmer(logT ~ 1 + (1 | hospital), data = data)
  var_cluster <- as.numeric(VarCorr(fit)$hospital[1])
  var_resid_theoretical <- pi^2 / 6
  icc_weibull <- var_cluster / (var_cluster + var_resid_theoretical)
  
  # Method 2 – CoxME Gaussian frailty
  mod_coxme <- tryCatch({
    coxme(Surv(eventtime, status) ~ 1 + (1|hospital), data = data)
  }, error = function(e) NULL)
  
  icc_coxme <- if (!is.null(mod_coxme)) {
    var_coxme <- VarCorr(mod_coxme)$hospital[[1]]
    var_coxme / (var_coxme + pi^2 / 6)
  } else { NA }
  
  # Method 3 – Cox with gamma frailty (non-parametric)
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
                    intTau.sigma2 = sigma2)$value - 1
    }
    
    icc_gamma_np <- fr.lognormal(k = 1, s = 1, sigma2 = var_gamma)
  } else { icc_gamma_np <- NA }
  
  # Method 4 – GLMM cloglog on discretized time
  q <- quantile(data$eventtime, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  q <- unique(q)
  tmin <- min(data$eventtime, na.rm = TRUE)
  q <- q[q > tmin]
  
  discretized_db <- survSplit(
    Surv(eventtime, status) ~ hospital,
    data = data,
    cut = q,
    start = "tstart",
    end = "tstop"
  )
  
  discretized_db$interval <- case_when(
    discretized_db$tstart == 0 ~ 1,
    discretized_db$tstart == q[1] ~ 2,
    length(q) > 1 & discretized_db$tstart == q[2] ~ 3,
    length(q) > 2 & discretized_db$tstart == q[3] ~ 4
  )
  
  mod_glmm <- tryCatch({
    glmer(
      status ~ as.factor(interval) + (1|hospital),
      data = discretized_db,
      family = binomial(link = "cloglog"),
      nAGQ = 5
    )
  }, error = function(e) NULL)
  
  icc_glmm <- if (!is.null(mod_glmm) && !isSingular(mod_glmm)) {
    var_glmm <- mod_glmm@theta^2
    var_glmm / (var_glmm + pi^2 / 6)
  } else { NA }
  
  # Method 5 – CoxME on discretized time
  mod_cox_exploded <- tryCatch({
    coxme(Surv(tstart, tstop, status) ~ 1 + (1|hospital), data = discretized_db)
  }, error = function(e) NULL)
  
  icc_cox_exploded <- if (!is.null(mod_cox_exploded)) {
    var_exploded <- VarCorr(mod_cox_exploded)$hospital[[1]]
    var_exploded / (var_exploded + pi^2 / 6)
  } else { NA }
  
  # Method 6 – Censoring indicators
  grand_mean <- mean(data$status, na.rm = TRUE)
  cluster_n <- tapply(data$status, data$hospital, length)
  cluster_mean <- tapply(data$status, data$hospital, mean)
  
  r <- length(cluster_n)
  MSB <- sum(cluster_n * (cluster_mean - grand_mean)^2) / (r - 1)
  
  cluster_mean_ind <- cluster_mean[data$hospital]
  MSW <- sum((data$status - cluster_mean_ind)^2) / (nrow(data) - r)
  
  m_bar <- mean(cluster_n)
  icc_censoring <- (MSB - MSW) / (MSB + (m_bar - 1) * MSW)
  
  # Method 7 – Observed event times only
  obs_data <- data[data$status == 1, ]
  obs_data$hospital <- as.factor(obs_data$hospital)
  obs_data$hospital <- droplevels(obs_data$hospital)
  
  if (nrow(obs_data) > 0 && length(unique(obs_data$hospital)) > 1) {
    cluster_means <- tapply(obs_data$eventtime, obs_data$hospital, mean)
    overall_mean <- mean(obs_data$eventtime)
    
    mi <- table(obs_data$hospital)
    r <- length(mi)
    N_star <- sum(mi)
    
    m_o <- (1 / (r - 1)) * (N_star - sum(mi^2) / N_star)
    
    MSB <- (1 / (r - 1)) * sum(mi * (cluster_means - overall_mean)^2)
    MSW <- (1 / (N_star - r)) * sum(tapply(obs_data$eventtime, obs_data$hospital,
                                            function(x) sum((x - mean(x))^2)))
    icc_event <- (MSB - MSW) / (MSB + (m_o - 1) * MSW)
  } else {
    icc_event <- NA
  }
  
  # Method 8 – Martingale residuals
  cox_fit <- coxph(Surv(eventtime, status) ~ 1, data = data)
  data$mart_resid <- residuals(cox_fit, type = "martingale")
  
  cluster_means_m <- data %>%
    group_by(hospital) %>%
    summarise(mean_resid = mean(mart_resid, na.rm = TRUE),
              n = n(), .groups = "drop")
  
  grand_mean_m <- mean(data$mart_resid, na.rm = TRUE)
  MSB_m <- sum(cluster_means_m$n * (cluster_means_m$mean_resid - grand_mean_m)^2) /
    (nrow(cluster_means_m) - 1)
  
  data <- left_join(data, cluster_means_m, by = "hospital")
  MSW_m <- sum((data$mart_resid - data$mean_resid)^2) /
    (nrow(data) - nrow(cluster_means_m))
  
  icc_martingale <- (MSB_m - MSW_m) / (MSB_m + (mean(cluster_means_m$n) - 1) * MSW_m)
  
  # Method 9 – Proportion of event (Xie & Waksman 2003)
  p_hat <- mean(data$status, na.rm = TRUE)
  data$dstatus <- data$status - p_hat
  
  agg <- by(data$dstatus, data$hospital, function(x) {
    s <- sum(x)
    ssq <- sum(x^2)
    list(K = length(x), pair_sum = (s^2 - ssq))
  })
  
  agg_list <- lapply(agg, unclass)
  Ks <- vapply(agg_list, `[[`, numeric(1), "K")
  pair_sums <- vapply(agg_list, `[[`, numeric(1), "pair_sum")
  
  numerator <- sum(pair_sums)
  denom_pairs <- sum(Ks * (Ks - 1))
  denominator <- p_hat * (1 - p_hat) * denom_pairs
  
  icc_proportion <- if (denominator != 0) numerator / denominator else NA
  
  # Method 10 – Weibull combined frailty (Oliveira 2016)
  icc_weibull_combined_fit <- tryCatch({
    coxme(Surv(eventtime, status) ~ 1 + (1|hospital), data = data)
  }, error = function(e) NULL)
  
  icc_weibull_combined <- if (!is.null(icc_weibull_combined_fit)) {
    D <- VarCorr(icc_weibull_combined_fit)$hospital[[1]]
    D / (D + 1)
  } else { NA }
  
  # Valore teorico (se presente)
  icc_theoretical <- if ("icc" %in% names(data)) {
    unique(data$icc)[1]
  } else { NA }
  
  # Output
  data.frame(
    Method = c(
      "Analytical Log-Weibull",
      "CoxME Gaussian frailty",
      "Cox gamma frailty (NP)",
      "GLMM logit discretized",
      "CoxME discretized time",
      "Censoring indicators",
      "Observed event times",
      "Martingale",
      "Weibull combined model",
      "Proportion of event",
      "Start ICC - 0"
    ),
    ICC = round(c(
      icc_weibull,
      icc_coxme,
      icc_gamma_np,
      icc_glmm,
      icc_cox_exploded,
      icc_censoring,
      icc_event,
      icc_martingale,
      icc_weibull_combined,
      icc_proportion,
      icc_theoretical
    ), 3)
  )
}


# =============================================================================
# 3) POWER FUNCTION (ottimizzata)
# =============================================================================

surv_power_function <- function(simula_coorte_fun, simula_args = list(), nsim = 100) {
  
  p_values_icc <- numeric(nsim)
  significant <- logical(nsim)
  prop_cens <- numeric(nsim)
  
  for (i in 1:nsim) {
    cohort <- do.call(simula_coorte_fun, simula_args)
    
    if (i == 1) {
      sigma_hosp <- unique(cohort$sigma_hosp)
      icc <- simula_args$icc
      num_hosp <- length(unique(cohort$hospital))
      num_pat_group_mean <- nrow(cohort) / num_hosp
      sample_size <- nrow(cohort)
    }
    
    prop_cens[i] <- mean(cohort$status == 0, na.rm = TRUE)
    
    mod_icc <- coxph(Surv(eventtime, status) ~ treat + cluster(hospital), data = cohort)
    p_val_icc <- summary(mod_icc)$coefficients["treat", "Pr(>|z|)"]
    p_values_icc[i] <- p_val_icc
    significant[i] <- p_val_icc < 0.05
  }
  
  power <- mean(significant)
  prop_cens_mean <- mean(prop_cens, na.rm = TRUE)
  
  tibble(
    nsim = nsim,
    icc = icc,
    num_hosp = num_hosp,
    num_pat_group_mean = num_pat_group_mean,
    sample_size = sample_size,
    power = power,
    prop_cens = prop_cens_mean
  )
}


# =============================================================================
# 4) SAMPLE SIZE FUNCTION (ottimizzata con binary search)
# =============================================================================

#' Ricerca sample size con binary search (più veloce del +1 iterativo)
sample_size_icc_binary <- function(icc, num_hosp, pop_treat_effect, 
                                   nsim = 500,
                                   min_n = 2, max_n = 300,
                                   target_power = 0.8,
                                   lambda = 0.1, gammas = 1, 
                                   cens = 2, balancing_mode = 1,
                                   simula_fun = simulate_survival_cohort_individual) {
  
  # Funzione per stimare power dato n per gruppo
  estimate_power <- function(npt) {
    significant <- logical(nsim)
    
    for (i in 1:nsim) {
      cohort <- tryCatch({
        simula_fun(
          num_hosp = num_hosp,
          sample_size = npt * num_hosp,
          icc = icc,
          pop_treat_effect = pop_treat_effect,
          lambda = lambda,
          gammas = gammas,
          cens = cens,
          balancing_mode = balancing_mode
        )
      }, error = function(e) NULL)
      
      if (is.null(cohort)) next
      
      mod <- tryCatch({
        coxph(Surv(eventtime, status) ~ treat + cluster(hospital), data = cohort)
      }, error = function(e) NULL)
      
      if (!is.null(mod)) {
        p_val <- summary(mod)$coefficients["treat", "Pr(>|z|)"]
        significant[i] <- p_val < 0.05
      }
    }
    
    mean(significant, na.rm = TRUE)
  }
  
  # Binary search
  low <- min_n
  high <- max_n
  result_n <- NA
  result_power <- NA
  
  # Prima verifica se max_n è sufficiente
  power_max <- estimate_power(high)
  if (power_max < target_power) {
    message(sprintf("Attenzione: power = %.2f con n = %d < target %.2f", 
                    power_max, high, target_power))
    result_n <- high
    result_power <- power_max
  } else {
    # Binary search
    while (high - low > 1) {
      mid <- floor((low + high) / 2)
      power_mid <- estimate_power(mid)
      
      message(sprintf("  n = %d, power = %.3f", mid, power_mid))
      
      if (power_mid >= target_power) {
        high <- mid
        result_n <- mid
        result_power <- power_mid
      } else {
        low <- mid
      }
    }
    
    # Affina con ultimo check
    if (is.na(result_n)) {
      result_n <- high
      result_power <- estimate_power(high)
    }
  }
  
  tibble(
    nsim = nsim,
    icc = icc,
    num_hosp = num_hosp,
    lambda = lambda,
    cens = cens,
    pop_treat_effect = pop_treat_effect,
    balancing_mode = balancing_mode,
    num_pat_group = result_n,
    sample_size = result_n * num_hosp,
    power = result_power
  )
}
