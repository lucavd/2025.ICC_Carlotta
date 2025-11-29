# =============================================================================
# 01_config.R - Configurazione centralizzata e riduzione statistica del carico
# =============================================================================
# Questo file definisce TUTTI i parametri di simulazione.
# Modifica QUI per cambiare scenari, n_rep, n_cores, ecc.

# --- Carica pacchetti ---
library(tidyverse)
library(furrr)
library(future)
library(progressr)
library(simsurv)
library(survival)
library(lme4)
library(coxme)
library(lhs)  # per Latin Hypercube Sampling (opzionale)

# --- Parallelizzazione ---
# Imposta il numero di core. Su server HPC, usa tutti quelli disponibili.
N_CORES <- parallel::detectCores() - 1
# N_CORES <- 110  # <-- decommentare e impostare manualmente per HPC

# Inizializza il piano future (da chiamare PRIMA di ogni script)
setup_parallel <- function(n_workers = N_CORES) {
  plan(multisession, workers = n_workers)
  handlers(global = TRUE)  # abilita progress bar
  message(sprintf("Piano parallelo attivato: %d workers", n_workers))
}

# --- Directory output ---
DIR_BASE <- "results_optimized"
DIR_PAPER1_IND <- file.path(DIR_BASE, "paper1_icc_individual")
DIR_PAPER1_HOSP <- file.path(DIR_BASE, "paper1_icc_hospital")
DIR_PAPER2_POWER_IND <- file.path(DIR_BASE, "paper2_power_individual")
DIR_PAPER2_POWER_HOSP <- file.path(DIR_BASE, "paper2_power_hospital")
DIR_PAPER2_SS <- file.path(DIR_BASE, "paper2_sample_size")

create_dirs <- function() {
  dirs <- c(DIR_PAPER1_IND, DIR_PAPER1_HOSP, 
            DIR_PAPER2_POWER_IND, DIR_PAPER2_POWER_HOSP, 
            DIR_PAPER2_SS)
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  invisible(NULL)
}

# =============================================================================
# RIDUZIONE STATISTICA DEL CARICO
# =============================================================================

# -----------------------------------------------------------------------------
# PAPER 1 - ICC Estimation
# -----------------------------------------------------------------------------
# ORIGINALE: n_rep = 2000, 768 scenari × 2 disegni = 1536 scenari
# RIDOTTO:   n_rep = 1000, grid più rada → ~384 scenari × 2 = 768 scenari
#            Risparmio: ~75%

PAPER1_N_REP <- 1000   # ORIGINALE: 2000. Ridotto a 1000 (sufficiente per IC 95%)
PAPER1_CHUNK_SIZE <- 200  # Salva checkpoint ogni 200 repliche

# Grid ORIGINALE (768 scenari per disegno):
# betas: 2, lambdas: 2, sample_sizes: 4, num_hosps: 3, iccs: 4, cens: 2, balancing: 2
# = 2*2*4*3*4*2*2 = 768

# Grid RIDOTTA (mantieni i casi più informativi):
PAPER1_PARAMS <- list(
  betas = c(-0.2, -0.5),              # HR 0.8, 0.6 (entrambi importanti)
  lambdas = c(0.115, 0.012),          # mediana 6m, 60m (entrambi)
  sample_sizes = c(100, 500, 2000),   # RIDOTTO: tolto 1000, 6000 → 3 livelli
  num_hosps = c(5, 30, 60),           # RIDOTTO: tolto 15 → 3 livelli
  iccs = c(0.01, 0.1, 0.3),           # RIDOTTO: tolto 0.06 → 3 livelli (basso, medio, alto)
  cens_values = c(0.5, 5),            # mantenuti (alta e bassa censura)
  balancing_modes = c(1, 2)           # mantenuti (bilanciato, sbilanciato)
)
# Nuovi scenari: 2*2*3*3*3*2*2 = 432 per disegno = 864 totali
# Con n_rep = 1000: 864,000 simulazioni (vs 3,072,000 originali) → -72%

# Genera grid Paper 1
make_paper1_scenarios <- function() {
  expand.grid(
    beta = PAPER1_PARAMS$betas,
    lambda = PAPER1_PARAMS$lambdas,
    sample_size = PAPER1_PARAMS$sample_sizes,
    num_hosp = PAPER1_PARAMS$num_hosps,
    icc = PAPER1_PARAMS$iccs,
    cens = PAPER1_PARAMS$cens_values,
    balancing_mode = PAPER1_PARAMS$balancing_modes,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# PAPER 2 - Power
# -----------------------------------------------------------------------------
# ORIGINALE: nsim = 2000, 2688 scenari × 2 disegni + DE = ~8000 runs
# RIDOTTO:   nsim = 1000, grid più rada, skip DE per ICC=0

PAPER2_POWER_NSIM <- 1000  # ORIGINALE: 2000
PAPER2_POWER_CHUNK_SIZE <- 200

# Grid ORIGINALE (2688 scenari per disegno):
# treat_effects: 2, lambdas: 2, icc_values: 7, num_hosps: 4, sample_sizes: 4, cens: 3, balancing: 2
# = 2*2*7*4*4*3*2 = 2688

# Grid RIDOTTA:
PAPER2_POWER_PARAMS <- list(
  treat_effects = c(-0.2, -0.5),
  lambdas = c(0.115, 0.012),
  icc_values = c(0, 0.01, 0.05, 0.1, 0.3),  # RIDOTTO: 5 livelli (tolti duplicati, 0.02, 0.03)
  num_hosps = c(5, 15, 60),                 # RIDOTTO: 3 livelli (tolto 30)
  sample_sizes = c(100, 500, 2000, 6000),   # mantenuti 4
  cens_values = c(0.5, 5),                  # RIDOTTO: 2 livelli (tolto 2)
  balancing_modes = c(1, 2)
)
# Nuovi scenari: 2*2*5*3*4*2*2 = 960 per disegno
# Con nsim = 1000: 960,000 sim per disegno (vs 5,376,000) → -82%

# NOTA: Per ICC = 0 non serve ricalcolare con DE (DE = 1)
# Questo elimina ~1/5 degli scenari DE

make_paper2_power_scenarios <- function() {
  expand.grid(
    pop_treat_effect = PAPER2_POWER_PARAMS$treat_effects,
    lambda = PAPER2_POWER_PARAMS$lambdas,
    icc = PAPER2_POWER_PARAMS$icc_values,
    num_hosp = PAPER2_POWER_PARAMS$num_hosps,
    sample_size = PAPER2_POWER_PARAMS$sample_sizes,
    cens = PAPER2_POWER_PARAMS$cens_values,
    balancing_mode = PAPER2_POWER_PARAMS$balancing_modes,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# PAPER 2 - Sample Size
# -----------------------------------------------------------------------------
# ORIGINALE: nsim = 2000 × ~40 step di ricerca = ~80,000 sim per scenario, 96 scenari
# OTTIMIZZATO: nsim = 500, binary search invece di +1, early stopping

PAPER2_SS_NSIM <- 500      # ORIGINALE: 2000. Ridotto (per ricerca iterativa basta meno)
PAPER2_SS_MAX_PAT <- 300   # ORIGINALE: 500. Ridotto
PAPER2_SS_TIMEOUT <- 180   # secondi timeout per scenario

# Grid ORIGINALE: 96 scenari
# Grid RIDOTTA (scenari più informativi):
PAPER2_SS_PARAMS <- list(
  pte_values = c(-0.2, -0.5),
  lambda_values = c(0.115, 0.012),
  n_hosp_values = c(5, 30, 60),       # ESTESO: aggiungo 30 per coprire meglio
  icc_values = c(0, 0.01, 0.05, 0.1), # 4 livelli
  cens_values = c(0.5, 5),
  balancing_modes = c(1, 2)
)
# Nuovi scenari: 2*2*3*4*2*2 = 192 scenari (più che originale, ma con meno sim/scenario)

make_paper2_ss_scenarios <- function() {
  expand.grid(
    pop_treat_effect = PAPER2_SS_PARAMS$pte_values,
    lambda = PAPER2_SS_PARAMS$lambda_values,
    num_hosp = PAPER2_SS_PARAMS$n_hosp_values,
    icc = PAPER2_SS_PARAMS$icc_values,
    cens = PAPER2_SS_PARAMS$cens_values,
    balancing_mode = PAPER2_SS_PARAMS$balancing_modes,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# OPZIONE: LATIN HYPERCUBE SAMPLING (alternativa a grid completa)
# =============================================================================
# Se vuoi esplorare lo spazio dei parametri con meno punti ma buona copertura

make_lhs_scenarios <- function(n_points = 200, seed = 42) {
  set.seed(seed)
  
  # Genera punti LHS normalizzati [0,1]
  lhs_matrix <- lhs::randomLHS(n_points, 7)
  
  # Mappa su valori reali
  scenarios <- data.frame(
    beta = ifelse(lhs_matrix[,1] < 0.5, -0.2, -0.5),
    lambda = ifelse(lhs_matrix[,2] < 0.5, 0.115, 0.012),
    sample_size = round(100 + lhs_matrix[,3] * (6000 - 100)),
    num_hosp = round(5 + lhs_matrix[,4] * (60 - 5)),
    icc = lhs_matrix[,5] * 0.3,  # continuo tra 0 e 0.3
    cens = ifelse(lhs_matrix[,6] < 0.5, 0.5, 5),
    balancing_mode = ifelse(lhs_matrix[,7] < 0.5, 1, 2)
  )
  
  # Arrotonda num_hosp a valori sensati
  scenarios$num_hosp <- pmax(5, scenarios$num_hosp)
  
  return(scenarios)
}

# =============================================================================
# SOMMARIO RIDUZIONE CARICO
# =============================================================================
print_reduction_summary <- function() {
  cat("\n========== SOMMARIO RIDUZIONE CARICO ==========\n\n")
  
  # Paper 1
  orig_p1 <- 2 * 2 * 4 * 3 * 4 * 2 * 2 * 2000 * 2  # scenari × n_rep × disegni

  new_p1 <- nrow(make_paper1_scenarios()) * PAPER1_N_REP * 2
  cat(sprintf("PAPER 1 (ICC estimation):\n"))
  cat(sprintf("  Originale: %s simulazioni\n", format(orig_p1, big.mark = ",")))
  cat(sprintf("  Ridotto:   %s simulazioni\n", format(new_p1, big.mark = ",")))
  cat(sprintf("  Risparmio: %.0f%%\n\n", (1 - new_p1/orig_p1) * 100))
  
  # Paper 2 Power
  orig_p2p <- 2 * 2 * 7 * 4 * 4 * 3 * 2 * 2000 * 2
  new_p2p <- nrow(make_paper2_power_scenarios()) * PAPER2_POWER_NSIM * 2
  cat(sprintf("PAPER 2 (Power):\n"))
  cat(sprintf("  Originale: %s simulazioni\n", format(orig_p2p, big.mark = ",")))
  cat(sprintf("  Ridotto:   %s simulazioni\n", format(new_p2p, big.mark = ",")))
  cat(sprintf("  Risparmio: %.0f%%\n\n", (1 - new_p2p/orig_p2p) * 100))
  
  # Paper 2 Sample Size (stima approssimativa)
  orig_p2s <- 96 * 2000 * 40  # scenari × nsim × step medi
  new_p2s <- nrow(make_paper2_ss_scenarios()) * PAPER2_SS_NSIM * 15  # meno step con binary search
  cat(sprintf("PAPER 2 (Sample Size):\n"))
  cat(sprintf("  Originale: ~%s simulazioni\n", format(orig_p2s, big.mark = ",")))
  cat(sprintf("  Ridotto:   ~%s simulazioni\n", format(new_p2s, big.mark = ",")))
  cat(sprintf("  Risparmio: ~%.0f%%\n\n", (1 - new_p2s/orig_p2s) * 100))
  
  cat("================================================\n")
}

# Esegui per vedere il sommario
# print_reduction_summary()
