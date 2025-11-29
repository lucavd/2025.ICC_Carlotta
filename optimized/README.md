# ICC Simulation Study - Pipeline Ottimizzata

## Struttura File

```
optimized/
├── 00_setup.R              # Installa pacchetti (esegui UNA volta)
├── 01_config.R             # Configurazione parametri e grid
├── 02_functions.R          # Funzioni simulazione + helper chunking
├── 03_paper1_icc.R         # Paper 1: confronto stimatori ICC
├── 04_paper2_power.R       # Paper 2: impatto ICC su potenza
├── 05_paper2_sample_size.R # Paper 2: sample size per power 80%
├── run_all.R               # Script master
└── README.md               # Questo file
```

## Quick Start

```bash
# 1. Setup (solo la prima volta)
cd optimized
Rscript 00_setup.R

# 2. Esegui tutto
Rscript run_all.R

# Oppure solo parti specifiche:
Rscript run_all.R paper1
Rscript run_all.R paper2
```

## Configurazione Server HPC

Modifica `01_config.R`:

```r
# Imposta il numero di core del server
N_CORES <- 110  # <-- modifica qui
```

## Riduzione Carico Computazionale

Rispetto agli script originali:

| Aspetto | Originale | Ottimizzato | Risparmio |
|---------|-----------|-------------|-----------|
| **Paper 1 n_rep** | 2000 | 1000 | 50% |
| **Paper 1 scenari** | 768×2 | 432×2 | 44% |
| **Paper 2 nsim** | 2000 | 1000 | 50% |
| **Paper 2 scenari** | 2688×2 | 960×2 | 64% |
| **Sample size** | +1 iterativo | Binary search | ~70% |
| **Totale stimato** | ~30+ ore su 110 core | ~5-8 ore | **75%** |

## Checkpoint e Ripresa

- Ogni scenario salva risultati parziali in `results_optimized/`
- Se lo script si interrompe, rilancialo: riprende da dove si è fermato
- I chunk (200 repliche ciascuno) sono salvati separatamente

## Modifica Parametri

In `01_config.R` puoi modificare:

```r
# Numero repliche Monte Carlo
PAPER1_N_REP <- 1000        # default: 1000, originale: 2000
PAPER2_POWER_NSIM <- 1000   # default: 1000, originale: 2000
PAPER2_SS_NSIM <- 500       # default: 500

# Grid parametri (esempio Paper 1)
PAPER1_PARAMS <- list(
  betas = c(-0.2, -0.5),
  lambdas = c(0.115, 0.012),
  sample_sizes = c(100, 500, 2000),  # puoi aggiungere/togliere
  num_hosps = c(5, 30, 60),
  iccs = c(0.01, 0.1, 0.3),
  cens_values = c(0.5, 5),
  balancing_modes = c(1, 2)
)
```

## Output

Tutti i risultati sono salvati in `results_optimized/`:

```
results_optimized/
├── paper1_icc_individual/
│   ├── scenario_001_chunk_001.rds  # chunk parziali
│   ├── scenario_001_summary.rds    # riassunto scenario
│   └── ICC_results_individual_FINAL.rds
├── paper1_icc_hospital/
│   └── ICC_results_hospital_FINAL.rds
├── paper2_power_individual/
│   └── power_results_individual_FINAL.rds
├── paper2_power_hospital/
│   └── power_results_hospital_FINAL.rds
├── paper2_sample_size/
│   └── sample_size_individual_FINAL.rds
├── paper1_ICC_ALL_RESULTS.rds       # combinato
├── paper2_POWER_ALL_RESULTS.rds     # combinato
└── paper2_SAMPLE_SIZE_ALL_RESULTS.rds
```

## Caricare Risultati

```r
library(tidyverse)

# Paper 1
icc_results <- readRDS("results_optimized/paper1_ICC_ALL_RESULTS.rds")

# Paper 2
power_results <- readRDS("results_optimized/paper2_POWER_ALL_RESULTS.rds")
ss_results <- readRDS("results_optimized/paper2_SAMPLE_SIZE_ALL_RESULTS.rds")
```

## Alternativa: Latin Hypercube Sampling

Se vuoi esplorare lo spazio parametri con ancora meno punti:

```r
source("01_config.R")

# Genera 200 punti LHS invece di grid completa
scenarios_lhs <- make_lhs_scenarios(n_points = 200)
```

## Troubleshooting

**Errore "cannot allocate vector"**: Riduci `N_CORES` o `PAPER1_N_REP`

**Scenario troppo lento**: I casi con `sample_size = 6000` e `num_hosp = 60` sono i più pesanti. Puoi escluderli dalla grid.

**Riprendere da interruzione**: Basta rilanciare lo stesso script; i file `.rds` esistenti vengono saltati.
