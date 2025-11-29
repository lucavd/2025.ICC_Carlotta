# =============================================================================
# 00_setup.R - Setup e configurazione per simulazioni ICC
# =============================================================================
# Esegui questo script UNA VOLTA per installare i pacchetti necessari
# e verificare che tutto funzioni.

# --- Pacchetti richiesti ---
required_packages <- c(

"tidyverse",
"furrr",
"future",
"progressr",
"simsurv",
"survival",
"lme4",
"coxme",
"lhs"
)

# Installa i pacchetti mancanti
install_if_missing <- function(pkg) {
if (!requireNamespace(pkg, quietly = TRUE)) {
install.packages(pkg, repos = "https://cloud.r-project.org")
}
}

invisible(lapply(required_packages, install_if_missing))

# --- Verifica installazione ---
cat("Verifica pacchetti:\n")
for (pkg in required_packages) {
  loaded <- require(pkg, character.only = TRUE, quietly = TRUE)
  cat(sprintf("  %s: %s\n", pkg, ifelse(loaded, "OK", "ERRORE")))
}

# --- Test rapido furrr ---
cat("\nTest parallelizzazione:\n")
plan(multisession, workers = 2)
test_result <- future_map(1:4, ~.x^2, .options = furrr_options(seed = TRUE))
cat(sprintf("  future_map test: %s\n", 
            ifelse(identical(unlist(test_result), c(1L, 4L, 9L, 16L)), "OK", "ERRORE")))
plan(sequential)

cat("\nSetup completato.\n")
