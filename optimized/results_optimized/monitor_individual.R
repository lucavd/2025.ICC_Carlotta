##############################################################################
# MONITOR SCENARI INDIVIDUAL MANCANTI - notifica Telegram ogni ora
##############################################################################

library(telegram.bot)

# Bot config
bot <- telegram.bot::Bot(token = "1089304740:AAHA_buNHuJ2XMTHvCzrPgFaiZDApafHCWo")
chat_id <- "118334609"

# Cartella da monitorare
dir_output <- "/home/user/2025.ICC_carlotta/optimized/results_optimized/individual_mancanti_partial"
total_scenarios <- 241

# Milestone da notificare
milestones <- c(100, 150, 200, total_scenarios)
notified <- c()

# Funzione per contare file completati
count_completed <- function() {
  files <- list.files(dir_output, pattern = "\\.rds$")
  length(files)
}

# Funzione per inviare messaggio
send_msg <- function(msg) {
  bot$sendMessage(chat_id = chat_id, text = msg)
  cat(sprintf("[%s] Inviato: %s\n", Sys.time(), msg))
}

# Messaggio iniziale
send_msg(sprintf("🚀 Monitor avviato! Totale scenari: %d\nCheck ogni ora.", total_scenarios))

# Loop principale
repeat {
  n_done <- count_completed()
  cat(sprintf("[%s] Completati: %d/%d\n", Sys.time(), n_done, total_scenarios))
  
  # Check milestone
  for (m in milestones) {
    if (n_done >= m && !(m %in% notified)) {
      if (m == total_scenarios) {
        send_msg(sprintf("✅ COMPLETATO! %d/%d scenari finiti!", n_done, total_scenarios))
        quit(save = "no")
      } else {
        send_msg(sprintf("📊 Progresso: %d/%d scenari completati (milestone %d)", n_done, total_scenarios, m))
      }
      notified <- c(notified, m)
    }
  }
  
  # Attendi 1 ora (3600 secondi)
  Sys.sleep(3600)
}
