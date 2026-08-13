# ==============================================================================
# 00 - FUNCIONES
# Autor: Facundo Betros
# ==============================================================================

# --- Limpieza de texto ---------------------------------------------------------
clean_text <- function(x) {
  x %>%
    as.character() %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    stringr::str_squish() %>%
    toupper()
}

# --- Limpieza de números -------------------------------------------------------
clean_numeric <- function(x) {
  suppressWarnings(
    x %>%
      as.character() %>%
      gsub("%", "", .) %>%
      gsub(",", ".", .) %>%
      trimws() %>%
      as.numeric()
  )
}

# --- Clasificación mtDNA -------------------------------------------------------
mtdna_native <- function(x) {
  dplyr::if_else(
    x %in% MTDNA_NATIVE,
    "SI",
    "NO"
  )
}

# --- Clasificación Y -----------------------------------------------------------
y_origin <- function(x) {
  dplyr::if_else(
    x %in% Y_NATIVE,
    "NATIVO",
    "NO_NATIVO"
  )
}

# --- Redondeo de Valores -----------------------------------------------------------
label_pct <- function(x) {
  ifelse(
    is.na(x),
    "",
    sprintf("%.1f%%", x)
  )
}