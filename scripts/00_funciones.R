# ==============================================================================
# 00 - FUNCIONES AUXILIARES
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
  if_else(x %in% c("A2", "B2", "C1", "D1"), "SI", "NO")
}

# --- Clasificación Y -----------------------------------------------------------
y_origin <- function(x) {
  if_else(x == "Q", "NATIVO", "NO_NATIVO")
}

# --- Paletas de colores --------------------------------------------------------
paleta_ancestria <- c(
  "Africano" = "#f0ad4e",
  "Europeo"  = "#377eb8",
  "Nativo"   = "#e41a1c"
)

paleta_continentes <- c(
  "Africa"   = "#f0ad4e",
  "America"  = "#e41a1c",
  "Asia"     = "#4daf4a",
  "Europe"   = "#377eb8",
  "Oceania"  = "#984ea3"
)
