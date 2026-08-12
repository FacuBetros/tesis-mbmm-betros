# ==============================================================================
# 00 - FUNCIONES AUXILIARES
# Autor: Facundo Betros
# ==============================================================================

# --- Librerías ----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  tidyverse,
  readxl,
  janitor,
  stringi,
  maps,
  mapdata,
  countrycode,
  ggtern,
  grid,
  gridExtra,
  irr,
  dplyr,
  stringr,
  tibble,
  FSA,
  ape,
  ggtree,
  dendextend,
  factoextra,
  ggplot2,
  ggforce,
  tidyr
)

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

paleta_hg <- c(
  "R" = "#377eb8", "Q" = "#e41a1c", "E1B1B" = "#31a354", "E1B1A" = "#ff7f00",
  "E1A" = "#a1d99b", "J" = "#756bb1", "I" = "#fd8d3c", "G" = "#bf812d",
  "L" = "#00bfc4", "T" = "#f781bf", "O" = "#984ea3", "A" = "#8c510a",
  "B" = "#b35806", "ND" = "grey50", "NO" = "grey75", "H" = "#ff1493"
)

paleta_haplotipos <- c(
  "Únicos" = "#1f6285",
  "2 Repeticiones" = "#e8752f",
  "3 Repeticiones" = "#1f7a2e",
  "4 Repeticiones" = "#2ea8df"
)

paleta_grupos <- c(
  "Europeo" = "#377eb8",
  "Mediterraneo" = "#756bb1",
  "Nativo Americano" = "#e41a1c",
  "Africano" = "#31a354",
  "Asiatico" = "#ff7f00",
  "Otros" = "grey70"
)

paleta_mt <- c(
  "A2"="#e41a1c",
  "B2"="#f781bf",
  "C1"="#a65628",
  "D1"="#984ea3"
)

paleta_sistemas <- c(
  "STRs Autosómicos"   = "#377eb8",
  "InDels Autosómicos" = "#4daf4a",
  "AIMs"               = "#e41a1c",
  "Cromosoma X"        = "#f0ad4e"
)

