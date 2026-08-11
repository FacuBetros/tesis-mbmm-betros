# ==============================================================================
# TESIS: GENÉTICA FORENSE Y POBLACIONAL
# Autor: Facundo Betros
# ==============================================================================

# 0. CONFIGURACIÓN --------------------------------------------------------------
set.seed(123)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  tidyverse, readxl, janitor, stringi, maps, mapdata, countrycode,
  ggtern, grid, gridExtra, irr, dplyr, stringr, tibble, FSA, ape,
  ggtree, dendextend, factoextra, ggplot2, here
)

# Ruta segura para GitHub
excel_file <- here::here("datos", "Tesis.xlsx")

if (!file.exists(excel_file)) {
  stop("CRÍTICO: No se encontró 'datos/Tesis.xlsx'. Verifica la carpeta.")
}

message("Cargando bases de datos...")
