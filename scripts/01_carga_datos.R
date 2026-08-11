# ==============================================================================
# 01 - CARGA DE DATOS
# ==============================================================================

set.seed(123)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  tidyverse, readxl, janitor, stringi, maps, mapdata, countrycode,
  ggtern, grid, gridExtra, irr, dplyr, stringr, tibble, FSA, ape,
  ggtree, dendextend, factoextra, ggplot2, here
)

excel_file <- here::here("datos", "Tesis.xlsx")

if (!file.exists(excel_file)) {
  stop("CRÍTICO: No se encontró 'datos/Tesis.xlsx'. Verifica la carpeta.")
}

message("Cargando bases de datos...")

muestras              <- read_excel(excel_file, sheet = "Muestras") %>% clean_names()
ancestria_sistemas    <- read_excel(excel_file, sheet = "Ancestria_Sistemas") %>% clean_names()
ancestria_integrada   <- read_excel(excel_file, sheet = "Ancestria_Integrada") %>% clean_names()
yhrd                  <- read_excel(excel_file, sheet = "YHRD") %>% clean_names()
microvariantes        <- read_excel(excel_file, sheet = "Microvariantes") %>% clean_names()
comparativas          <- read_excel(excel_file, sheet = "Poblaciones_Comparativas") %>% clean_names()

message(">> Datos cargados correctamente.")
