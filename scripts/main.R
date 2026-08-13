# ==============================================================================
# MAIN - EJECUCIÓN COMPLETA DE LA TESIS
# Autor: Facundo Betros
# ==============================================================================
set.seed(123) # Para reproducibilidad

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
  FSA,
  ape,
  ggtree,
  dendextend,
  factoextra,
  ggforce,
  hierfstat,
  pheatmap,
  patchwork
)


message(">> Iniciando pipeline completo de la tesis...")

source(here::here("scripts", "02_constantes.R"))
source(here::here("scripts", "03_funciones.R"))
source(here::here("scripts", "01_carga_datos.R"))
source(here::here("scripts", "04_objetivo1.R"))
source(here::here("scripts", "05_objetivo2.R"))
source(here::here("scripts", "06_objetivo3.R"))
source(here::here("scripts", "07_objetivo4.R"))
source(here::here("scripts", "08_objetivo5.R"))

message(">> Pipeline completado correctamente.")
