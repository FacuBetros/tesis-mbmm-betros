# ==============================================================================
# MAIN - EJECUCIÓN COMPLETA DE LA TESIS
# Autor: Facundo Betros
# ==============================================================================

message(">> Iniciando pipeline completo de la tesis...")

source(here::here("scripts", "00_funciones.R"))
source(here::here("scripts", "01_carga_datos.R"))
source(here::here("scripts", "02_objetivo1.R"))
source(here::here("scripts", "03_objetivo2.R"))
source(here::here("scripts", "04_objetivo3.R"))
source(here::here("scripts", "05_objetivo4.R"))
source(here::here("scripts", "06_objetivo5.R"))

message(">> Pipeline completado correctamente.")
