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


# ------------------------------------------------------------------------------
# 1. PRE-PROCESAMIENTO DE MUESTRAS
# ------------------------------------------------------------------------------

# Limpieza de texto y creación de variables auxiliares
muestras <- muestras %>%
  mutate(
    id          = clean_text(id),
    hg_y        = clean_text(hg_y),
    cat_y       = clean_text(cat_y),
    mt_dna      = clean_text(mt_dna),
    shg_nevgen  = clean_text(shg_nevgen),
    shg_hprg    = clean_text(shg_hprg),
    mt_dna_nat  = mtdna_native(mt_dna),
    origen_y    = y_origin(cat_y),
    grupo_y     = case_when(
      cat_y == "Q" ~ "NATIVO",
      cat_y %in% c("R", "I", "G") ~ "EUROPEO",
      cat_y %in% c("J", "E1B1B", "T") ~ "MEDITERRANEO",
      cat_y %in% c("A", "B", "E1A", "E1B1A") ~ "AFRICANO",
      cat_y %in% c("L", "O", "H") ~ "ASIATICO",
      TRUE ~ "OTROS"
    )
  )

# Lista de loci Y-STR
loci_y <- c(
  "dys576","dys389i","dys448","dys389ii","dys19","dys391","dys481","dys549",
  "dys533","dys438","dys437","dys570","dys635","dys390","dys439","dys392",
  "dys643","dys393","dys458","dys385a","dys385b","dys456","gatah4"
)

# Unir haplotipo completo
muestras <- muestras %>% unite(haplotipo, all_of(loci_y), sep = "-", remove = FALSE)
