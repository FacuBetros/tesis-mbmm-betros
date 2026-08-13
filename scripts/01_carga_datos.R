# ==============================================================================
# 01 - CARGA DE DATOS
# ==============================================================================
excel_file <- here::here("datos", "Tesis.xlsx")

if (!file.exists(excel_file)) {
  stop("CRÍTICO: No se encontró 'datos/Tesis.xlsx'. Verifica la carpeta.")
}

message("Cargando bases de datos...")

sheets <- c(
  muestras = "Muestras",
  ancestria_sistemas = "Ancestria_Sistemas",
  ancestria_integrada = "Ancestria_Integrada",
  yhrd = "YHRD",
  microvariantes = "Microvariantes",
  comparativas = "Poblaciones_Comparativas"
)

datos <- purrr::map(
  sheets,
  ~ read_excel(excel_file, sheet = .x) %>%
    clean_names()
)

muestras <- datos$muestras
ancestria_sistemas  <- datos$ancestria_sistemas
ancestria_integrada <- datos$ancestria_integrada
yhrd                <- datos$yhrd
microvariantes      <- datos$microvariantes
comparativas        <- datos$comparativas

message(">> Datos cargados correctamente.")


# ------------------------------------------------------------------------------
# 1. PRE-PROCESAMIENTO DE MUESTRAS
# ------------------------------------------------------------------------------

#Validacion de columnas
required_cols <- c(
  "id",
  "hg_y",
  "cat_y",
  "mt_dna"
)

missing_cols <- setdiff(required_cols, names(muestras))

if(length(missing_cols) > 0) {
  stop(
    "Faltan columnas: ",
    paste(missing_cols, collapse = ", ")
  )
}

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
    origen_y    = y_origin(cat_y)
    )

# Unir haplotipo completo
missing_loci <- setdiff(loci_y, names(muestras))

if(length(missing_loci) > 0){
  stop(
    "Loci faltantes: ",
    paste(missing_loci, collapse = ", ")
  )
}

muestras <- muestras %>% unite(haplotipo, all_of(loci_y), sep = "-", remove = FALSE)
