# ==============================================================================
# 00 - CONSTANTES
# ==============================================================================

# Lista de loci Y-STR
loci_y <- c(
  "dys576","dys389i","dys448","dys389ii","dys19","dys391","dys481","dys549",
  "dys533","dys438","dys437","dys570","dys635","dys390","dys439","dys392",
  "dys643","dys393","dys458","dys385a","dys385b","dys456","gatah4"
)

marcadores_minimos <- c(
  "dys19",
  "dys389i",
  "dys389ii",
  "dys390",
  "dys391",
  "dys392",
  "dys393",
  "dys385a",
  "dys385b"
)

MTDNA_NATIVE <- c("A2", "B2", "C1", "D1")

Y_NATIVE <- c("Q")

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
  "4+ Repeticiones" = "#2ea8df"
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

paleta_completa <- c(
  "Nativo"     = "#e41a1c",
  "Europeo"    = "#377eb8",
  "Africano"   = "#f0ad4e",
  "Haplogrupo" = "white"
)

paleta_texto <- c(
  "Nativo"     = "white",
  "Europeo"    = "white",
  "Africano"   = "black",
  "Haplogrupo" = "black"
)