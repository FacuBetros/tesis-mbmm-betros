# ==============================================================================
# 06 - OBJETIVO 5: MESTIZAJE Y TRAYECTORIAS GENÓMICAS
# Autor: Facundo Betros
# ==============================================================================
#
# Requiere haber ejecutado previamente:
#   00funciones.R
#   01cargadatos.R
#
# Objetos esperados:
#   ancestria_integrada
#   clean_numeric()
#   paleta_ancestria
#
# Salidas:
#   outputs/tablas/
#   outputs/figuras/
# ==============================================================================

message(">> Ejecutando Objetivo 5...")

# ------------------------------------------------------------------------------
# 0. CONFIGURACIÓN GENERAL
# ------------------------------------------------------------------------------

set.seed(123)

dir.create(here::here("outputs", "tablas"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("outputs", "figuras"), showWarnings = FALSE, recursive = TRUE)

# Umbral mínimo para mostrar porcentajes dentro de gráficos.
# Si querés gráficos más limpios, subir a 6, 7 u 8.
# Si querés mostrar más detalle, bajar a 3 o 4.
UMBRAL_ETIQUETA_PCT <- 5

# Límite de caracteres para etiquetas de haplogrupos dentro de barras.
MAX_CHARS_HAPLO <- 12

# ------------------------------------------------------------------------------
# Funciones auxiliares específicas del Objetivo 5
# ------------------------------------------------------------------------------

pctificar <- function(x) {
  x_num <- clean_numeric(x)
  
  dplyr::if_else(
    !is.na(x_num) & x_num <= 1,
    x_num * 100,
    x_num
  )
}

label_pct <- function(x) {
  x_num <- as.numeric(x)
  
  dplyr::if_else(
    !is.na(x_num),
    paste0(round(x_num, 1), "%"),
    ""
  )
}

label_haplo <- function(x, max_chars = MAX_CHARS_HAPLO) {
  stringr::str_trunc(
    dplyr::coalesce(as.character(x), "S/D"),
    width = max_chars
  )
}

clasificar_y_general <- function(x) {
  x_chr <- toupper(dplyr::coalesce(as.character(x), ""))
  
  dplyr::case_when(
    stringr::str_detect(x_chr, "^Q") ~ "NATIVO",
    TRUE ~ "NO_NATIVO"
  )
}

clasificar_y_plot <- function(x) {
  x_chr <- toupper(dplyr::coalesce(as.character(x), ""))
  
  dplyr::case_when(
    stringr::str_detect(x_chr, "^Q") ~ "Nativo",
    stringr::str_detect(x_chr, "^(A|B|E)") ~ "Africano",
    x_chr == "" ~ NA_character_,
    TRUE ~ "Europeo"
  )
}

clasificar_mt_general <- function(x) {
  x_chr <- toupper(dplyr::coalesce(as.character(x), ""))
  
  dplyr::case_when(
    stringr::str_detect(x_chr, "^[ABCD]") ~ "NATIVO",
    TRUE ~ "NO_NATIVO"
  )
}

# Paletas del Objetivo 5.
# Se mantiene compatibilidad con las paletas generales de 00funciones.
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

# ------------------------------------------------------------------------------
# 1. PREPARACIÓN Y LIMPIEZA DE DATOS
# ------------------------------------------------------------------------------

message("   - Preparando datos de mestizaje y trayectorias...")

datos_mestizaje <- ancestria_integrada %>%
  mutate(
    pct_au_nat = pctificar(au_nativo),
    pct_au_eur = pctificar(au_europeo),
    pct_au_afr = pctificar(au_africano),
    
    pct_cx_nat = pctificar(cx_na),
    pct_cx_eur = pctificar(cx_eu),
    pct_cx_afr = pctificar(cx_af),
    
    origen_y  = clasificar_y_general(cr_y),
    origen_mt = clasificar_mt_general(mt_dna),
    
    delta_nativo = pct_cx_nat - pct_au_nat
  ) %>%
  filter(
    !is.na(pct_au_nat),
    !is.na(pct_cx_nat)
  )

write.csv(
  datos_mestizaje,
  here::here("outputs", "tablas", "Obj5_Datos_Mestizaje.csv"),
  row.names = FALSE
)

message(paste("   - Número de muestras analizadas:", nrow(datos_mestizaje)))

# ------------------------------------------------------------------------------
# 2. ANÁLISIS ESTADÍSTICO INICIAL
# ------------------------------------------------------------------------------

message("   - Calculando pruebas estadísticas AU vs CX...")

cor_au_cx <- cor.test(
  datos_mestizaje$pct_au_nat,
  datos_mestizaje$pct_cx_nat,
  method = "spearman",
  exact = FALSE
)

capture.output(
  cor_au_cx,
  file = here::here("outputs", "tablas", "Obj5_Correlacion_AU_vs_CX.txt")
)

wilcox_au_cx <- wilcox.test(
  datos_mestizaje$pct_au_nat,
  datos_mestizaje$pct_cx_nat,
  paired = TRUE,
  exact = FALSE
)

capture.output(
  wilcox_au_cx,
  file = here::here("outputs", "tablas", "Obj5_Wilcoxon_AU_vs_CX.txt")
)

wilcox_delta <- wilcox.test(
  datos_mestizaje$delta_nativo,
  mu = 0,
  exact = FALSE
)

capture.output(
  wilcox_delta,
  file = here::here("outputs", "tablas", "Obj5_Wilcoxon_Delta_vs_0.txt")
)

# ------------------------------------------------------------------------------
# 3. TRAYECTORIAS GENÓMICAS INDIVIDUALES
# ------------------------------------------------------------------------------

message("   - Generando gráfico de trayectorias genómicas individuales...")

df_procesado <- datos_mestizaje %>%
  rowwise() %>%
  mutate(
    max_au = max(c(pct_au_afr, pct_au_eur, pct_au_nat), na.rm = TRUE),
    
    Origen_AU = case_when(
      max_au == pct_au_afr ~ "Africano",
      max_au == pct_au_eur ~ "Europeo",
      TRUE ~ "Nativo"
    ),
    
    max_cx = max(c(pct_cx_afr, pct_cx_eur, pct_cx_nat), na.rm = TRUE),
    
    Origen_CX = case_when(
      max_cx == pct_cx_afr ~ "Africano",
      max_cx == pct_cx_eur ~ "Europeo",
      TRUE ~ "Nativo"
    ),
    
    Origen_mt = if_else(origen_mt == "NATIVO", "Nativo", "No Nativo"),
    Origen_mt_plot = if_else(
      Origen_mt == "Nativo",
      "mtDNA Nativo",
      "mtDNA No Nativo"
    ),
    
    Origen_Y = clasificar_y_plot(cr_y)
  ) %>%
  ungroup()

df_long_tray <- df_procesado %>%
  mutate(
    Pct_mt = 50,
    Pct_Y  = 50,
    Pct_AU = max_au,
    Pct_CX = max_cx
  ) %>%
  select(
    id,
    Origen_mt_plot,
    Origen_mt,
    Origen_AU,
    Origen_CX,
    Origen_Y,
    Pct_mt,
    Pct_AU,
    Pct_CX,
    Pct_Y
  ) %>%
  pivot_longer(
    cols = matches("^(Origen|Pct)_(mt|AU|CX|Y)$"),
    names_to = c(".value", "Paso_Code"),
    names_sep = "_"
  ) %>%
  mutate(
    Origen_Paso = factor(
      Origen,
      levels = c("Africano", "Europeo", "Nativo", "No Nativo")
    ),
    
    Y_Plot = case_when(
      Paso_Code == "mt" & Origen_Paso == "Nativo" ~ 130,
      Paso_Code == "mt" & Origen_Paso == "No Nativo" ~ 280,
      
      Paso_Code == "Y" & Origen_Paso == "Nativo" ~ 50,
      Paso_Code == "Y" & Origen_Paso == "Europeo" ~ 200,
      Paso_Code == "Y" & Origen_Paso == "Africano" ~ 350,
      
      Origen_Paso == "Nativo" ~ Pct + runif(n(), -3, 3),
      Origen_Paso %in% c("Europeo", "No Nativo") ~ 150 + Pct + runif(n(), -3, 3),
      Origen_Paso == "Africano" ~ 300 + Pct + runif(n(), -3, 3),
      
      TRUE ~ NA_real_
    ),
    
    Paso_Num = case_when(
      Paso_Code == "mt" ~ 1,
      Paso_Code == "AU" ~ 2,
      Paso_Code == "CX" ~ 3,
      Paso_Code == "Y"  ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(Paso_Num),
    !is.na(Y_Plot)
  )

df_labels_izq <- tibble(
  x = 0.85,
  y = c(280, 130),
  label = c("mtDNA No Nativo", "mtDNA Nativo")
)

df_labels_der <- tibble(
  x = 4.15,
  y = c(350, 200, 50),
  label = c(
    "Hg Y Africano (A, B, E)",
    "Hg Y Europeo / No nativo",
    "Hg Y Nativo (Q)"
  )
)

p_trayectorias <- ggplot() +
  annotate(
    "rect",
    xmin = 0.5,
    xmax = 4.5,
    ymin = -5,
    ymax = 115,
    fill = "#e41a1c",
    alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = 0.5,
    xmax = 4.5,
    ymin = 145,
    ymax = 265,
    fill = "#377eb8",
    alpha = 0.08
  ) +
  annotate(
    "rect",
    xmin = 0.5,
    xmax = 4.5,
    ymin = 295,
    ymax = 415,
    fill = "#f0ad4e",
    alpha = 0.08
  ) +
  geom_vline(
    xintercept = c(1, 2, 3, 4),
    color = "grey70",
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_line(
    data = df_long_tray,
    aes(
      x = Paso_Num,
      y = Y_Plot,
      group = id,
      color = Origen_mt_plot
    ),
    alpha = 0.6,
    linewidth = 0.8
  ) +
  geom_text(
    data = df_labels_izq,
    aes(x = x, y = y, label = label),
    hjust = 1,
    fontface = "bold",
    size = 4,
    color = "black"
  ) +
  geom_text(
    data = df_labels_der,
    aes(x = x, y = y, label = label),
    hjust = 0,
    fontface = "bold",
    size = 4,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "mtDNA Nativo" = "#e41a1c",
      "mtDNA No Nativo" = "grey50"
    )
  ) +
  scale_x_continuous(
    breaks = 1:4,
    labels = c("ADNmt", "Autosómico", "Cromosoma X", "Hg Y Paterno"),
    expand = expansion(add = c(0.8, 0.8))
  ) +
  labs(
    title = "Trayectorias Genómicas Individuales",
    subtitle = "Líneas coloreadas por origen materno. Franjas coloreadas por origen paterno.",
    x = "Sistema Analizado",
    y = "Dominancia de Ancestría",
    color = "Origen Materno"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    plot.title = element_text(face = "bold")
  )

ggsave(
  here::here("outputs", "figuras", "Obj5_Trayectorias_Genomicas.png"),
  p_trayectorias,
  width = 28,
  height = 18,
  units = "cm",
  dpi = 300
)

# ------------------------------------------------------------------------------
# 4. GRÁFICO 1: COMPOSICIÓN GLOBAL PROMEDIO AU VS CX
# ------------------------------------------------------------------------------

message("     > Generando Gráfico 1: Composición global promedio con %...")

df_comp_parental <- datos_mestizaje %>%
  summarise(
    AU_Nativo   = mean(pct_au_nat, na.rm = TRUE),
    AU_Europeo  = mean(pct_au_eur, na.rm = TRUE),
    AU_Africano = mean(pct_au_afr, na.rm = TRUE),
    CX_Nativo   = mean(pct_cx_nat, na.rm = TRUE),
    CX_Europeo  = mean(pct_cx_eur, na.rm = TRUE),
    CX_Africano = mean(pct_cx_afr, na.rm = TRUE)
  ) %>%
  pivot_longer(
    everything(),
    names_to = c("Sistema", "Ancestria"),
    names_sep = "_"
  ) %>%
  mutate(
    Sistema = recode(
      Sistema,
      "AU" = "Autosómico",
      "CX" = "Cromosoma X"
    ),
    Ancestria = factor(
      Ancestria,
      levels = c("Africano", "Europeo", "Nativo")
    ),
    Etiqueta = label_pct(value)
  )

p_comp_parental <- ggplot(
  df_comp_parental,
  aes(x = Sistema, y = value, fill = Ancestria)
) +
  geom_col(width = 0.5, linewidth = 0.2) +
  geom_text(
    data = df_comp_parental %>% filter(value >= 0.03), # Si tus datos van de 0 a 100, cambiá 0.03 por 3
    aes(label = Etiqueta, color = Ancestria),
    position = position_stack(vjust = 0.5),
    fontface = "bold",
    size = 4.2
  ) +
  geom_text(
    data = df_comp_parental %>% filter(value < 3), # Si tus datos van de 0 a 100, cambiá 0.03 por 3
    aes(label = Etiqueta),
    position = position_stack(vjust = 0.5),
    vjust = 0, 
    color = "black",
    fontface = "bold",
    size = 4.2
  ) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = paleta_completa) +
  scale_color_manual(values = paleta_texto, guide = "none") +
  labs(
    title = "Composición Parental Promedio por Sistema Genético",
    x = "Sistema",
    y = "Proporción Promedio (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  here::here("outputs", "figuras", "Obj5_G1_Comp_Total.png"),
  p_comp_parental,
  width = 18,
  height = 14,
  units = "cm",
  dpi = 300
)

# ------------------------------------------------------------------------------
# 5. GRÁFICO 2: BARRAS INDIVIDUALES CON HAPLOGRUPOS REALES
# ------------------------------------------------------------------------------

message("     > Generando Gráfico 2: Perfil individual con haplogrupos...")

df_barras_ind <- datos_mestizaje %>%
  arrange(id) %>%
  mutate(
    id = factor(id, levels = unique(id)),
    haplo_mt = coalesce(as.character(mt_dna), "S/D"),
    haplo_y  = coalesce(as.character(cr_y), "S/D")
  )

df_mt_g2 <- df_barras_ind %>%
  select(id, Valor = haplo_mt) %>%
  mutate(
    Sistema = "1_ADNmt",
    Pct = 100,
    Fill = "Haplogrupo"
  )

df_y_g2 <- df_barras_ind %>%
  select(id, Valor = haplo_y) %>%
  mutate(
    Sistema = "4_Hg Y",
    Pct = 100,
    Fill = "Haplogrupo"
  )

df_au_g2 <- df_barras_ind %>%
  select(
    id,
    Nativo = pct_au_nat,
    Europeo = pct_au_eur,
    Africano = pct_au_afr
  ) %>%
  pivot_longer(
    -id,
    names_to = "Fill",
    values_to = "Pct"
  ) %>%
  mutate(
    Sistema = "2_Autosómico",
    Valor = ""
  )

df_cx_g2 <- df_barras_ind %>%
  select(
    id,
    Nativo = pct_cx_nat,
    Europeo = pct_cx_eur,
    Africano = pct_cx_afr
  ) %>%
  pivot_longer(
    -id,
    names_to = "Fill",
    values_to = "Pct"
  ) %>%
  mutate(
    Sistema = "3_Cromosoma X",
    Valor = ""
  )

df_combo_g2 <- bind_rows(
  df_mt_g2,
  df_au_g2,
  df_cx_g2,
  df_y_g2
) %>%
  mutate(
    Pct = as.numeric(Pct),
    Fill = factor(
      Fill,
      levels = c("Haplogrupo", "Africano", "Europeo", "Nativo")
    ),
    Etiqueta = case_when(
      Fill == "Haplogrupo" ~ label_haplo(Valor),
      Pct >= UMBRAL_ETIQUETA_PCT ~ paste0(round(Pct, 1), "%"),
      TRUE ~ ""
    )
  )

p_ind_bars <- ggplot(
  df_combo_g2,
  aes(x = Pct, y = id, fill = Fill)
) +
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.2
  ) +
  geom_text(
    data = df_combo_g2 %>% filter(Etiqueta != ""),
    aes(label = Etiqueta, color = Fill),
    position = position_stack(vjust = 0.5),
    fontface = "bold",
    size = 2.5,
    check_overlap = TRUE
  ) +
  scale_fill_manual(values = paleta_completa) +
  scale_color_manual(values = paleta_texto, guide = "none") +
  facet_grid(
    . ~ Sistema,
    scales = "free_x",
    labeller = labeller(
      Sistema = c(
        "1_ADNmt" = "ADNmt",
        "2_Autosómico" = "Autosomas",
        "3_Cromosoma X" = "Cromosoma X",
        "4_Hg Y" = "Hg Y"
      )
    )
  ) +
  labs(
    title = "Composición Genómica Individual Completa",
    x = "Porcentaje / Haplogrupo Real",
    y = "Muestra (ID)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 5.5),
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

ggsave(
  here::here("outputs", "figuras", "Obj5_G2_Barras_Individuales.png"),
  p_ind_bars,
  width = 25,
  height = 35,
  units = "cm",
  dpi = 300
)

# ------------------------------------------------------------------------------
# 6. GRÁFICO 3: BARRAS CATEGORIZADAS
# ------------------------------------------------------------------------------

message("     > Generando Gráfico 3: Perfil categorizado...")

df_combo_cat <- df_combo_g2 %>%
  left_join(
    df_procesado %>%
      select(id, origen_mt, origen_y, Origen_AU, Origen_CX),
    by = "id"
  ) %>%
  mutate(
    Categoria = case_when(
      as.character(id) == "57" ~ "1. Mayoritariamente Nativo",
      
      origen_y == "NATIVO" &
        origen_mt == "NATIVO" &
        Origen_AU == "Nativo" &
        Origen_CX == "Nativo" ~ "1. Mayoritariamente Nativo",
      
      origen_y == "NO_NATIVO" &
        origen_mt == "NO_NATIVO" &
        Origen_AU == "Europeo" &
        Origen_CX == "Europeo" ~ "3. Mayoritariamente Europeo",
      
      TRUE ~ "2. Patrón de Mezcla"
    ),
    
    Categoria = factor(
      Categoria,
      levels = c(
        "1. Mayoritariamente Nativo",
        "2. Patrón de Mezcla",
        "3. Mayoritariamente Europeo"
      )
    )
  )

p_ind_cat <- ggplot(
  df_combo_cat,
  aes(x = Pct, y = id, fill = Fill)
) +
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.2
  ) +
  geom_text(
    data = df_combo_cat %>% filter(Etiqueta != ""),
    aes(label = Etiqueta, color = Fill),
    position = position_stack(vjust = 0.5),
    fontface = "bold",
    size = 2.5,
    check_overlap = TRUE
  ) +
  facet_grid(
    Categoria ~ Sistema,
    scales = "free",
    space = "free_y",
    labeller = labeller(
      Sistema = c(
        "1_ADNmt" = "ADNmt",
        "2_Autosómico" = "Autosomas",
        "3_Cromosoma X" = "Cromosoma X",
        "4_Hg Y" = "Hg Y"
      )
    )
  ) +
  scale_fill_manual(values = paleta_completa) +
  scale_color_manual(values = paleta_texto, guide = "none") +
  labs(
    title = "Estructura Individual según Contribución Mayoritaria",
    x = "Porcentaje / Haplogrupo Real",
    y = "Muestra (ID)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 5.5),
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(angle = 0, face = "bold")
  )

ggsave(
  here::here("outputs", "figuras", "Obj5_G3_Barras_Categorizadas.png"),
  p_ind_cat,
  width = 28,
  height = 35,
  units = "cm",
  dpi = 300
)

# ------------------------------------------------------------------------------
# 7. GRÁFICO 4: MODELO TEÓRICO CON MADRE Y PADRE P0
# ------------------------------------------------------------------------------

message("     > Generando Gráfico 4: Modelo teórico con padres iniciales...")

df_teorico_parents <- tibble(
  Etiqueta = c("P0 Madre", "P0 Padre"),
  AU_Nativo = c(100, 0),
  AU_Europeo = c(0, 100),
  CX_Nativo = c(100, 0),
  CX_Europeo = c(0, 100)
)

df_teorico_filiales <- tibble(
  Generacion = 1:8
) %>%
  mutate(
    Etiqueta = paste0("F", Generacion),
    AU_Nativo  = (1 / 2)^Generacion * 100,
    AU_Europeo = 100 - AU_Nativo,
    CX_Nativo  = (1 / 2)^(Generacion - 1) * 100,
    CX_Europeo = 100 - CX_Nativo
  )

df_teo_long <- bind_rows(
  df_teorico_parents,
  df_teorico_filiales %>% select(-Generacion)
) %>%
  mutate(
    Etiqueta = factor(
      Etiqueta,
      levels = c("P0 Madre", "P0 Padre", paste0("F", 1:8))
    )
  ) %>%
  pivot_longer(
    -Etiqueta,
    names_to = c("Sistema", "Ancestria"),
    names_sep = "_"
  ) %>%
  mutate(
    Ancestria = factor(
      Ancestria,
      levels = c("Europeo", "Nativo")
    ),
    Sistema = recode(
      Sistema,
      "AU" = "Autosómico",
      "CX" = "Cromosoma X"
    ),
    Etiqueta_Pct = label_pct(value)
  )

p_modelo <- ggplot(
  df_teo_long,
  aes(x = Sistema, y = value, fill = Ancestria)
) +
  geom_col(
    position = "stack",
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  facet_grid(. ~ Etiqueta) +
  scale_fill_manual(values = paleta_completa) +
  geom_text(
    data = df_teo_long %>% filter(Etiqueta_Pct != ""),
    aes(label = Etiqueta_Pct, color = Ancestria),
    position = position_stack(vjust = 0.5),
    fontface = "bold",
    size = 2.5,
    check_overlap = TRUE
  ) +
  scale_color_manual(values = paleta_texto, guide = "none") +
  labs(
    title = "Modelo Teórico de Retrocruza Asimétrica Continua",
    subtitle = "Madre inicial 100% nativa y padre inicial 100% europeo, seguido de cruza con varones europeos",
    x = "Sistema Genético",
    y = "Proporción Esperada (%)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey20"),
    strip.text = element_text(color = "white", face = "bold", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  here::here("outputs", "figuras", "Obj5_G4_Modelo_Teorico.png"),
  p_modelo,
  width = 28,
  height = 13,
  units = "cm",
  dpi = 300
)

# ------------------------------------------------------------------------------
# 8. GRÁFICO 5: AGRUPACIÓN EMPÍRICA VS MODELO TEÓRICO
# ------------------------------------------------------------------------------

message("     > Generando Gráfico 5: Similitud Empírica vs Teórica filtrada al 10%...")

df_teo_calc <- df_teorico_filiales %>%
  select(
    F_Name = Etiqueta,
    AU_Teo = AU_Nativo,
    CX_Teo = CX_Nativo
  )

muestras_distancia <- datos_mestizaje %>%
  filter(
    origen_mt == "NATIVO",
    origen_y == "NO_NATIVO"
  ) %>%
  mutate(.join_key = 1) %>%
  inner_join(
    df_teo_calc %>% mutate(.join_key = 1),
    by = ".join_key"
  ) %>%
  select(-.join_key) %>%
  filter(
    abs(pct_au_nat - AU_Teo) <= 10,
    abs(pct_cx_nat - CX_Teo) <= 10
  ) %>%
  mutate(
    Distancia = sqrt(
      (pct_au_nat - AU_Teo)^2 +
        (pct_cx_nat - CX_Teo)^2
    )
  ) %>%
  group_by(id) %>%
  slice_min(Distancia, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    id,
    F_Name,
    pct_au_nat,
    pct_au_eur,
    pct_au_afr,
    pct_cx_nat,
    pct_cx_eur,
    pct_cx_afr,
    Distancia
  )

write.csv(
  muestras_distancia,
  here::here("outputs", "tablas", "Obj5_Muestras_Agrupadas_Modelo_Teorico.csv"),
  row.names = FALSE
)

if (nrow(muestras_distancia) > 0) {
  
  df_similitud <- bind_rows(
    df_teorico_filiales %>%
      select(
        F_Name = Etiqueta,
        AU_Nativo,
        AU_Europeo,
        CX_Nativo,
        CX_Europeo
      ) %>%
      mutate(
        id = paste0("[ REF ", F_Name, " ]"),
        Es_Ref = TRUE,
        AU_Africano = 0,
        CX_Africano = 0
      ) %>%
      rename(
        pct_au_nat = AU_Nativo,
        pct_au_eur = AU_Europeo,
        pct_cx_nat = CX_Nativo,
        pct_cx_eur = CX_Europeo,
        pct_au_afr = AU_Africano,
        pct_cx_afr = CX_Africano
      ),
    
    muestras_distancia %>%
      select(-Distancia) %>%
      mutate(Es_Ref = FALSE)
  ) %>%
    pivot_longer(
      cols = starts_with("pct_"),
      names_to = c("Sistema", "Ancestria"),
      names_pattern = "pct_(.*)_(.*)"
    ) %>%
    filter(F_Name %in% unique(muestras_distancia$F_Name)) %>%
    mutate(
      Ancestria = case_when(
        Ancestria == "nat" ~ "Nativo",
        Ancestria == "eur" ~ "Europeo",
        TRUE ~ "Africano"
      ),
      Ancestria = factor(
        Ancestria,
        levels = c("Africano", "Europeo", "Nativo")
      ),
      Sistema = if_else(
        Sistema == "au",
        "Autosómico",
        "Cromosoma X"
      ),
      Etiqueta = label_pct(value)
    ) %>%
    arrange(F_Name, desc(Es_Ref), id) %>%
    mutate(
      id = factor(id, levels = unique(id))
    )
  
  p_similitud <- ggplot(
    df_similitud,
    aes(x = value, y = id, fill = Ancestria)
  ) +
    geom_col(
      position = "stack",
      width = 0.8,
      color = "black",
      linewidth = 0.2
    ) +
    geom_text(
      data = df_similitud %>% filter(Etiqueta != ""),
      aes(label = Etiqueta, color = Ancestria),
      position = position_stack(vjust = 0.5),
      fontface = "bold",
      size = 2.4,
      check_overlap = TRUE
    ) +
    facet_grid(
      F_Name ~ Sistema,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_fill_manual(values = paleta_completa) +
    scale_color_manual(values = paleta_texto, guide = "none") +
    labs(
      title = "Alineación Empírica al Modelo Teórico de Retrocruza (Criterio ±10%)",
      subtitle = "Solo se muestran las muestras agrupadas matemáticamente a su filial de referencia",
      x = "Proporción (%)",
      y = "Muestra / Referencia"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 7),
      strip.text.x = element_text(face = "bold"),
      strip.text.y = element_text(face = "bold", angle = 0),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(
    here::here("outputs", "figuras", "Obj5_G5_Agrupacion_Teorica.png"),
    p_similitud,
    width = 25,
    height = 30,
    units = "cm",
    dpi = 300
  )
  
} else {
  
  message("     > No hubo muestras que cumplieran el criterio ±10%. Se omite el Gráfico 5.")
  
}

# ------------------------------------------------------------------------------
# 1. PREPARACIÓN DE DATOS Y VARIABLES (INCLUYENDO ADNmt)
# ------------------------------------------------------------------------------

# Tomamos la base integrada y preparamos las variables para los test
datos_estadistica <- ancestria_integrada %>%
  mutate(
    # Aseguramos que las variables de ancestría sean numéricas
    au_nativo = as.numeric(au_nativo),
    cx_nativo = as.numeric(cx_na),
    
    # Calculamos el Delta (Retención diferencial en el cromosoma X)
    delta_nativo = cx_nativo - au_nativo,
    
    # Clasificación Binaria del Linaje Paterno (Cromosoma Y)
    # Asumimos que los nativos empiezan con "Q"
    origen_y = case_when(
      is.na(cr_y) ~ NA_character_,
      str_starts(cr_y, "^Q") ~ "NATIVO",
      TRUE ~ "ALOCTONO"
    ),
    
    # Clasificación Binaria del Linaje Materno (ADNmt)
    # Los haplogrupos nativos americanos fundadores son A, B, C, D
    origen_mt = case_when(
      is.na(mt_dna) ~ NA_character_,
      str_starts(mt_dna, "^[ABCD]") ~ "NATIVO",
      TRUE ~ "ALOCTONO"
    )
  ) %>%
  # Filtramos solo aquellos individuos con información genómica completa
  filter(!is.na(au_nativo) & !is.na(cx_nativo))

message(paste("Muestras válidas para el análisis:", nrow(datos_estadistica)))

# ------------------------------------------------------------------------------
# 2. TEST 1: CORRELACIÓN ENTRE ANCESTRÍA AUTOSÓMICA Y CROMOSOMA X
# ------------------------------------------------------------------------------
message("\n--- TEST 1: Correlación de Spearman (AU_na vs CX_na) ---")

test_spearman <- cor.test(
  datos_estadistica$au_nativo, 
  datos_estadistica$cx_nativo, 
  method = "spearman"
)
print(test_spearman)

# ------------------------------------------------------------------------------
# 3. TEST 2: DIFERENCIA DE MAGNITUD ENTRE SISTEMAS (MISMOS INDIVIDUOS)
# ------------------------------------------------------------------------------
message("\n--- TEST 2: Prueba de Wilcoxon Emparejada (CX_na > AU_na) ---")

test_wilcox_paired <- wilcox.test(
  datos_estadistica$cx_nativo, 
  datos_estadistica$au_nativo, 
  paired = TRUE
)
print(test_wilcox_paired)

# ------------------------------------------------------------------------------
# 4. TEST 3: ANÁLISIS DEL DELTA DE RETENCIÓN DIFERENCIAL (Delta != 0)
# ------------------------------------------------------------------------------
message("\n--- TEST 3: Prueba de Wilcoxon de 1 Muestra (Delta = 0) ---")

test_delta_cero <- wilcox.test(
  datos_estadistica$delta_nativo, 
  mu = 0
)
print(test_delta_cero)

# ------------------------------------------------------------------------------
# 5. TEST 4: IMPACTO DEL LINAJE PATERNO EN EL DELTA (Q vs ALÓCTONOS)
# ------------------------------------------------------------------------------
message("\n--- TEST 4: Prueba de U de Mann-Whitney (Delta según Origen Y) ---")

# Filtramos los NA en el origen Y para el test
datos_y_validos <- datos_estadistica %>% filter(!is.na(origen_y))

test_y_delta <- wilcox.test(
  delta_nativo ~ origen_y, 
  data = datos_y_validos
)
print(test_y_delta)

# ------------------------------------------------------------------------------
# 6. TEST 5: IMPACTO DEL LINAJE MATERNO EN EL DELTA (ABCD vs ALÓCTONOS)
# ------------------------------------------------------------------------------
message("\n--- TEST 5: Prueba de U de Mann-Whitney (Delta según Origen ADNmt) ---")

# Filtramos los NA en el origen mtDNA para el test
datos_mt_validos <- datos_estadistica %>% filter(!is.na(origen_mt))

# Solo se corre si hay variabilidad (ej. si no son TODOS nativos)
if(length(unique(datos_mt_validos$origen_mt)) > 1) {
  test_mt_delta <- wilcox.test(
    delta_nativo ~ origen_mt, 
    data = datos_mt_validos
  )
  print(test_mt_delta)
} else {
  message("No hay suficiente varianza en origen_mt (quizás todos son linajes nativos) para correr el test.")
}

# ------------------------------------------------------------------------------
# 10. GRÁFICOS COMPLEMENTARIOS
# ------------------------------------------------------------------------------

message("   - Generando scatter plot de correlación y boxplot de deltas...")

graf_cor <- ggplot(
  datos_mestizaje,
  aes(x = pct_au_nat, y = pct_cx_nat)
) +
  geom_point(
    size = 3,
    alpha = 0.8,
    colour = "#377eb8"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "#e41a1c",
    linewidth = 1
  ) +
  labs(
    title = "Relación entre ancestría nativa autosómica y de cromosoma X",
    x = "Ancestría nativa autosómica (%)",
    y = "Ancestría nativa cromosoma X (%)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  here::here("outputs", "figuras", "Obj5_Figura_Correlacion_AU_CX.png"),
  graf_cor,
  width = 18,
  height = 12,
  units = "cm",
  dpi = 300
)

graf_delta <- ggplot(
  datos_mestizaje,
  aes(x = origen_y, y = delta_nativo, fill = origen_y)
) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(
    width = 0.1,
    alpha = 0.7,
    size = 2
  ) +
  scale_fill_manual(
    values = c(
      "NATIVO" = "#e41a1c",
      "NO_NATIVO" = "#377eb8"
    )
  ) +
  labs(
    title = "Exceso de ancestría nativa retenida en cromosoma X",
    x = "Origen del cromosoma Y",
    y = "Delta nativo (CX - AU)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  here::here("outputs", "figuras", "Obj5_Figura_Delta_OrigenY.png"),
  graf_delta,
  width = 16,
  height = 12,
  units = "cm",
  dpi = 300
)

message(">> Objetivo 5 completado.")