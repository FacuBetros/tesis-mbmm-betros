# ==============================================================================
# 06 - OBJETIVO 5: MESTIZAJE Y TRAYECTORIAS GENÓMICAS
# Autor: Facundo Betros
# ==============================================================================

message(">> Ejecutando Objetivo 5...")

# ------------------------------------------------------------------------------
# 1. PREPARACIÓN DE DATOS
# ------------------------------------------------------------------------------

message("   - Preparando datos de mestizaje...")

datos_mestizaje <- ancestria_integrada %>%
  left_join(
    muestras %>% select(id, cat_y),
    by = "id"
  ) %>%
  mutate(
    au_nativo = clean_numeric(au_nativo),
    cx_nativo = clean_numeric(cx_nativo),
    
    origen_y = if_else(cat_y == "Q", "NATIVO", "NO_NATIVO"),
    
    delta_nativo = cx_nativo - au_nativo
  ) %>%
  filter(
    !is.na(au_nativo),
    !is.na(cx_nativo)
  )

write.csv(
  datos_mestizaje,
  here::here("outputs", "tablas", "Obj5_Datos_Mestizaje.csv"),
  row.names = FALSE
)

message(paste("   - Número de muestras analizadas:", nrow(datos_mestizaje)))

# ------------------------------------------------------------------------------
# 2. CORRELACIÓN ENTRE AUTOSOMAS Y CROMOSOMA X
# ------------------------------------------------------------------------------

message("   - Calculando correlación Spearman AU vs CX...")

cor_au_cx <- cor.test(
  datos_mestizaje$au_nativo,
  datos_mestizaje$cx_nativo,
  method = "spearman"
)

capture.output(
  cor_au_cx,
  file = here::here("outputs", "tablas", "Obj5_Correlacion_AU_vs_CX.txt")
)

# ------------------------------------------------------------------------------
# 3. WILCOXON PAREADO (AU vs CX)
# ------------------------------------------------------------------------------

message("   - Ejecutando Wilcoxon pareado AU vs CX...")

wilcox_au_cx <- wilcox.test(
  datos_mestizaje$au_nativo,
  datos_mestizaje$cx_nativo,
  paired = TRUE,
  exact = FALSE
)

capture.output(
  wilcox_au_cx,
  file = here::here("outputs", "tablas", "Obj5_Wilcoxon_AU_vs_CX.txt")
)

# ------------------------------------------------------------------------------
# 4. DELTA DE RETENCIÓN NATIVA
# ------------------------------------------------------------------------------

message("   - Calculando delta de ancestría nativa...")

summary_delta <- datos_mestizaje %>%
  summarise(
    Promedio = mean(delta_nativo),
    Mediana  = median(delta_nativo),
    SD       = sd(delta_nativo),
    Min      = min(delta_nativo),
    Max      = max(delta_nativo)
  )

write.csv(
  summary_delta,
  here::here("outputs", "tablas", "Obj5_Delta_Resumen.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. GRÁFICO ALUVIAL (Trayectorias Genómicas Individuales)
# ------------------------------------------------------------------------------

message("   - Generando gráfico aluvial...")

df_aluvial <- datos_mestizaje %>%
  mutate(
    Origen_AU = case_when(
      au_nativo >= 0.5 ~ "Nativo",
      TRUE ~ "No Nativo"
    ),
    Origen_CX = case_when(
      cx_nativo >= 0.5 ~ "Nativo",
      TRUE ~ "No Nativo"
    ),
    Origen_Y = if_else(cat_y == "Q", "Nativo", "No Nativo"),
    Frecuencia = 1
  ) %>%
  mutate(
    Origen_AU = factor(Origen_AU, levels = c("Nativo", "No Nativo")),
    Origen_CX = factor(Origen_CX, levels = c("Nativo", "No Nativo")),
    Origen_Y  = factor(Origen_Y,  levels = c("Nativo", "No Nativo"))
  )

p_aluvial <- ggplot(df_aluvial,
                    aes(y = Frecuencia,
                        axis1 = Origen_Y,
                        axis2 = Origen_AU,
                        axis3 = Origen_CX)) +
  geom_alluvium(aes(fill = Origen_Y), width = 1/8, alpha = 0.8, color = "white") +
  geom_stratum(width = 1/4, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Nativo" = "#e41a1c", "No Nativo" = "grey60")) +
  scale_x_discrete(limits = c("Y", "Autosomas", "Cromosoma X")) +
  labs(title = "Trayectoria Genómica Individual",
       y = NULL, fill = "Origen Y") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  here::here("outputs", "figuras", "Obj5_Aluvial_Individual.png"),
  p_aluvial, width = 25, height = 18, units = "cm", dpi = 300
)

# ------------------------------------------------------------------------------
# 6. TRAYECTORIAS PARALELAS (Parallel Coordinates)
# ------------------------------------------------------------------------------

message("   - Generando trayectorias paralelas...")

df_parallel <- datos_mestizaje %>%
  mutate(
    AU = au_nativo,
    CX = cx_nativo,
    Y  = if_else(cat_y == "Q", 1, 0)
  ) %>%
  select(id, AU, CX, Y)

p_parallel <- ggplot(df_parallel) +
  geom_line(aes(x = c("Y", "Autosomas", "Cromosoma X"),
                y = c(Y, AU, CX),
                group = id,
                color = AU),
            alpha = 0.5) +
  scale_color_gradient(low = "#f0ad4e", high = "#e41a1c") +
  labs(title = "Trayectorias Genómicas (Parallel Coordinates)",
       x = NULL, y = "Proporción") +
  theme_bw(base_size = 14)

ggsave(
  here::here("outputs", "figuras", "Obj5_Trayectorias_Paralelas.png"),
  p_parallel, width = 25, height = 18, units = "cm", dpi = 300
)

message(">> Objetivo 5 completado.")
