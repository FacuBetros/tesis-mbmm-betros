# ==============================================================================
# 06 - OBJETIVO 5: MESTIZAJE Y TRAYECTORIAS GENÓMICAS
# Autor: Facundo Betros
# ==============================================================================

message(">> Ejecutando Objetivo 5...")

# --- Asegurar que existan los directorios de salida ---
dir.create(here::here("outputs", "tablas"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("outputs", "figuras"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. PREPARACIÓN Y LIMPIEZA DE DATOS
# ------------------------------------------------------------------------------

message("   - Preparando datos de mestizaje y trayectorias...")

datos_mestizaje <- ancestria_integrada %>%
  mutate(
    # Normalizamos porcentajes (0 a 100)
    pct_au_nat = if_else(clean_numeric(au_nativo) <= 1, clean_numeric(au_nativo) * 100, clean_numeric(au_nativo)),
    pct_au_eur = if_else(clean_numeric(au_europeo) <= 1, clean_numeric(au_europeo) * 100, clean_numeric(au_europeo)),
    pct_au_afr = if_else(clean_numeric(au_africano) <= 1, clean_numeric(au_africano) * 100, clean_numeric(au_africano)),
    
    pct_cx_nat = if_else(clean_numeric(cx_na) <= 1, clean_numeric(cx_na) * 100, clean_numeric(cx_na)),
    pct_cx_eur = if_else(clean_numeric(cx_eu) <= 1, clean_numeric(cx_eu) * 100, clean_numeric(cx_eu)),
    pct_cx_afr = if_else(clean_numeric(cx_af) <= 1, clean_numeric(cx_af) * 100, clean_numeric(cx_af)),
    
    # Clasificación de Linajes
    origen_y  = if_else(coalesce(as.character(cr_y), "") == "Q", "NATIVO", "NO_NATIVO"),
    origen_mt = if_else(str_detect(toupper(coalesce(as.character(mt_dna), "")), "^[ABCD]"), "NATIVO", "NO_NATIVO"),
    
    # Delta
    delta_nativo = pct_cx_nat - pct_au_nat
  ) %>%
  filter(!is.na(pct_au_nat), !is.na(pct_cx_nat))

write.csv(
  datos_mestizaje,
  here::here("outputs", "tablas", "Obj5_Datos_Mestizaje.csv"),
  row.names = FALSE
)

message(paste("   - Número de muestras analizadas:", nrow(datos_mestizaje)))

# Paletas unificadas
paleta_completa <- c(
  "Nativo"      = "#e41a1c", 
  "Europeo"     = "#377eb8", 
  "Africano"    = "#f0ad4e",
  "Haplogrupo"  = "white"
)

paleta_texto <- c(
  "Nativo"      = "white", 
  "Europeo"     = "white", 
  "Africano"    = "black",
  "Haplogrupo"  = "black"
)

# ------------------------------------------------------------------------------
# 2. ANÁLISIS ESTADÍSTICO INICIAL (Correlación y Wilcoxon)
# ------------------------------------------------------------------------------

message("   - Calculando pruebas estadísticas AU vs CX...")

cor_au_cx <- cor.test(datos_mestizaje$pct_au_nat, datos_mestizaje$pct_cx_nat, method = "spearman")
capture.output(cor_au_cx, file = here::here("outputs", "tablas", "Obj5_Correlacion_AU_vs_CX.txt"))

wilcox_au_cx <- wilcox.test(datos_mestizaje$pct_au_nat, datos_mestizaje$pct_cx_nat, paired = TRUE, exact = FALSE)
capture.output(wilcox_au_cx, file = here::here("outputs", "tablas", "Obj5_Wilcoxon_AU_vs_CX.txt"))

wilcox_delta <- wilcox.test(datos_mestizaje$delta_nativo, mu = 0, exact = FALSE)
capture.output(wilcox_delta, file = here::here("outputs", "tablas", "Obj5_Wilcoxon_Delta_vs_0.txt"))

# ------------------------------------------------------------------------------
# 3. TRAYECTORIAS GENÓMICAS INDIVIDUALES
# ------------------------------------------------------------------------------

message("   - Generando gráfico de trayectorias genómicas individuales...")

df_procesado <- datos_mestizaje %>% 
  rowwise() %>% 
  mutate(
    max_au = max(c(pct_au_afr, pct_au_eur, pct_au_nat), na.rm = TRUE), 
    Origen_AU = case_when(max_au == pct_au_afr ~ "Africano", max_au == pct_au_eur ~ "Europeo", TRUE ~ "Nativo"), 
    max_cx = max(c(pct_cx_afr, pct_cx_eur, pct_cx_nat), na.rm = TRUE), 
    Origen_CX = case_when(max_cx == pct_cx_afr ~ "Africano", max_cx == pct_cx_eur ~ "Europeo", TRUE ~ "Nativo"), 
    Origen_mt = if_else(origen_mt == "NATIVO", "Nativo", "No Nativo"),
    Origen_mt_plot = if_else(Origen_mt == "Nativo", "mtDNA Nativo", "mtDNA No Nativo"),
    Origen_Y = case_when(
      coalesce(as.character(cr_y), "") == "Q" ~ "Nativo",
      str_detect(toupper(coalesce(as.character(cr_y), "")), "^E|^A|^B") ~ "Africano",
      TRUE ~ "Europeo"
    )
  ) %>% 
  ungroup() 

df_long_tray <- df_procesado %>% 
  mutate(Pct_mt = 50, Pct_Y = 50, Pct_AU = max_au, Pct_CX = max_cx) %>% 
  select(id, Origen_mt_plot, Origen_mt, Origen_AU, Origen_CX, Origen_Y, Pct_mt, Pct_AU, Pct_CX, Pct_Y) %>% 
  pivot_longer(cols = matches("^(Origen|Pct)_(mt|AU|CX|Y)$"), names_to = c(".value", "Paso_Code"), names_sep = "_") %>% 
  mutate(
    Origen_Paso = factor(Origen, levels = c("Africano", "Europeo", "Nativo", "No Nativo")), 
    # Anclajes estratégicos para evitar confusiones de banda
    Y_Plot = case_when(
      Paso_Code == "mt" & Origen_Paso == "Nativo" ~ 130,          
      Paso_Code == "mt" & Origen_Paso == "No Nativo" ~ 280,       
      Paso_Code == "Y" & Origen_Paso == "Nativo" ~ 50,
      Paso_Code == "Y" & Origen_Paso == "Europeo" ~ 200,
      Paso_Code == "Y" & Origen_Paso == "Africano" ~ 350,
      Origen_Paso == "Nativo" ~ Pct + runif(n(), -3, 3),
      Origen_Paso %in% c("Europeo", "No Nativo") ~ 150 + Pct + runif(n(), -3, 3), 
      Origen_Paso == "Africano" ~ 300 + Pct + runif(n(), -3, 3),          
      TRUE ~ 150 + Pct
    ), 
    Paso_Num = case_when(Paso_Code == "mt" ~ 1, Paso_Code == "AU" ~ 2, Paso_Code == "CX" ~ 3, Paso_Code == "Y"  ~ 4)
  ) %>% filter(!is.na(Paso_Num), !is.na(Y_Plot))

# Etiquetas ancladas exactamente en los vacíos inter-bandas
df_labels_izq <- tibble(x = 0.85, y = c(280, 130), label = c("mtDNA No Nativo", "mtDNA Nativo"))
df_labels_der <- tibble(x = 4.15, y = c(350, 200, 50), label = c("Hg Y Africano (E, etc.)", "Hg Y Europeo (R, J, I, etc.)", "Hg Y Nativo (Q)")) 

p_trayectorias <- ggplot() + 
  # Franjas coloreadas según el destino del Hg Y
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -5, ymax = 115, fill = "#e41a1c", alpha = 0.08) +   
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 145, ymax = 265, fill = "#377eb8", alpha = 0.08) +  
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 295, ymax = 415, fill = "#f0ad4e", alpha = 0.08) +  
  
  geom_vline(xintercept = c(1, 2, 3, 4), color = "grey70", size = 1, linetype = "dashed") + 
  
  geom_line(data = df_long_tray, aes(x = Paso_Num, y = Y_Plot, group = id, color = Origen_mt_plot), alpha = 0.6, size = 0.8) + 
  
  geom_text(data = df_labels_izq, aes(x = x, y = y, label = label), hjust = 1, fontface = "bold", size = 4, color = "black") + 
  geom_text(data = df_labels_der, aes(x = x, y = y, label = label), hjust = 0, fontface = "bold", size = 4, color = "black") + 
  
  scale_color_manual(values = c("mtDNA Nativo" = "#e41a1c", "mtDNA No Nativo" = "grey50")) + 
  scale_x_continuous(breaks = 1:4, labels = c("ADNmt", "Autosómico", "Cromosoma X", "Hg Y Paterno"), expand = expansion(add = c(0.8, 0.8))) + 
  labs(title = "Trayectorias Genómicas Individuales", subtitle = "Líneas coloreadas por origen materno. Franjas coloreadas por origen paterno.",
       x = "Sistema Analizado", y = "Dominancia de Ancestría", color = "Origen Materno") + 
  theme_minimal(base_size = 14) + 
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold", color = "black", size = 12), plot.title = element_text(face = "bold"))

ggsave(here::here("outputs", "figuras", "Obj5_Trayectorias_Genomicas.png"), p_trayectorias, width = 28, height = 18, units = "cm", dpi = 300)


# ------------------------------------------------------------------------------
# GRAFICO 1: Composición de Autosomas y Cromosoma X Total (Con Porcentajes)
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
  pivot_longer(everything(), names_to = c("Sistema", "Ancestria"), names_sep = "_") %>%
  mutate(
    Sistema   = recode(Sistema, "AU" = "Autosómico", "CX" = "Cromosoma X"),
    Ancestria = factor(Ancestria, levels = c("Africano", "Europeo", "Nativo"))
  )

p_comp_parental <- ggplot(df_comp_parental, aes(x = Sistema, y = value, fill = Ancestria)) +
  geom_col(width = 0.5, color = "black", size = 0.2) +
  geom_text(aes(
    label = ifelse(round(value, 1) > 0, paste0(round(value, 1), "%"), ""),
    color = Ancestria
  ), position = position_stack(vjust = 0.5), fontface = "bold", size = 4.5) +
  scale_fill_manual(values = paleta_completa) +
  scale_color_manual(values = paleta_texto, guide = "none") +
  labs(title = "Composición Ancestral Promedio por Sistema Genético",
       subtitle = "Comparación entre Autosomas y Cromosoma X para el total de la muestra",
       x = "Sistema Genético", y = "Proporción Promedio (%)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(here::here("outputs", "figuras", "Obj5_G1_Comp_Total.png"), p_comp_parental, width = 18, height = 14, units = "cm", dpi = 300)


# ------------------------------------------------------------------------------
# GRAFICO 2: Barras Individuales con Haplogrupos Reales
# ------------------------------------------------------------------------------
message("     > Generando Gráfico 2: Perfil individual con haplogrupos...")

df_barras_ind <- datos_mestizaje %>%
  arrange(id) %>%
  mutate(id = factor(id, levels = unique(id)), haplo_mt = coalesce(as.character(mt_dna), "S/D"), haplo_y  = coalesce(as.character(cr_y), "S/D"))

df_mt_g2 <- df_barras_ind %>% select(id, Valor = haplo_mt) %>% mutate(Sistema = "1_ADNmt", Pct = 100, Fill = "Haplogrupo")
df_y_g2  <- df_barras_ind %>% select(id, Valor = haplo_y)  %>% mutate(Sistema = "4_Hg Y",  Pct = 100, Fill = "Haplogrupo")

df_au_g2 <- df_barras_ind %>% select(id, Nativo=pct_au_nat, Europeo=pct_au_eur, Africano=pct_au_afr) %>% 
  pivot_longer(-id, names_to = "Fill", values_to = "Pct") %>% mutate(Sistema = "2_Autosómico", Valor = "")

df_cx_g2 <- df_barras_ind %>% select(id, Nativo=pct_cx_nat, Europeo=pct_cx_eur, Africano=pct_cx_afr) %>% 
  pivot_longer(-id, names_to = "Fill", values_to = "Pct") %>% mutate(Sistema = "3_Cromosoma X", Valor = "")

df_combo_g2 <- bind_rows(df_mt_g2, df_au_g2, df_cx_g2, df_y_g2) %>%
  mutate(Fill = factor(Fill, levels = c("Haplogrupo", "Africano", "Europeo", "Nativo")))

p_ind_bars <- ggplot(df_combo_g2, aes(x = Pct, y = id, fill = Fill)) +
  geom_col(width = 0.8, color = ifelse(df_combo_g2$Fill == "Haplogrupo", "grey70", "black"), size = 0.2) +
  geom_text(aes(
    label = case_when(
      Fill == "Haplogrupo" ~ Valor,
      round(as.numeric(Pct), 1) > 0 ~ paste0(round(as.numeric(Pct), 1), "%"),
      TRUE ~ ""
    ),
    color = Fill
  ), position = position_stack(vjust = 0.5), fontface = "bold", size = 2.8) +
  scale_fill_manual(values = paleta_completa) +
  scale_color_manual(values = paleta_texto, guide = "none") +
  facet_grid(. ~ Sistema, scales = "free_x", labeller = labeller(Sistema = c("1_ADNmt"="ADNmt", "2_Autosómico"="Autosomas", "3_Cromosoma X"="Cromosoma X", "4_Hg Y"="Hg Y"))) +
  labs(title = "Composición Genómica Individual Completa", x = "Porcentaje / Haplogrupo Real", y = "Muestra (ID)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 5.5), panel.grid = element_blank(), legend.position = "bottom", strip.text = element_text(face="bold"))

ggsave(here::here("outputs", "figuras", "Obj5_G2_Barras_Individuales.png"), p_ind_bars, width = 25, height = 35, units = "cm", dpi = 300)


# ------------------------------------------------------------------------------
# GRAFICO 3: Barras Categorizadas (Muestra 57 en Nativo)
# ------------------------------------------------------------------------------
message("     > Generando Gráfico 3: Perfil categorizado con Muestra 57 corregida...")

df_combo_cat <- df_combo_g2 %>%
  left_join(df_procesado %>% select(id, origen_mt, origen_y, Origen_AU, Origen_CX), by = "id") %>%
  mutate(
    Categoria = case_when(
      str_detect(as.character(id), "57") ~ "1. Mayoritariamente Nativo", 
      origen_y == "NATIVO" & origen_mt == "NATIVO" & Origen_AU == "Nativo" & Origen_CX == "Nativo" ~ "1. Mayoritariamente Nativo",
      origen_y == "NO_NATIVO" & origen_mt == "NO_NATIVO" & Origen_AU == "Europeo" & Origen_CX == "Europeo" ~ "3. Mayoritariamente Europeo",
      TRUE ~ "2. Patrón de Mezcla"
    )
  )

p_ind_cat <- ggplot(df_combo_cat, aes(x = Pct, y = id, fill = Fill)) +
  geom_col(width = 0.8, color = ifelse(df_combo_cat$Fill == "Haplogrupo", "grey70", "black"), size = 0.2) +
  geom_text(aes(
    label = case_when(
      Fill == "Haplogrupo" ~ Valor,
      round(as.numeric(Pct), 1) > 0 ~ paste0(round(as.numeric(Pct), 1), "%"),
      TRUE ~ ""
    ),
    color = Fill
  ), position = position_stack(vjust = 0.5), fontface = "bold", size = 2.8) +
  facet_grid(Categoria ~ Sistema, scales = "free", space = "free_y", labeller = labeller(Sistema = c("1_ADNmt"="ADNmt", "2_Autosómico"="Autosomas", "3_Cromosoma X"="Cromosoma X", "4_Hg Y"="Hg Y"))) +
  scale_fill_manual(values = paleta_completa) +
  scale_color_manual(values = paleta_texto, guide = "none") +
  labs(title = "Estructura Individual según Contribución Mayoritaria", x = "Porcentaje / Haplogrupo Real", y = "Muestra (ID)") +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 5.5), panel.grid = element_blank(), legend.position = "bottom", strip.text.x = element_text(face="bold"), strip.text.y = element_text(angle=0, face="bold"))

ggsave(here::here("outputs", "figuras", "Obj5_G3_Barras_Categorizadas.png"), p_ind_cat, width = 28, height = 35, units = "cm", dpi = 300)


# ------------------------------------------------------------------------------
# GRAFICO 4: Modelo Teórico Incluyendo Madre (P0) y Padre (P0)
# ------------------------------------------------------------------------------
message("     > Generando Gráfico 4: Modelo teórico con padres iniciales...")

df_teorico_parents <- tibble(
  Etiqueta  = c("P0 Madre", "P0 Padre"),
  AU_Nativo = c(100, 0), AU_Europeo= c(0, 100),
  CX_Nativo = c(100, 0), CX_Europeo= c(0, 100)
)

df_teorico_filiales <- tibble(Generacion = 1:8) %>%
  mutate(
    Etiqueta   = paste0("F", Generacion),
    AU_Nativo  = (1/2)^Generacion * 100,
    AU_Europeo = 100 - AU_Nativo,
    CX_Nativo  = (1/2)^(Generacion - 1) * 100,
    CX_Europeo = 100 - CX_Nativo
  )

df_teo_long <- bind_rows(df_teorico_parents, df_teorico_filiales %>% select(-Generacion)) %>%
  mutate(Etiqueta = factor(Etiqueta, levels = c("P0 Madre", "P0 Padre", paste0("F", 1:8)))) %>%
  pivot_longer(-Etiqueta, names_to = c("Sistema", "Ancestria"), names_sep = "_") %>%
  mutate(Ancestria = factor(Ancestria, levels = c("Europeo", "Nativo")), Sistema = recode(Sistema, "AU" = "Autosómico", "CX" = "Cromosoma X"))

p_modelo <- ggplot(df_teo_long, aes(x = Sistema, y = value, fill = Ancestria)) +
  geom_col(position = "stack", width = 0.7, color = "black", size = 0.2) +
  facet_grid(. ~ Etiqueta) +
  scale_fill_manual(values = paleta_completa) +
  geom_text(aes(
    label = ifelse(round(value, 1) > 0, paste0(round(value, 1), "%"), ""),
    color = Ancestria
  ), position = position_stack(vjust = 0.5), fontface = "bold", size = 2.8) +
  scale_color_manual(values = paleta_texto, guide = "none") +
  labs(title = "Modelo Teórico de Retrocruza Asimétrica Continua",
       subtitle = "Inicia con Madre 100% Nativa y Padre 100% Europeo, seguido de cruza con varones europeos (F1 a F8)",
       x = "Sistema Genético", y = "Proporción Esperada (%)") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_rect(fill = "grey20"), 
        strip.text = element_text(color = "white", face = "bold", size = 9), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here::here("outputs", "figuras", "Obj5_G4_Modelo_Teorico.png"), p_modelo, width = 28, height = 13, units = "cm", dpi = 300)


# ------------------------------------------------------------------------------
# GRAFICO 5: Agrupación Empírica vs Teórica en Barras Horizontales
# ------------------------------------------------------------------------------
message("     > Generando Gráfico 5: Similitud Empírica vs Teórica filtrada al 10%...")

df_teo_calc <- df_teorico_filiales %>% select(F_Name = Etiqueta, AU_Teo = AU_Nativo, CX_Teo = CX_Nativo)

# Sólo se agrupan las muestras que cumplan el criterio estricto de +-10%
muestras_distancia <- datos_mestizaje %>%
  filter(origen_mt == "NATIVO", origen_y == "NO_NATIVO") %>%
  cross_join(df_teo_calc) %>%
  filter(
    abs(pct_au_nat - AU_Teo) <= 10,
    abs(pct_cx_nat - CX_Teo) <= 10
  ) %>%
  mutate(Distancia = sqrt((pct_au_nat - AU_Teo)^2 + (pct_cx_nat - CX_Teo)^2)) %>%
  group_by(id) %>% 
  slice_min(Distancia, n = 1, with_ties = FALSE) %>% 
  ungroup() %>%
  select(id, F_Name, pct_au_nat, pct_au_eur, pct_au_afr, pct_cx_nat, pct_cx_eur, pct_cx_afr)

df_similitud <- bind_rows(
  df_teorico_filiales %>% select(F_Name = Etiqueta, AU_Nativo, AU_Europeo, CX_Nativo, CX_Europeo) %>% 
    mutate(id = paste0("[ REF ", F_Name, " ]"), Es_Ref = TRUE, AU_Africano=0, CX_Africano=0) %>%
    rename(pct_au_nat=AU_Nativo, pct_au_eur=AU_Europeo, pct_cx_nat=CX_Nativo, pct_cx_eur=CX_Europeo, pct_au_afr=AU_Africano, pct_cx_afr=CX_Africano),
  muestras_distancia %>% mutate(Es_Ref = FALSE)
) %>%
  pivot_longer(cols = starts_with("pct_"), names_to = c("Sistema", "Ancestria"), names_pattern = "pct_(.*)_(.*)") %>%
  filter(F_Name %in% unique(muestras_distancia$F_Name)) %>% # Mostrar solo las F con coincidencias
  mutate(
    Ancestria = case_when(Ancestria == "nat" ~ "Nativo", Ancestria == "eur" ~ "Europeo", TRUE ~ "Africano"),
    Ancestria = factor(Ancestria, levels = c("Africano", "Europeo", "Nativo")),
    Sistema   = if_else(Sistema == "au", "Autosómico", "Cromosoma X")
  ) %>%
  arrange(F_Name, desc(Es_Ref), id) %>% mutate(id = factor(id, levels = unique(id)))

p_similitud <- ggplot(df_similitud, aes(x = value, y = id, fill = Ancestria)) +
  geom_col(position = "stack", width = 0.8, color = "black", size = 0.2) +
  geom_text(aes(
    label = ifelse(round(value, 1) > 0, paste0(round(value, 1), "%"), ""),
    color = Ancestria
  ), position = position_stack(vjust = 0.5), fontface = "bold", size = 2.6) +
  facet_grid(F_Name ~ Sistema, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = paleta_completa) +
  scale_color_manual(values = paleta_texto, guide = "none") +
  labs(title = "Alineación Empírica al Modelo Teórico de Retrocruza (Criterio ±10%)",
       subtitle = "Solo se muestran las muestras agrupadas matemáticamente a su Filial de Referencia",
       x = "Proporción (%)", y = "Muestra / Referencia") +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 7, face = ifelse(str_detect(levels(df_similitud$id), "REF"), "bold", "plain")),
        strip.text.x = element_text(face = "bold"), strip.text.y = element_text(face = "bold", angle = 0),
        panel.grid = element_blank(), plot.title = element_text(face = "bold"))

ggsave(here::here("outputs", "figuras", "Obj5_G5_Agrupacion_Teorica.png"), p_similitud, width = 25, height = 30, units = "cm", dpi = 300)

# ------------------------------------------------------------------------------
# 6. ANÁLISIS ESTADÍSTICO DEL DELTA (Diferencia de Cero e Influencia Paterna)
# ------------------------------------------------------------------------------
message("   - Calculando estadísticos para el Delta de retención nativa...")

# Wilcoxon Delta vs 0
wilcox_delta <- wilcox.test(datos_mestizaje$delta_nativo, mu = 0, exact = FALSE)
capture.output(wilcox_delta, file = here::here("outputs", "tablas", "Obj5_Wilcoxon_Delta_vs_0.txt"))

# Wilcoxon Delta ~ Origen Paterno
if(length(unique(datos_mestizaje$origen_y)) > 1){
  wilcox_delta_y <- wilcox.test(delta_nativo ~ origen_y, data = datos_mestizaje, exact = FALSE)
  capture.output(wilcox_delta_y, file = here::here("outputs", "tablas", "Obj5_Wilcoxon_Delta_OrigenY.txt"))
} else {
  wilcox_delta_y <- NULL
}

# Tabla Resumen Objetivo 5
tabla_obj5 <- tibble::tibble(
  Analisis = c("Spearman AU_Nativo vs CX_Nativo", 
               "Wilcoxon pareado AU_Nativo vs CX_Nativo",
               "Wilcoxon Delta (CX-AU) vs 0", 
               "Wilcoxon Delta según origen paterno"),
  Estadistico = c(round(as.numeric(cor_au_cx$estimate), 3), 
                  round(as.numeric(wilcox_au_cx$statistic), 3),
                  round(as.numeric(wilcox_delta$statistic), 3),
                  ifelse(is.null(wilcox_delta_y), NA, round(as.numeric(wilcox_delta_y$statistic), 3))),
  p_valor = c(cor_au_cx$p.value, 
              wilcox_au_cx$p.value, 
              wilcox_delta$p.value,
              ifelse(is.null(wilcox_delta_y), NA, wilcox_delta_y$p.value))
)

write.csv(tabla_obj5, here::here("outputs", "tablas", "Obj5_Resumen_Estadistico.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. GRÁFICOS COMPLEMENTARIOS: CORRELACIÓN Y BOXPLOT
# ------------------------------------------------------------------------------
message("   - Generando scatter plot de correlación y boxplot de deltas...")

# Gráfico de Correlación
graf_cor <- ggplot(datos_mestizaje, aes(x = au_nativo, y = cx_nativo)) +
  geom_point(size = 3, alpha = 0.8, colour = "#377eb8") +
  geom_smooth(method = "lm", se = FALSE, colour = "#e41a1c", linewidth = 1) +
  labs(title = "Relación entre ancestría nativa autosómica y de cromosoma X",
       x = "Ancestría nativa autosómica", y = "Ancestría nativa cromosoma X") +
  theme_minimal(base_size = 12)

ggsave(here::here("outputs", "figuras", "Obj5_Figura_Correlacion_AU_CX.png"), 
       graf_cor, width = 18, height = 12, units = "cm", dpi = 300)

# Gráfico Boxplot del Delta
graf_delta <- ggplot(datos_mestizaje, aes(x = origen_y, y = delta_nativo, fill = origen_y)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.7) +
  scale_fill_manual(values = c("NATIVO" = "#e41a1c", "NO_NATIVO" = "#377eb8")) +
  labs(title = "Exceso de ancestría nativa retenida en cromosoma X",
       x = "Origen del cromosoma Y", y = "Δ Nativo (CX - AU)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(here::here("outputs", "figuras", "Obj5_Figura_Delta_OrigenY.png"), 
       graf_delta, width = 16, height = 12, units = "cm", dpi = 300)

message(">> Objetivo 5 completado.")