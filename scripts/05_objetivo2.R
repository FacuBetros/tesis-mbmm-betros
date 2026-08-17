# ==============================================================================
# 03 - OBJETIVO 2: Predecir y caracterizar los haplogrupos del cromosoma Y a partir de los haplotipos obtenidos en la población analizada.
# Autor: Facundo Betros
# ==============================================================================

message(">> Ejecutando Objetivo 2...")

# ------------------------------------------------------------------------------
# 1. FRECUENCIA DE SUBHAPLOGRUPOS
# ------------------------------------------------------------------------------

message("   - Calculando frecuencias de haplogrupos...")

frecuencias_subhg <- muestras %>%
  count(cat_y, hg_y, sort = TRUE) %>%
  ungroup() %>%
  mutate(frecuencia = round(100 * n / sum(n), 2))

p_hg <- ggplot(frecuencias_subhg, aes(x = reorder(hg_y, n), y = n, fill = cat_y)) +
  geom_col(color = "black", width = 0.8) +
  geom_text(aes(label = paste0(n, " (", frecuencia, "%)")),
            hjust = -0.1, size = 3.7, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = paleta_hg) +
  labs(title = "Frecuencia de Subhaplogrupos del Cromosoma Y",
       x = "Subhaplogrupo", y = "Frecuencia") +
  expand_limits(y = max(frecuencias_subhg$n) * 1.25) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"))

ggsave(
  here::here("outputs", "figuras", "Obj2.1_Frecuencia_Subhaplogrupos.png"),
  plot = p_hg, width = 24, height = 14, units = "cm", dpi = 300
)

write.csv(
  frecuencias_subhg,
  here::here("outputs", "tablas", "Obj2_Frecuencias_Subhaplogrupos.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 2. CONCORDANCIA ENTRE PREDICTORES
# ------------------------------------------------------------------------------

message("   - Calculando concordancia NevGen vs Whit Athey...")

datos_concordancia_loci <- muestras %>%
  mutate(
    Tipo_Resolucion = case_when(
      Kit_Origen == "PowerPlex Y12" ~ "Perfiles de 12 loci",
      Kit_Origen == "Yfiler 17"     ~ "Perfiles de 17 loci",
      Kit_Origen == "PowerPlex Y23" ~ "Perfiles de 23 loci",
      TRUE                          ~ "Otros"
    ),
    Coincide = if_else(shg_nevgen == shg_hprg, "Concordante", "Discrepante")
  ) %>%
  filter(
    !is.na(shg_nevgen), !is.na(shg_hprg),
    shg_nevgen != "ND", shg_hprg != "ND",
    shg_nevgen != "NA", shg_hprg != "NA"
  )

tabla_concordancia_resolucion <- datos_concordancia_loci %>%
  count(Tipo_Resolucion, Coincide) %>%
  group_by(Tipo_Resolucion) %>%
  mutate(
    Total_Grupo = sum(n),
    Porcentaje = round(100 * n / Total_Grupo, 2)
  ) %>%
  ungroup()

write.csv(
  tabla_concordancia_resolucion,
  here::here("outputs", "tablas", "Obj2_Concordancia_Predictores.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 3. DIVERSIDAD GENÉTICA (Nei, Haplotipos únicos, Fisher)
# ------------------------------------------------------------------------------

message("   - Calculando diversidad genética...")

# Diversidad por haplogrupo
diversidad_hg <- muestras %>%
  group_by(cat_y) %>%
  summarise(
    N = n(),
    Haplotipos = n_distinct(haplotipo),
    HD = {
      frecuencias <- table(haplotipo) / n()
      n_ind <- n()
      if(n_ind > 1) (n_ind / (n_ind - 1)) * (1 - sum(frecuencias^2)) else 0
    }
  ) %>%
  arrange(desc(N))

write.csv(
  diversidad_hg,
  here::here("outputs", "tablas", "Obj2_Diversidad_Haplogrupos.csv"),
  row.names = FALSE
)

# Haplotipos únicos por haplogrupo
haplo_hg <- muestras %>%
  group_by(cat_y, haplotipo) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(unico = if_else(n == 1, 1, 0)) %>%
  group_by(cat_y) %>%
  summarise(
    Total = n(),
    Unicos = sum(unico),
    Prop_Unicos = Unicos / Total
  ) %>%
  arrange(desc(Prop_Unicos))

write.csv(
  haplo_hg,
  here::here("outputs", "tablas", "Obj2_Haplotipos_Unicos_por_HG.csv"),
  row.names = FALSE
)

# Fisher Q vs R
tabla_qr <- muestras %>%
  group_by(cat_y, haplotipo) %>%
  filter(cat_y %in% c("Q", "R")) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(categoria = if_else(n == 1, "Unico", "Repetido")) %>%
  count(cat_y, categoria) %>%
  pivot_wider(names_from = categoria, values_from = n, values_fill = 0)

mat_qr <- tabla_qr %>% tibble::column_to_rownames("cat_y") %>% as.matrix()
resultado_fisher <- fisher.test(mat_qr)

write.csv(
  tabla_qr,
  here::here("outputs", "tablas", "Obj2_Fisher_Q_vs_R_tabla.csv"),
  row.names = FALSE
)

capture.output(
  resultado_fisher,
  file = here::here("outputs", "tablas", "Obj2_Fisher_Q_vs_R_resultado.txt")
)

# ------------------------------------------------------------------------------
# HEATMAP 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 4. DISTANCIAS GENÉTICAS ENTRE HAPLOGRUPOS (FST Y SLATKIN)
# ------------------------------------------------------------------------------

message("   - Calculando FST pareado y distancia de Slatkin...")

# Mantener solamente haplogrupos con tamaño muestral suficiente
haplos_validos <- muestras %>%
  count(cat_y) %>%
  filter(
    n >= 5,
    cat_y !="NO", 
    cat_y !="ND"
    ) %>%
  pull(cat_y)

# Construcción de la matriz para FST
datos_fst <- muestras %>%
  filter(cat_y %in% haplos_validos) %>%
  select(cat_y, all_of(marcadores_minimos)) %>%
  mutate(
    across(
      all_of(marcadores_minimos),
      ~ suppressWarnings(as.numeric(.))
    )
  ) %>%
  mutate(
    across(
      all_of(marcadores_minimos),
      ~ tidyr::replace_na(., 0)
    )
  ) %>%
  as.data.frame()

message(
  "Muestras utilizadas para FST: ",
  nrow(datos_fst)
)

# Conversión para haploides
datos_haploides <- datos_fst

for(i in 2:ncol(datos_haploides)){
  
  datos_haploides[[i]] <- as.numeric(
    paste0(
      datos_haploides[[i]],
      datos_haploides[[i]]
    )
<<<<<<< HEAD:scripts/05_objetivo2.R
=======
  ) +
    geom_point(
      size = 3,
      alpha = 0.8
    ) +
    stat_ellipse(
      aes(fill = Haplogrupo),
      geom = "polygon",
      alpha = 0.15,
      color = NA,
      level = 0.95
    ) +
    scale_color_manual(values = paleta_hg) +
    scale_fill_manual(values = paleta_hg) +
    theme_bw(base_size = 12) +
    labs(
      title = "PCoA de perfiles Y-STR",
      x = "Coordenada Principal 1",
      y = "Coordenada Principal 2",
      color = "Haplogrupo",
      fill = "Haplogrupo"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
  
  ggsave(
    here::here(
      "outputs",
      "figuras",
      "Obj2.2_PCoA_YSTR_Elipses_Tesis.pdf"
    ),
    p_pca2,
    width = 10,
    height = 8,
    dpi = 600
  )
  
  message("   -> Gráfico PCoA guardado correctamente.")
  
} else {
  
  warning(
    "No hay suficientes perfiles completos para calcular el PCoA."
>>>>>>> 8f0b8e67cff1d400351dfb1cfdc1c6074f0900e3:scripts/03_objetivo2.R
  )
  
}

# Cálculo del FST
matriz_fst <- hierfstat::genet.dist(
  datos_haploides,
  method = "WC84"
)

matriz_fst_cuadrada <- as.matrix(matriz_fst)

matriz_fst_cuadrada[
  matriz_fst_cuadrada < 0
] <- 0

# Distancia linealizada de Slatkin
matriz_slatkin <-
  matriz_fst_cuadrada /
  (1 - matriz_fst_cuadrada)

# Exportar matrices
write.csv(
  matriz_fst_cuadrada,
  here::here(
    "outputs",
    "tablas",
    "Obj2_FST_Pareado.csv"
  )
)

write.csv(
  matriz_slatkin,
  here::here(
    "outputs",
    "tablas",
    "Obj2_FST_Slatkin.csv"
  )
)

# ------------------------------------------------------------------------------
# HEATMAP TRIPARTITO
# ------------------------------------------------------------------------------
hd_haplogrupos <- diversidad_hg$HD
names(hd_haplogrupos) <- diversidad_hg$cat_y

matriz_fusion <- matriz_fst_cuadrada

matriz_fusion[
  lower.tri(matriz_fusion)
] <- matriz_slatkin[
  lower.tri(matriz_slatkin)
]

for(hg in intersect(
  rownames(matriz_fusion),
  names(hd_haplogrupos)
)){
  
  matriz_fusion[hg, hg] <- hd_haplogrupos[hg]
  
}
df_heatmap <- as.data.frame(
  as.table(matriz_fusion)
)

colnames(df_heatmap) <- c(
  "Haplogrupo_Y",
  "Haplogrupo_X",
  "Valor"
)

orden_haplos <- sort(
  unique(
    as.character(df_heatmap$Haplogrupo_X)
  )
)

df_heatmap$Haplogrupo_Y <- factor(
  df_heatmap$Haplogrupo_Y,
  levels = rev(orden_haplos)
)

df_heatmap$Haplogrupo_X <- factor(
  df_heatmap$Haplogrupo_X,
  levels = orden_haplos
)

p_tripartito <- ggplot(
  df_heatmap,
  aes(
    x = Haplogrupo_X,
    y = Haplogrupo_Y,
    fill = Valor
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = ifelse(
        is.na(Valor),
        "-",
        sprintf("%.3f", Valor)
      )
    ),
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_gradientn(
    colors = c(
      "#4575b4",
      "white",
      "#d73027"
    ),
    na.value = "grey90",
    name = "Distancia"
  ) +
  labs(
    title = "Distancias Genéticas entre Haplogrupos Y",
    subtitle = "Diagonal: Diversidad intrahaplogrupo (HD de Nei) | Superior: FST | Inferior: FST linealizado de Slatkin",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

ggsave(
  here::here(
    "outputs",
    "figuras",
    "Obj2.2_Heatmap_FST_Slatkin.png"
  ),
  plot = p_tripartito,
  width = 10,
  height = 8,
  dpi = 300
)

message("   -> Heatmap FST/Slatkin guardado correctamente.")

# ------------------------------------------------------------------------------
# 5. ÍNDICE KAPPA (NevGen vs Whit Athey)
# ------------------------------------------------------------------------------

message("   - Calculando índice Kappa...")

datos_kappa <- muestras %>%
  select(shg_nevgen, shg_hprg) %>%
  filter(
    !is.na(shg_nevgen), !is.na(shg_hprg),
    shg_nevgen != "ND", shg_hprg != "ND",
    shg_nevgen != "NA", shg_hprg != "NA"
  )

resultado_kappa <- irr::kappa2(datos_kappa)

capture.output(
  resultado_kappa,
  file = here::here("outputs", "tablas", "Obj2_Kappa_NevGen_vs_HPRG.txt")
)

# ------------------------------------------------------------------------------
# 6. ANÁLISIS DE COORDENADAS PRINCIPALES (PCoA)
# ------------------------------------------------------------------------------

# El análisis PCoA fue realizado en GenAlEx utilizando el haplotipo mínimo
# y se presenta en la sección de resultados correspondiente.

message(">> Objetivo 2 completado.")
