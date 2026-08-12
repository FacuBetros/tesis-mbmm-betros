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
  here::here("outputs", "figuras", "Obj2_Frecuencia_Subhaplogrupos.png"),
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
# 4. ANÁLISIS DE COORDENADAS PRINCIPALES (PCoA) DE PERFILES Y-STR
# ------------------------------------------------------------------------------

message("   - Calculando PCoA a partir de los perfiles Y-STR...")

datos_pcoa <- muestras %>%
  dplyr::select(all_of(loci_y), cat_y) %>%
  mutate(
    across(
      all_of(loci_y),
      ~ suppressWarnings(as.numeric(.))
    )
  ) %>%
  tidyr::drop_na()

if (nrow(datos_pcoa) >= 3) {
  
  matriz_ystr <- datos_pcoa %>%
    dplyr::select(-cat_y) %>%
    mutate(across(everything(), as.numeric))
  
  distancias <- dist(
    as.matrix(matriz_ystr),
    method = "euclidean"
  )
  
  pcoa <- cmdscale(
    distancias,
    eig = TRUE,
    k = 2
  )
  
  df_pca <- data.frame(
    PC1 = pcoa$points[, 1],
    PC2 = pcoa$points[, 2],
    Haplogrupo = datos_pcoa$cat_y
  )
  
  write.csv(
    df_pca,
    here::here(
      "outputs",
      "tablas",
      "Obj2_PCoA_Coordenadas.csv"
    ),
    row.names = FALSE
  )
  
  p_pca2 <- ggplot(
    df_pca,
    aes(
      x = PC1,
      y = PC2,
      color = Haplogrupo
    )
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
      "Obj2_PCoA_YSTR_Elipses_Tesis.pdf"
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
  )
  
}

# ------------------------------------------------------------------------------
# 4. ÍNDICE KAPPA (NevGen vs Whit Athey)
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

message(">> Objetivo 2 completado.")
