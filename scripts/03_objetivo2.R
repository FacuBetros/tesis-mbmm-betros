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
# 2. DIVERSIDAD GENÉTICA (Nei, Haplotipos únicos, Fisher)
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
# 3. ÍNDICE KAPPA (NevGen vs Whit Athey)
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
