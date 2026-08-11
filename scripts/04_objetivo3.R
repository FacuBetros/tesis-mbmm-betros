# ==============================================================================
# 04 - OBJETIVO 3: POBLACIONES COMPARATIVAS
# Autor: Facundo Betros
# ==============================================================================

message(">> Ejecutando Objetivo 3...")

# ------------------------------------------------------------------------------
# 1. GRÁFICO TERNARIO (Q - R - E)
# ------------------------------------------------------------------------------

message("   - Generando gráfico ternario...")

comparativas <- comparativas %>%
  mutate(
    q = clean_numeric(q),
    r = clean_numeric(r),
    e = clean_numeric(e)
  )

p_ternario <- ggtern(data = comparativas, aes(x = e, y = q, z = r)) +
  geom_point(aes(color = country), size = 3) +
  geom_point(
    data = comparativas %>% filter(region == "AMBA (Tesis)"),
    color = "black", shape = 18, size = 4
  ) +
  Tlab("Q (Nativo)") + Llab("E (Africano)") + Rlab("R (Europeo)") +
  labs(title = "Comparación de Haplogrupos Paternos", color = "País") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  here::here("outputs", "figuras", "Obj3_Haplogrupos_Ternario.png"),
  p_ternario, width = 18, height = 16, units = "cm", dpi = 300
)

# ------------------------------------------------------------------------------
# 2. MA-PLOT (Argentina vs Mundo usando YHRD)
# ------------------------------------------------------------------------------

message("   - Generando MA-Plot...")

yhrd_plot <- yhrd %>%
  left_join(muestras %>% select(id, cat_y), by = "id") %>%
  mutate(
    origen_linaje = case_when(
      cat_y == "Q" ~ "Q",
      cat_y %in% c("R", "E1B1B") ~ "R_E1B1B",
      cat_y == "E1B1A" ~ "E1B1A",
      TRUE ~ "OTROS"
    )
  )

# Ajuste de ceros para log2
min_m <- min(yhrd_plot$frec_mun[yhrd_plot$frec_mun > 0], na.rm = TRUE) / 2
min_a <- min(yhrd_plot$frec_arg[yhrd_plot$frec_arg > 0], na.rm = TRUE) / 2

yhrd_plot <- yhrd_plot %>%
  mutate(
    frec_mun_adj = ifelse(frec_mun == 0, min_m, frec_mun),
    frec_arg_adj = ifelse(frec_arg == 0, min_a, frec_arg),
    log2fc = log2(frec_arg_adj / frec_mun_adj)
  )

p_ma <- ggplot(yhrd_plot, aes(x = frec_mun_adj, y = log2fc, color = origen_linaje)) +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  geom_hline(yintercept = c(-3, 3), linetype = "dashed", color = "grey60") +
  geom_point(data = yhrd_plot %>% filter(origen_linaje == "OTROS"),
             alpha = 0.15, size = 2) +
  geom_point(data = yhrd_plot %>% filter(origen_linaje != "OTROS"),
             alpha = 0.9, size = 3.5) +
  scale_x_log10() +
  scale_color_manual(values = c(
    "Q" = "#e41a1c",
    "R_E1B1B" = "#377eb8",
    "E1B1A" = "#f0ad4e",
    "OTROS" = "grey70"
  )) +
  labs(
    title = "Frecuencia Relativa de Haplotipos Argentinos en YHRD",
    x = "Frecuencia Mundial",
    y = "Log2(Argentina / Mundo)",
    color = "Linaje"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(
  here::here("outputs", "figuras", "Obj3_MA_Plot.png"),
  p_ma, width = 28, height = 16, units = "cm", dpi = 300
)

# ------------------------------------------------------------------------------
# 3. CHI-CUADRADO (Q vs R × Frecuente vs Raro)
# ------------------------------------------------------------------------------

message("   - Ejecutando Chi-cuadrado...")

yhrd_chi <- yhrd_plot %>%
  filter(cat_y %in% c("Q", "R")) %>%
  mutate(
    Clasif_Haplogrupo = if_else(cat_y == "Q", "Haplogrupo Q", "Haplogrupo R"),
    Categoria_Frecuencia = if_else(
      frec_arg > median(frec_arg, na.rm = TRUE),
      "Frecuente",
      "No Frecuente / Raro"
    )
  )

tabla_contingencia_yhrd <- table(
  yhrd_chi$Clasif_Haplogrupo,
  yhrd_chi$Categoria_Frecuencia
)

write.csv(
  as.data.frame.matrix(tabla_contingencia_yhrd),
  here::here("outputs", "tablas", "Obj3_ChiCuadrado_Tabla.csv"),
  row.names = TRUE
)

resultado_chi_yhrd <- chisq.test(tabla_contingencia_yhrd)

capture.output(
  resultado_chi_yhrd,
  file = here::here("outputs", "tablas", "Obj3_ChiCuadrado_Resultado.txt")
)

message(">> Objetivo 3 completado.")
