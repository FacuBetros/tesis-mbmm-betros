# ==============================================================================
# 02 - OBJETIVO 1: Identificar y caracterizar los haplotipos del cromosoma Y presentes en la muestra poblacional del Área Metropolitana de Buenos Aires.
# Autor: Facundo Betros
# ==============================================================================

message(">> Ejecutando Objetivo 1...")

# Crear directorios de salida si no existen (evita errores al clonar de GitHub)
dir.create(here::here("outputs", "tablas"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("outputs", "figuras"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. PARÁMETROS FORENSES POR KIT
# ------------------------------------------------------------------------------

message("   - Calculando parámetros forenses...")

# Definición de kits
marcadores_y12 <- c("dys19","dys389i","dys389ii","dys390","dys391","dys392","dys393",
                    "dys385a","dys385b","dys438","dys437","dys439")

marcadores_y17 <- c(marcadores_y12, "dys448","dys635","dys458","dys456","gatah4")

marcadores_y23 <- c(marcadores_y17, "dys576","dys481","dys549","dys533","dys570","dys643")

extra_y17 <- c(
  "dys448","dys635","dys458",
  "dys456","gatah4"
)

extra_y23 <- c(
  "dys576","dys481","dys549",
  "dys533","dys570","dys643"
)

# Clasificación de muestras según los marcadores observados
muestras <- muestras %>%
  rowwise() %>%
  mutate(
    
    tiene_y17 = any(
      !is.na(c_across(all_of(extra_y17))) &
        !(c_across(all_of(extra_y17)) %in% c("", "NA", "0", 0))
    ),
    
    tiene_y23 = any(
      !is.na(c_across(all_of(extra_y23))) &
        !(c_across(all_of(extra_y23)) %in% c("", "NA", "0", 0))
    ),
    
    sin_datos = all(
      is.na(c_across(all_of(loci_y))) |
        c_across(all_of(loci_y)) %in% c("", "NA")
    ),
    
    Kit_Origen = case_when(
      sin_datos ~ "Sin Datos",
      tiene_y23 ~ "PowerPlex Y23",
      tiene_y17 ~ "Yfiler 17",
      TRUE ~ "PowerPlex Y12"
    )
    
  ) %>%
  ungroup()

# Función para calcular parámetros forenses
calcular_forenses <- function(df, marcadores, nombre_kit) {
  df %>%
    filter(Kit_Origen == nombre_kit) %>%
    unite("haplotipo_kit", all_of(marcadores), sep = "-", remove = FALSE) %>%
    count(haplotipo_kit, name = "Conteo") %>%
    mutate(
      N_total = sum(Conteo),
      Frecuencia = Conteo / N_total,
      Frec_Cuadrado = Frecuencia^2
    ) %>%
    summarise(
      Kit = nombre_kit,
      N_Analizado = max(N_total),
      Ht_Distintos = n(),
      Ht_Unicos = sum(Conteo == 1),
      DC = Ht_Distintos / N_Analizado,
      MP = sum(Frec_Cuadrado),
      PD = 1 - MP,
      HD = ifelse(
        N_Analizado > 1,
        (N_Analizado / (N_Analizado - 1)) *
          (1 - sum(Frec_Cuadrado)),
        NA_real_
      )
    )
}

tabla_parametros_forenses <- bind_rows(
  calcular_forenses(muestras, marcadores_y12, "PowerPlex Y12"),
  calcular_forenses(muestras, marcadores_y17, "Yfiler 17"),
  calcular_forenses(muestras, marcadores_y23, "PowerPlex Y23")
)

write.csv(tabla_parametros_forenses,
          here::here("outputs", "tablas", "Obj1_Parametros_Forenses.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. HAPLOTIPOS ÚNICOS Y REPETIDOS
# ------------------------------------------------------------------------------

message("   - Analizando haplotipos únicos y repetidos...")

haplotipos_plot <- muestras %>%
  count(haplotipo, name = "repeticiones") %>%
  mutate(
    categoria = case_when(
      repeticiones == 1 ~ "Únicos",
      repeticiones == 2 ~ "2 Repeticiones",
      repeticiones == 3 ~ "3 Repeticiones",
      repeticiones >= 4 ~ "4+ Repeticiones"
    ),
    grupo = if_else(categoria == "Únicos", "Únicos", "Repetidos")
  ) %>%
  count(grupo, categoria, name = "n") %>%
  group_by(grupo) %>%
  mutate(porcentaje = n / sum(n)) %>%
  ungroup()

p_haplotipos <- ggplot(haplotipos_plot, aes(x = grupo, y = porcentaje, fill = categoria)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3) +
  scale_fill_manual(values = paleta_haplotipos) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Distribución de Haplotipos", x = NULL, y = "Porcentaje") +
  theme_bw(base_size = 12)

ggsave(
  here::here("outputs", "figuras", "Obj1.1_Haplotipos_Unicos_Repetidos.png"),
  p_haplotipos, width = 16, height = 12, units = "cm", dpi = 300
)

# ------------------------------------------------------------------------------
# 4. MICROVARIANTES
# ------------------------------------------------------------------------------

message("   - Procesando microvariantes...")

microvariantes_muestras <- muestras %>%
  select(id, all_of(loci_y)) %>%
  pivot_longer(cols = all_of(loci_y), names_to = "locus", values_to = "alelo") %>%
  filter(str_detect(alelo, "\\."))

microvariantes_resumen <- microvariantes_muestras %>%
  count(locus, alelo, sort = TRUE) %>%
  mutate(frecuencia = round(100 * n / sum(n), 2))

write.csv(
  microvariantes_resumen,
  here::here("outputs", "tablas", "Obj1_Microvariantes_Resumen.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. BUBBLE CHART DE MICROVARIANTES
# ------------------------------------------------------------------------------

message("   - Generando Bubble Chart...")

microvariantes_plot <- microvariantes %>%
  mutate(
    country = factor(country, levels = sort(unique(country))),
    microvariant = factor(microvariant, levels = sort(unique(microvariant)))
  )

p_bubble <- ggplot(microvariantes_plot, aes(x = country, y = microvariant)) +
  geom_point(aes(color = continent), size = 5, alpha = 0.9) +
  geom_text(aes(label = cases), color = "white", size = 2, fontface = "bold") +
  scale_color_manual(values = paleta_continentes) +
  labs(title = "Distribución Mundial de Microvariantes Y-STR",
       x = "País", y = "Microvariante") +
  theme_bw(base_size = 12)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # <--- Aquí está la magia
  )

ggsave(
  here::here("outputs", "figuras", "Obj1.2_BubbleChart.png"),
  p_bubble, width = 30, height = 12, units = "cm", dpi = 300
)

# ------------------------------------------------------------------------------
# 6. MAPA FILOGEOGRÁFICO
# ------------------------------------------------------------------------------

message("   - Generando mapa filogeográfico...")

world_map <- map_data("world")

microvariantes_mapa <- microvariantes %>%
  mutate(region_map = countrycode(country, origin = "country.name", destination = "country.name"))

p_mapa <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "grey92", color = "white", linewidth = 0.1) +
  geom_map(data = microvariantes_mapa, map = world_map,
           aes(fill = cases, map_id = region_map),
           color = "white", linewidth = 0.1) +
  scale_fill_gradient(low = "#fde0dd", high = "#cb181d", name = "Casos") +
  facet_wrap(~ microvariant, ncol = 3) +
  labs(title = "Distribución Filogeográfica de Microvariantes Y-STR") +
  theme_void()

ggsave(
  here::here("outputs", "figuras", "Obj1.3_Mapa_Filogeografico.png"),
  p_mapa, width = 28, height = 18, units = "cm", dpi = 300
)

message(">> Objetivo 1 completado.")
