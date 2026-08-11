# ==============================================================================
# 05 - OBJETIVO 4: ANCESTRÍA GENÓMICA INTEGRAL
# Autor: Facundo Betros
# ==============================================================================

message(">> Ejecutando Objetivo 4...")

# ------------------------------------------------------------------------------
# 0. PREPARACIÓN DE DATOS
# ------------------------------------------------------------------------------

message("   - Preparando datos de ancestría...")

ancestria_sistemas_hg <- ancestria_sistemas %>%
  left_join(
    muestras %>% select(id, cat_y, grupo_y, mt_dna, mt_dna_nat),
    by = c("sample_id" = "id")
  ) %>%
  mutate(
    hg = case_when(
      cat_y == "E1B1A" ~ "E1B1A",
      cat_y == "E1B1B" ~ "E1B1B",
      cat_y == "I"     ~ "I",
      cat_y == "J"     ~ "J",
      cat_y == "Q"     ~ "Q",
      cat_y == "R"     ~ "R",
      TRUE ~ "Otros"
    )
  )

# ------------------------------------------------------------------------------
# 1. ANCESTRÍA PROMEDIO POR SISTEMA
# ------------------------------------------------------------------------------

message("   - Calculando ancestría promedio por sistema...")

ancestria_promedio <- tibble(
  Sistema = c("STR", "STR", "STR",
              "InDel", "InDel", "InDel",
              "AIMs", "AIMs", "AIMs",
              "Cromosoma X", "Cromosoma X", "Cromosoma X"),
  Componente = c("Africano", "Europeo", "Nativo",
                 "Africano", "Europeo", "Nativo",
                 "Africano", "Europeo", "Nativo",
                 "Africano", "Europeo", "Nativo"),
  Valor = c(
    mean(ancestria_sistemas_hg$str_af, na.rm = TRUE),
    mean(ancestria_sistemas_hg$str_eu, na.rm = TRUE),
    mean(ancestria_sistemas_hg$str_na, na.rm = TRUE),
    
    mean(ancestria_sistemas_hg$ind_af, na.rm = TRUE),
    mean(ancestria_sistemas_hg$ind_eu, na.rm = TRUE),
    mean(ancestria_sistemas_hg$ind_na, na.rm = TRUE),
    
    mean(ancestria_sistemas_hg$aims_af, na.rm = TRUE),
    mean(ancestria_sistemas_hg$aims_eu, na.rm = TRUE),
    mean(ancestria_sistemas_hg$aims_na, na.rm = TRUE),
    
    mean(ancestria_sistemas_hg$cx_af, na.rm = TRUE),
    mean(ancestria_sistemas_hg$cx_eu, na.rm = TRUE),
    mean(ancestria_sistemas_hg$cx_nativo, na.rm = TRUE)
  )
)

p_ancestria <- ggplot(ancestria_promedio, aes(x = Sistema, y = Valor, fill = Componente)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(Valor * 100, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = paleta_ancestria) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Ancestría Promedio por Sistema de Marcadores",
       x = NULL, y = "Proporción", fill = "Componente") +
  theme_bw(base_size = 14)

ggsave(
  here::here("outputs", "figuras", "Obj4_Ancestria_Promedio.png"),
  p_ancestria, width = 18, height = 12, units = "cm", dpi = 300
)

# ------------------------------------------------------------------------------
# 2. FUNCIONES PARA GRÁFICOS STRUCTURE Y PROMEDIO POR HAPLOGRUPO
# ------------------------------------------------------------------------------

grafico_structure <- function(datos, af, eu, na, titulo, archivo){
  datos_long <- datos %>%
    select(sample_id, all_of(c(af, eu, na))) %>%
    pivot_longer(cols = c(af, eu, na),
                 names_to = "componente",
                 values_to = "valor") %>%
    mutate(componente = recode(componente,
                               !!af := "Africano",
                               !!eu := "Europeo",
                               !!na := "Nativo"))
  
  orden <- datos %>% arrange(.data[[na]]) %>% pull(sample_id)
  datos_long <- datos_long %>% mutate(sample_id = factor(sample_id, levels = unique(orden)))
  
  p <- ggplot(datos_long, aes(x = sample_id, y = valor, fill = componente)) +
    geom_col(width = 1) +
    scale_fill_manual(values = paleta_ancestria) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = titulo, x = NULL, y = NULL, fill = "Componente") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "bottom")
  
  ggsave(here::here("outputs", "figuras", archivo),
         p, width = 32, height = 10, units = "cm", dpi = 300)
}

grafico_promedio_hg <- function(datos, af, eu, na, titulo, archivo){
  promedio <- datos %>%
    filter(hg %in% c("E1B1A", "E1B1B", "I", "J", "Q", "R")) %>%
    group_by(hg) %>%
    summarise(
      Africano = mean(.data[[af]], na.rm = TRUE),
      Europeo  = mean(.data[[eu]], na.rm = TRUE),
      Nativo   = mean(.data[[na]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-hg, names_to = "Componente", values_to = "Valor")
  
  promedio$hg <- factor(promedio$hg, levels = c("E1B1A", "E1B1B", "I", "J", "Q", "R"))
  
  p <- ggplot(promedio, aes(x = hg, y = Valor, fill = Componente)) +
    geom_col(color = "white") +
    geom_text(aes(label = paste0(round(Valor * 100), "%")),
              position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 5) +
    scale_fill_manual(values = paleta_ancestria) +
    scale_y_continuous(labels = scales::percent_format(),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = titulo, x = "Haplogrupo Y", y = NULL, fill = "Componente") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "bottom")
  
  ggsave(here::here("outputs", "figuras", archivo),
         p, width = 18, height = 12, units = "cm", dpi = 300)
}

# ------------------------------------------------------------------------------
# 3. EJECUCIÓN DE GRÁFICOS STRUCTURE Y PROMEDIO
# ------------------------------------------------------------------------------

message("   - Generando gráficos STRUCTURE y promedios por haplogrupo...")

grafico_structure(ancestria_sistemas_hg, "str_af", "str_eu", "str_na",
                  "Estructura Genética: STRs", "Obj4_1_Structure_STRs.png")

grafico_structure(ancestria_sistemas_hg, "ind_af", "ind_eu", "ind_na",
                  "Estructura Genética: InDels", "Obj4_1_Structure_InDels.png")

grafico_structure(ancestria_sistemas_hg, "aims_af", "aims_eu", "aims_na",
                  "Estructura Genética: AIMs", "Obj4_1_Structure_AIMs.png")

grafico_structure(ancestria_sistemas_hg, "cx_af", "cx_eu", "cx_nativo",
                  "Estructura Genética: Cromosoma X", "Obj4_1_Structure_CromosomaX.png")

grafico_promedio_hg(ancestria_sistemas_hg, "str_af", "str_eu", "str_na",
                    "Ancestría Promedio: STRs", "Obj4_2_Promedio_STRs.png")

grafico_promedio_hg(ancestria_sistemas_hg, "ind_af", "ind_eu", "ind_na",
                    "Ancestría Promedio: InDels", "Obj4_2_Promedio_InDels.png")

grafico_promedio_hg(ancestria_sistemas_hg, "aims_af", "aims_eu", "aims_na",
                    "Ancestría Promedio: AIMs", "Obj4_2_Promedio_AIMs.png")

grafico_promedio_hg(ancestria_sistemas_hg, "cx_af", "cx_eu", "cx_nativo",
                    "Ancestría Promedio: Cromosoma X", "Obj4_2_Promedio_CromosomaX.png")

# ------------------------------------------------------------------------------
# 4. INTERACCIÓN Y - mtDNA
# ------------------------------------------------------------------------------

message("   - Analizando interacción entre linajes Y y mtDNA...")

df_mt_plot <- ancestria_sistemas_hg %>%
  mutate(
    mtDNA_clado = substr(mt_dna, 1, 2)
  ) %>%
  filter(mtDNA_clado %in% c("A2", "B2", "C1", "D1")) %>%
  count(hg, mtDNA_clado)

paleta_mt <- c("A2"="#e41a1c", "B2"="#f781bf", "C1"="#a65628", "D1"="#984ea3")

p_mt_y <- ggplot(df_mt_plot, aes(x = hg, y = n, fill = mtDNA_clado)) +
  geom_col(position = "fill", color = "black", width = 0.7) +
  geom_text(aes(label = n),
            position = position_fill(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.5) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = paleta_mt, name = "Clado mtDNA") +
  labs(title = "Interacción de Linajes: Haplogrupo Y vs mtDNA Nativo",
       x = NULL, y = NULL) +
  theme_bw(base_size = 14)

ggsave(
  here::here("outputs", "figuras", "Obj4_Interaccion_Y_mtDNA.png"),
  p_mt_y, width = 18, height = 12, units = "cm", dpi = 300
)

message(">> Objetivo 4 completado.")
