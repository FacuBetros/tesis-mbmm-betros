# ==============================================================================
# 05 - OBJETIVO 4: Evaluar la asociación entre los haplogrupos del cromosoma Y 
# y los componentes de ancestría estimados mediante marcadores de herencia biparental...
# Autor: Facundo Betros
# ==============================================================================

message(">> Ejecutando Objetivo 4...")

# ------------------------------------------------------------------------------
# 0. PREPARACIÓN DE DATOS
# ------------------------------------------------------------------------------

message("   - Preparando datos de ancestría...")

ancestria_sistemas_hg <- ancestria_sistemas %>%
  left_join(
    muestras %>% select(id, cat_y, grupo_y, mt_dna_nat),
    by = c("sample_id" = "id")
  ) %>%
  mutate(
    # DESAGRUPADO: Tomamos el haplogrupo tal cual viene, sin agrupar en "Otros"
    hg = str_to_upper(trimws(cat_y))
  ) %>%
  # Limpiamos posibles celdas vacías para no tener una barra de "NA" en los gráficos
  mutate(hg = if_else(hg %in% c("", "N/D", "ND", NA), NA_character_, hg))

# ------------------------------------------------------------------------------
# 1. ANCESTRÍA PROMEDIO POR SISTEMA
# ------------------------------------------------------------------------------

message("   - Calculando ancestría promedio por sistema...")

ancestria_promedio <- tibble(
  Sistema = rep(c("STR", "InDel", "AIMs", "Cromosoma X"), each = 3),
  Componente = rep(c("Africano", "Europeo", "Nativo"), times = 4),
  Valor = c(
    mean(ancestria_sistemas$str_af, na.rm = TRUE), mean(ancestria_sistemas$str_eu, na.rm = TRUE), mean(ancestria_sistemas$str_na, na.rm = TRUE),
    mean(ancestria_sistemas$ind_af, na.rm = TRUE), mean(ancestria_sistemas$ind_eu, na.rm = TRUE), mean(ancestria_sistemas$ind_na, na.rm = TRUE),
    mean(ancestria_sistemas$aims_af, na.rm = TRUE), mean(ancestria_sistemas$aims_eu, na.rm = TRUE), mean(ancestria_sistemas$aims_na, na.rm = TRUE),
    mean(ancestria_sistemas$cx_af, na.rm = TRUE), mean(ancestria_sistemas$cx_eu, na.rm = TRUE), mean(ancestria_sistemas$cx_na, na.rm = TRUE)
  )
) %>%
  mutate(Sistema = factor(Sistema, levels = c("STR", "InDel", "AIMs", "Cromosoma X")))

p_ancestria <- ggplot(ancestria_promedio, aes(x = Sistema, y = Valor, fill = Componente)) +
  geom_col(linewidth = 0.3) +
  geom_text(aes(label = ifelse(Valor >= 0.03, paste0(round(Valor * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3) +
  geom_text(aes(label = ifelse(Valor < 0.03, paste0(round(Valor * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            vjust = -1, # Mueve el texto hacia la derecha, fuera de la barra
            color = "black", fontface = "bold", size = 3) +
  scale_fill_manual(values = paleta_ancestria) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Ancestría Promedio por Sistema de Marcadores",
       x = NULL, y = "Proporción", fill = "Componente") +
  theme_bw(base_size = 12)

ggsave(here::here("outputs", "figuras", "Obj4.1_Ancestria_Promedio.png"),
       p_ancestria, width = 18, height = 12, units = "cm", dpi = 300)

# ------------------------------------------------------------------------------
# 2. FUNCIONES PARA GRÁFICOS STRUCTURE Y PROMEDIO POR HAPLOGRUPO
# ------------------------------------------------------------------------------

grafico_structure <- function(datos, af, eu, na, titulo, archivo){
  datos_sistema <- datos %>%
    filter(!is.na(.data[[af]]) & !is.na(.data[[eu]]) & !is.na(.data[[na]]))
  
  if (nrow(datos_sistema) == 0) return(NULL)
  
  datos_long <- datos_sistema %>%
    select(sample_id, any_of(c(af, eu, na))) %>%
    pivot_longer(cols = any_of(c(af, eu, na)),
                 names_to = "componente", values_to = "valor") %>%
    mutate(componente = case_when(
      componente == af ~ "Africano", componente == eu ~ "Europeo",
      componente == na ~ "Nativo", TRUE ~ componente
    ))
  
  orden <- datos_sistema %>% arrange(sample_id) %>% pull(sample_id)
  datos_long <- datos_long %>% mutate(sample_id = factor(sample_id, levels = unique(orden)))
  
  p <- ggplot(datos_long, aes(x = sample_id, y = valor, fill = componente)) +
    geom_col(width = 0.85) +
    scale_fill_manual(values = paleta_ancestria) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    labs(title = titulo, x = NULL, y = NULL, fill = "Componente") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
  
  n_muestras <- length(unique(datos_sistema$sample_id))
  ancho_proporcional <- max(12, n_muestras * 0.15 + 4)
  
  ggsave(here::here("outputs", "figuras", archivo),
         p, width = ancho_proporcional, height = 10, units = "cm", dpi = 300)
}

grafico_promedio_hg <- function(datos, af, eu, na, titulo, archivo){
  promedio <- datos %>%
    filter(!is.na(hg)) %>%
    filter(!is.na(.data[[af]]) & !is.na(.data[[eu]]) & !is.na(.data[[na]])) %>%
    group_by(hg) %>%
    summarise(
      Africano = mean(.data[[af]], na.rm = TRUE),
      Europeo  = mean(.data[[eu]], na.rm = TRUE),
      Nativo   = mean(.data[[na]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-hg, names_to = "Componente", values_to = "Valor")
  
  if (nrow(promedio) == 0) return(NULL)
  
  # Ordenar alfabéticamente todos los haplogrupos encontrados
  promedio$hg <- factor(promedio$hg, levels = sort(unique(promedio$hg)))
  
  p <- ggplot(promedio, aes(x = hg, y = Valor, fill = Componente)) +
    geom_col(color = "white") +
    geom_text(aes(label = ifelse(Valor >= 0.03, paste0(round(Valor * 100, 1), "%"), "")),
              position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 3) +
    geom_text(aes(label = ifelse(Valor < 0.03, paste0(round(Valor * 100, 1), "%"), "")),
              position = position_stack(vjust = 0.5),
              vjust = -1, # Mueve el texto hacia la derecha, fuera de la barra
              color = "black", fontface = "bold", size = 3) +
    scale_fill_manual(values = paleta_ancestria) +
    scale_y_continuous(labels = scales::percent_format(),
                       expand = expansion(mult = c(0, 0.08))) +
    coord_cartesian(clip = "off") +
    labs(title = titulo, x = "Haplogrupo Y", y = NULL, fill = "Componente") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), # Giramos el texto 45°
          legend.position = "bottom")
  
  # Ancho dinámico para acomodar todos los haplogrupos
  n_hg <- length(unique(promedio$hg))
  ancho_fig <- max(16, n_hg * 1.5 + 3)
  
  ggsave(here::here("outputs", "figuras", archivo),
         p, width = ancho_fig, height = 12, units = "cm", dpi = 300)
}

# ------------------------------------------------------------------------------
# 3. EJECUCIÓN DE GRÁFICOS STRUCTURE Y PROMEDIO
# ------------------------------------------------------------------------------

message("   - Generando gráficos STRUCTURE y promedios por haplogrupo...")

grafico_structure(ancestria_sistemas_hg, "str_af", "str_eu", "str_na", "Estructura Genética: STRs", "Obj4_1_Structure_STRs.png")
grafico_structure(ancestria_sistemas_hg, "ind_af", "ind_eu", "ind_na", "Estructura Genética: InDels", "Obj4_1_Structure_InDels.png")
grafico_structure(ancestria_sistemas_hg, "aims_af", "aims_eu", "aims_na", "Estructura Genética: AIMs", "Obj4_1_Structure_AIMs.png")
grafico_structure(ancestria_sistemas_hg, "cx_af", "cx_eu", "cx_na", "Estructura Genética: Cromosoma X", "Obj4_1_Structure_CromosomaX.png")

grafico_promedio_hg(ancestria_sistemas_hg, "str_af", "str_eu", "str_na", "Ancestría Promedio: STRs", "Obj4_2_Promedio_STRs.png")
grafico_promedio_hg(ancestria_sistemas_hg, "ind_af", "ind_eu", "ind_na", "Ancestría Promedio: InDels", "Obj4_2_Promedio_InDels.png")
grafico_promedio_hg(ancestria_sistemas_hg, "aims_af", "aims_eu", "aims_na", "Ancestría Promedio: AIMs", "Obj4_2_Promedio_AIMs.png")
grafico_promedio_hg(ancestria_sistemas_hg, "cx_af", "cx_eu", "cx_na", "Ancestría Promedio: Cromosoma X", "Obj4_2_Promedio_CromosomaX.png")


# ------------------------------------------------------------------------------
# 4. BoxPlot Ancestria Nativa x Sistema x Hg CrY
# ------------------------------------------------------------------------------
# 0. Inspección rápida: conteo de NAs por haplogrupo y sistema
na_resumen <- ancestria_sistemas_hg %>%
  group_by(cat_y) %>%
  summarise(
    n = n(),
    n_str_na = sum(is.na(str_na)),
    n_ind_na = sum(is.na(ind_na)),
    n_aims_na = sum(is.na(aims_na)),
    n_cx_na = sum(is.na(cx_na))
  ) %>%
  arrange(n_str_na + n_ind_na + n_aims_na + n_cx_na)

print(na_resumen)   # revisá esto: te muestra qué haplogrupos tienen NAs

# 1. haplos_validos - Opción B (al menos 1 individuo no-NA por sistema)
haplos_validos <- ancestria_sistemas_hg %>%
  group_by(cat_y) %>%
  summarise(
    ok = all(
      c(
        any(!is.na(str_na)),
        any(!is.na(ind_na)),
        any(!is.na(aims_na)),
        any(!is.na(cx_na))
      )
    )
  ) %>%
  filter(ok) %>%
  pull(cat_y)

# Mensaje de control
if(length(haplos_validos) == 0) {
  stop("No se encontraron haplogrupos válidos bajo la condición seleccionada. Revisá 'na_resumen' para ver por qué.")
} else {
  message("Haplogrupos válidos encontrados: ", paste(haplos_validos, collapse = ", "))
}

# 2. Renombrar columnas ANTES del pivot y crear df
df_box_multi <- ancestria_sistemas_hg %>%
  filter(cat_y %in% haplos_validos) %>%
  rename(
    "STRs Autosómicos"   = str_na,
    "InDels Autosómicos" = ind_na,
    "AIMs"               = aims_na,
    "Cromosoma X"        = cx_na
  ) %>%
  pivot_longer(
    cols = c("STRs Autosómicos", "InDels Autosómicos", "AIMs", "Cromosoma X"),
    names_to = "Sistema",
    values_to = "Nativo"
  ) %>%
  mutate(
    Sistema = factor(Sistema, levels = c("STRs Autosómicos", "InDels Autosómicos", "AIMs", "Cromosoma X")),
    cat_y = factor(cat_y)
  )

# 3. Verificaciones antes de graficar
message("Filas en df_box_multi: ", nrow(df_box_multi))
message("Niveles Sistema: ", paste(levels(df_box_multi$Sistema), collapse = " | "))
message("Niveles paleta (deben coincidir): ", paste(names(paleta_sistemas), collapse = " | "))

# 4. Si hay filas, graficar; si no, avisar y no intentar scale_fill_manual
if(nrow(df_box_multi) == 0) {
  stop("df_box_multi tiene 0 filas. No se puede graficar. Revisá 'na_resumen' y la lógica de 'haplos_validos'.")
} else {
  paleta_sistemas <- c(
    "STRs Autosómicos"   = "#377eb8",
    "InDels Autosómicos" = "#4daf4a",
    "AIMs"               = "#e41a1c",
    "Cromosoma X"        = "#f0ad4e"
  )
  
  p_box_multi <- ggplot(df_box_multi, aes(x = cat_y, y = Nativo, fill = Sistema)) +
    geom_boxplot(alpha = 0.85, outlier.alpha = 0.4, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = paleta_sistemas) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Distribución de Ancestría Nativa por Sistema y Linaje Paterno",
         x = "Haplogrupo del Cromosoma Y", y = "Proporción de Ancestría Nativa Americana", fill = "Sistema de Marcadores") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.position = "right")
  
  # Guardar
  ggsave(filename = here::here("outputs", "figuras", "Obj4_Boxplot_Multisistema.png"),
         plot = p_box_multi, width = 28, height = 16, units = "cm", dpi = 300)
}

# ------------------------------------------------------------------------------
# 4. INTERACCIÓN Y - mtDNA
# ------------------------------------------------------------------------------

message("   - Analizando interacción entre linajes Y y mtDNA...")

col_mt <- "mt_dna"

if (!is.na(col_mt)) {
  df_mt_plot <- ancestria_sistemas_hg %>%
    filter(!is.na(hg)) %>%
    mutate(mtDNA_clado = substr(.data[[col_mt]], 1, 2)) %>%
    filter(mtDNA_clado %in% c("A2", "B2", "C1", "D1")) %>%
    count(hg, mtDNA_clado)
  
  paleta_mt <- c("A2"="#e41a1c", "B2"="#f781bf", "C1"="#a65628", "D1"="#984ea3")
  
  p_mt_y <- ggplot(df_mt_plot, aes(x = hg, y = n, fill = mtDNA_clado)) +
    geom_col(position = "fill", width = 0.7) +
    geom_text(aes(label = n), position = position_fill(vjust = 0.5), color = "white", fontface = "bold", size = 3.5) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = paleta_mt, name = "Clado mtDNA") +
    labs(title = "Interacción de Linajes: Haplogrupo Y vs mtDNA Nativo", x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) # Texto en ángulo
  
  n_hg_mt <- length(unique(df_mt_plot$hg))
  ancho_mt <- max(18, n_hg_mt * 1.5 + 4)
  
  ggsave(here::here("outputs", "figuras", "Obj4_Interaccion_Y_mtDNA.png"),
         p_mt_y, width = ancho_mt, height = 12, units = "cm", dpi = 300)
}

# ------------------------------------------------------------------------------
# 5. ASOCIACIÓN ESTADÍSTICA: HAPLOGRUPOS Y vs. ANCESTRÍA GENÓMICA
# ------------------------------------------------------------------------------

message("   - Calculando pruebas estadísticas (Kruskal-Wallis, Spearman y Dunn)...")

# ------------------------------------------------------------------------------
# 5.1 Pruebas de Kruskal-Wallis por Sistema de Marcadores
# ------------------------------------------------------------------------------

# 1. STRs autosómicos (str_na)
datos_kw_str <- ancestria_sistemas_hg %>%
  filter(!is.na(cat_y), !is.na(str_na)) %>%
  mutate(cat_y = as.factor(cat_y))

kruskal_str <- kruskal.test(str_na ~ cat_y, data = datos_kw_str)

# 2. InDels autosómicos (ind_na)
datos_kw_ind <- ancestria_sistemas_hg %>%
  filter(!is.na(cat_y), !is.na(ind_na)) %>%
  mutate(cat_y = as.factor(cat_y))

kruskal_ind <- kruskal.test(ind_na ~ cat_y, data = datos_kw_ind)

# 3. AIMs-SNPs autosómicos (aims_na)
datos_kw_aims <- ancestria_sistemas_hg %>%
  filter(!is.na(cat_y), !is.na(aims_na)) %>%
  mutate(cat_y = as.factor(cat_y))

kruskal_aims <- kruskal.test(aims_na ~ cat_y, data = datos_kw_aims)

# 4. Cromosoma X (cx_na)
datos_kw_cx <- ancestria_sistemas_hg %>%
  filter(!is.na(cat_y), !is.na(cx_na)) %>%
  mutate(cat_y = as.factor(cat_y))

kruskal_cx <- kruskal.test(cx_na ~ cat_y, data = datos_kw_cx)

# Imprimir resúmenes en consola
print("--- Kruskal-Wallis: STRs ---"); print(kruskal_str)
print("--- Kruskal-Wallis: InDels ---"); print(kruskal_ind)
print("--- Kruskal-Wallis: AIMs ---"); print(kruskal_aims)
print("--- Kruskal-Wallis: Cromosoma X ---"); print(kruskal_cx)


# ------------------------------------------------------------------------------
# 5.2 Comparaciones Múltiples Post-Hoc (Dunn para STRs)
# ------------------------------------------------------------------------------
# Se realiza únicamente para STRs ya que fue el único sistema con p < 0.05 global
if (requireNamespace("FSA", quietly = TRUE)) {
  dunn_resultado <- FSA::dunnTest(str_na ~ cat_y, data = datos_kw_str, method = "bh")
  
  # Crear carpeta si no existe y exportar
  dir.create(here::here("outputs", "tablas"), recursive = TRUE, showWarnings = FALSE)
  write.csv(dunn_resultado$res, 
            here::here("outputs", "tablas", "Obj4_Dunn_PostHoc_HaplogruposY.csv"), 
            row.names = FALSE)
  message("     [+] Resultados de Post-Hoc Dunn exportados correctamente.")
} else {
  warning("El paquete 'FSA' no está instalado. No se pudo ejecutar el Post-Hoc de Dunn.")
}


# ------------------------------------------------------------------------------
# 5.3 Correlación de Spearman (Origen Paterno Q vs. Ancestría)
# ------------------------------------------------------------------------------

# Crear variable binaria de origen paterno (1 = Q / Nativo, 0 = Otros / Alóctono)
datos_corr <- ancestria_sistemas_hg %>%
  mutate(origen_paterno_num = if_else(cat_y == "Q", 1, 0)) %>%
  filter(!is.na(origen_paterno_num))

# 1. Correlación para STRs
cor_str <- NULL
if ("str_na" %in% names(datos_corr)) {
  datos_str <- datos_corr %>% filter(!is.na(str_na))
  cor_str <- cor.test(datos_str$origen_paterno_num, datos_str$str_na, method = "spearman", exact = FALSE)
}

# 2. Correlación para InDels
cor_ind <- NULL
if ("ind_na" %in% names(datos_corr)) {
  datos_ind <- datos_corr %>% filter(!is.na(ind_na))
  cor_ind <- cor.test(datos_ind$origen_paterno_num, datos_ind$ind_na, method = "spearman", exact = FALSE)
}

# 3. Correlación para AIMs-SNPs
cor_aims <- NULL
if ("aims_na" %in% names(datos_corr)) {
  datos_aims <- datos_corr %>% filter(!is.na(aims_na))
  cor_aims <- cor.test(datos_aims$origen_paterno_num, datos_aims$aims_na, method = "spearman", exact = FALSE)
}

# 4. Correlación para Cromosoma X
cor_cx <- NULL
if ("cx_na" %in% names(datos_corr)) {
  datos_cx <- datos_corr %>% filter(!is.na(cx_na))
  cor_cx <- cor.test(datos_cx$origen_paterno_num, datos_cx$cx_na, method = "spearman", exact = FALSE)
}

# Imprimir correlaciones en consola
print("--- Spearman: STRs ---"); print(cor_str)
print("--- Spearman: InDels ---"); print(cor_ind)
print("--- Spearman: AIMs ---"); print(cor_aims)
print("--- Spearman: Cromosoma X ---"); print(cor_cx)


# ------------------------------------------------------------------------------
# 5.4 Exportación de Resumen Estadístico Consolidado (CSV)
# ------------------------------------------------------------------------------

tabla_resultados <- tibble::tibble(
  Analisis = c(
    rep("Kruskal-Wallis (Ancestría Nativa ~ Haplogrupo Y)", 4),
    rep("Spearman (Origen Paterno Q vs Ancestría Nativa)", 4)
  ),
  Sistema = c(
    "STRs Autosómicos", "InDels Autosómicos", "AIMs-SNPs Autosómicos", "Cromosoma X",
    "STRs Autosómicos", "InDels Autosómicos", "AIMs-SNPs Autosómicos", "Cromosoma X"
  ),
  Estadistico_Nombre = c(
    rep("Chi-cuadrado (χ²)", 4),
    rep("Rho (rho)", 4)
  ),
  Estadistico_Valor = c(
    round(as.numeric(kruskal_str$statistic), 3),
    round(as.numeric(kruskal_ind$statistic), 3),
    round(as.numeric(kruskal_aims$statistic), 3),
    round(as.numeric(kruskal_cx$statistic), 3),
    if (!is.null(cor_str)) round(as.numeric(cor_str$estimate), 3) else NA,
    if (!is.null(cor_ind)) round(as.numeric(cor_ind$estimate), 3) else NA,
    if (!is.null(cor_aims)) round(as.numeric(cor_aims$estimate), 3) else NA,
    if (!is.null(cor_cx)) round(as.numeric(cor_cx$estimate), 3) else NA
  ),
  p_valor = c(
    kruskal_str$p.value,
    kruskal_ind$p.value,
    kruskal_aims$p.value,
    kruskal_cx$p.value,
    if (!is.null(cor_str)) cor_str$p.value else NA,
    if (!is.null(cor_ind)) cor_ind$p.value else NA,
    if (!is.null(cor_aims)) cor_aims$p.value else NA,
    if (!is.null(cor_cx)) cor_cx$p.value else NA
  )
)

# Exportar tabla completa a CSV
write.csv(tabla_resultados, 
          here::here("outputs", "tablas", "Obj4_Resultados_Estadisticos.csv"), 
          row.names = FALSE)

message("     [+] Tabla consolidada de resultados estadísticos exportada a CSV.")
message(">> Objetivo 4 completado.")