# tesis-mbmm-betros
Tesis – Genética Forense y Poblacional
Análisis de linajes paternos, ancestría genómica y poblaciones comparativas
Autor: Facundo Betros
Repositorio: Análisis completo reproducible en R

---
📌 Descripción del proyecto
Este repositorio contiene el análisis completo de mi tesis de grado sobre genética forense y poblacional, enfocado en:
• Haplogrupos del cromosoma Y
• Parámetros forenses (DC, MP, PD, HD)
• Microvariantes Y-STR
• Comparación con bases internacionales (YHRD)
• Ancestría autosómica, AIMs, InDels y cromosoma X
• Interacción entre linajes uniparentales (Y–mtDNA)
• Estadísticas de mestizaje
• Trayectorias genómicas individuales
Todo el análisis está implementado en R, con scripts modulares y exportación automática de figuras y tablas.

📁 Estructura del repositorio
├── datos/
│   └── Tesis.xlsx
│
├── scripts/
│   ├── 00_funciones.R
│   ├── 01_carga_datos.R
│   ├── 02_objetivo1.R
│   ├── 03_objetivo2.R
│   ├── 04_objetivo3.R
│   ├── 05_objetivo4.R
│   ├── 06_objetivo5.R
│   └── main.R
│
├── outputs/
│   ├── figuras/
│   └── tablas/
│
└── README.md

🚀 Cómo ejecutar el análisis completo
1. Clonar el repositorio:
git clone https://github.com/<tu_usuario>/<tu_repo>.git
2. Abrir el proyecto en RStudio
3. Instalar dependencias (se hace automáticamente desde los scripts)
4. Ejecutar: source("scripts/main.R")
Esto generará todas las figuras y tablas dentro de outputs/

📊 Resultados generados
Objetivo 1
• Parámetros forenses por kit
• Haplotipos únicos y repetidos
• Microvariantes Y-STR
• Bubble chart
• Mapa filogeográfico
• Concordancia NevGen vs Whit Athey
Objetivo 2
• Frecuencia de haplogrupos
• Diversidad genética
• Haplotipos únicos por haplogrupo
• Test de Fisher
• Índice Kappa
Objetivo 3
• Gráfico ternario (Q–R–E)
• MA-Plot (Argentina vs Mundo)
• Chi-cuadrado
Objetivo 4
• Structure plots (STR, InDel, AIMs, CX)
• Ancestría promedio por haplogrupo
• Interacción Y–mtDNA
Objetivo 5
• Correlación AU vs CX
• Wilcoxon pareado
• Delta de ancestría
• Gráfico aluvial
• Trayectorias paralelas

🧬 Dependencias principales
• tidyverse
• readxl
• janitor
• ggtern
• irr
• ape
• ggtree
• factoextra
• ggplot2
• here
Todas se cargan automáticamente desde los scripts.

Para consultas sobre el análisis, metodología o implementación:
Lic. Facundo Betros
Ciudad Autónoma de Buenos Aires, Argentina
GitHub: FacuBetros


