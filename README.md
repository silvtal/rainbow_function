# Rainbow Function - Análisis de Microbioma

## Descripción

Este repositorio contiene una función genérica en R para crear visualizaciones de comparación "rainbow" entre diferentes tratamientos en estudios de microbioma. La función permite analizar la abundancia relativa de taxones entre múltiples condiciones experimentales, destacando patrones de invasión microbiana y abundancia.

## Características Principales

- **Visualización Rainbow**: Esquema de colores que ordena taxones por abundancia
- **Detección de Invasores**: Identifica taxones que aparecen en nuevos tratamientos pero no en el baseline
- **Comparación Multi-tratamiento**: Permite comparar múltiples condiciones simultáneamente
    - Se hace la media de abundancias de cada taxón en las muestras de cada condición.
    - Se puede modificar `treatment_col` y poner la columna de muestras. En ese caso, hay una columna por muestra. En `central_treatment` se indicaría la muestra a partir de la cual definir los invasores (muestra invadida).
- **Etiquetado Inteligente**: Opción para etiquetar los taxones más abundantes e invasores
- **Flexibilidad**: Parámetros configurables para diferentes tipos de análisis

## Estructura del Proyecto

```
rainbow_function/
├── funcion_rainbow_generica.R    # Función principal
├── abundance_table_wide_SGB.tsv  # Datos de abundancia (taxones x muestras)
├── metadata_FULL.tsv             # Metadatos de muestras
├── rainbow_examples/             # Ejemplos de plots generados
└── README.md                     # Este archivo
```

## Uso de la Función

### Parámetros Principales

```r
rainbow_comparison_plot(
  abundance_data,      # Tabla de abundancias
  metadata,            # Metadatos
  treatment_col,       # Columna de tratamientos
  sample_col = "Sample", # Columna de muestras
  center_treatment,    # Tratamiento central (baseline)
  highlight_invaders = FALSE, # Destacar invasores
  show_labels = FALSE, # Mostrar etiquetas
  output_dir = NULL    # Directorio de salida
)
```

### Parámetros Avanzados

- `label_top_abundant`: Número de taxones más abundantes a etiquetar (default: 15)
- `label_invaders`: Número de invasores a etiquetar (default: 10)
- `taxon_level`: Texto para etiquetas y títulos (opcional)


## Datos de ejemplo

### Abundance Table (`abundance_table_wide_SGB.tsv`)
- **Filas**: Taxones (SGBs) identificados por `clade_name`
- **Columnas**: Muestras identificadas por `sn[ID]`
- **Valores**: Abundancia relativa de cada taxón en cada muestra

### Metadata (`metadata_FULL.tsv`)
- **Sample**: ID de la muestra
- **Treatment**: Tipo de tratamiento
- **Group**: Grupo experimental
- **Day**: Día del experimento
- **Experiment**: Número de experimento
- **Cage, Line, Sex**: Variables adicionales

**Nota**: Los datos han sido modificados para preservar la privacidad:
- Números de muestra shuffleados aleatoriamente
- Nombres de tratamientos anonimizados
- Ruido añadido a las abundancias (factor: 0.15)

### Ejemplos de Uso

#### 1. Comparación Básica
```r
# Cargar datos
metadata <- read_tsv("metadata_FULL.tsv", show_col_types = FALSE)
abundance_data <- read_tsv("abundance_table_wide_SGB.tsv", show_col_types = FALSE)

# Filtrar experimento 2
metadata_exp2 <- metadata %>% filter(Experiment == 2)

# Crear plot básico
p1 <- rainbow_comparison_plot(
  abundance_data = abundance_data,
  metadata = metadata_exp2,
  treatment_col = "Treatment",
  center_treatment = "Treatment9",
  output_dir = "rainbow_examples"
)
```

#### 2. Destacar Invasores
```r
# Plot con invasores destacados
p2 <- rainbow_comparison_plot(
  abundance_data = abundance_data,
  metadata = metadata_exp2,
  treatment_col = "Treatment",
  center_treatment = "Treatment9",
  highlight_invaders = TRUE,
  show_labels = TRUE,
  output_dir = "rainbow_examples"
)
```
