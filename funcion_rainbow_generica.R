# funcion_rainbow_generica.R
# Función genérica para crear plots de comparación rainbow entre tratamientos
# 
# === DESCRIPCIÓN ===
# Esta función crea plots de comparación entre múltiples tratamientos usando
# un esquema de colores rainbow para ordenar los taxones por abundancia.
# 
# === PARÁMETROS ===
# - abundance_data: tabla de abundancias (taxones como filas, muestras como columnas)
# - metadata: tabla de metadatos con información de tratamientos
# - treatment_col: nombre de la columna en metadata que indica el tratamiento
# - sample_col: nombre de la columna en metadata que indica el ID de muestra (default: "Sample")
# NOTA: Las filas del abundance_data deben corresponder a taxones, las columnas a muestras
# - center_treatment: tratamiento que va en el centro del plot (para ordenamiento)
# - taxon_level: texto para etiquetas y títulos (opcional, default: NULL)
# - output_dir: directorio de salida (opcional)
# - show_labels: mostrar etiquetas en el plot (default: FALSE)
# - label_top_abundant: número de taxones más abundantes a etiquetar (default: 15)
# - label_invaders: número de invasores a etiquetar (default: 10)

library(tidyverse)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(scales)


# --- FUNCIÓN PRINCIPAL DE COMPARACIÓN RAINBOW ---------------------------------
rainbow_comparison_plot <- function(
  abundance_data,
  metadata,
  treatment_col,
  sample_col = "Sample",
  center_treatment = NULL,
  highlight_invaders = FALSE,
  taxon_level = NULL,
  output_dir = NULL,
  show_labels = FALSE,
  label_top_abundant = 15,
  label_invaders = 10
) {
  
  # Validar parámetros
  if (!treatment_col %in% colnames(metadata)) {
    stop("La columna de tratamiento '", treatment_col, "' no existe en los metadatos")
  }
  
  if (!sample_col %in% colnames(metadata)) {
    stop("La columna de muestra '", sample_col, "' no existe en los metadatos")
  }
  
  # Obtener tratamientos únicos
  treatments <- unique(metadata[[treatment_col]])
  treatments <- treatments[!is.na(treatments)]
  
  if (length(treatments) < 2) {
    stop("Se necesitan al menos 2 tratamientos para la comparación")
  }
  
  # Si se especifica center_treatment, reordenar tratamientos
  if (!is.null(center_treatment) && center_treatment %in% treatments) {
    other_treatments <- setdiff(treatments, center_treatment)
    treatments <- c(other_treatments[1], center_treatment, other_treatments[-1])
  }
  
  cat("Tratamientos a comparar:", paste(treatments, collapse = ", "), "\n")
  cat("Orden de tratamientos:", paste(treatments, collapse = " -> "), "\n")
  
  # Preparar datos de abundancia
  # Las filas deben ser taxones, las columnas muestras
  
  # Obtener columnas de muestras para cada tratamiento
  treatment_samples <- list()
  treatment_cols <- list()
  
  for (treatment in treatments) {
    treatment_metadata <- metadata %>%
      filter(!!sym(treatment_col) == treatment)
    
    treatment_samples[[treatment]] <- treatment_metadata[[sample_col]]
    treatment_cols[[treatment]] <- paste0("sn", treatment_samples[[treatment]])
    
    cat("  ", treatment, ":", length(treatment_samples[[treatment]]), "muestras\n")
  }
  
  # Verificar que haya muestras suficientes
  if (any(sapply(treatment_samples, length) == 0)) {
    stop("Algunos tratamientos no tienen muestras")
  }
  
  # Calcular abundancias medias para cada tratamiento
  # Primero asegurar que todas las columnas de muestras sean numéricas
  all_sample_cols <- unlist(treatment_cols)
  abundance_means <- abundance_data %>%
    mutate(
      across(all_of(all_sample_cols[all_sample_cols %in% colnames(abundance_data)]), 
             ~as.numeric(replace_na(., 0)))
    )
  
  # Agregar columnas de medias por tratamiento
  for (treatment in treatments) {
    cols_to_use <- treatment_cols[[treatment]]
    cols_to_use <- cols_to_use[cols_to_use %in% colnames(abundance_means)]
    
    # Filtrar solo columnas numéricas
    numeric_cols <- cols_to_use[sapply(abundance_means[cols_to_use], is.numeric)]
    
    cat("  -> Tratamiento:", treatment, "| Columnas encontradas:", length(cols_to_use), 
        "| Columnas numéricas:", length(numeric_cols), "\n")
    
    if (length(numeric_cols) > 0) {
      abundance_means[[paste0(treatment, "_mean")]] <- rowMeans(
        select(abundance_means, all_of(numeric_cols)), 
        na.rm = TRUE
      )
    } else {
      cat("  -> Warning: No hay columnas numéricas para", treatment, ", asignando 0\n")
      abundance_means[[paste0(treatment, "_mean")]] <- 0
    }
  }
  
  # Filtrar datos válidos
  mean_cols <- paste0(treatments, "_mean")
  abundance_means <- abundance_means %>%
    filter(
      if_all(all_of(mean_cols), ~!is.na(.) & !is.infinite(.))
    )
  
  # Clasificar taxones según la lógica de invasión
  # El tratamiento central (center_treatment) es el baseline para determinar invasores
  center_treatment_col <- paste0(center_treatment, "_mean")
  other_treatments_cols <- mean_cols[mean_cols != center_treatment_col]
  
  # Clasificar taxones según la lógica de invasión
  abundance_means <- abundance_means %>%
    mutate(
      is_invader = !!sym(center_treatment_col) == 0 & 
                   rowSums(select(., all_of(other_treatments_cols)) > 0) > 0,
      is_resident = !!sym(center_treatment_col) > 0,
      center_treatment_rank = rank(-!!sym(center_treatment_col), ties.method = "min"),
      total_mean = rowSums(select(., all_of(mean_cols))),
      invader_rank = rank(-total_mean, ties.method = "min")
    )
  
  # Preparar datos para el plot (formato largo)
  plot_data_long <- abundance_means %>%
    pivot_longer(
      cols = all_of(mean_cols),
      names_to = "condition",
      values_to = "abundance"
    ) %>%
    mutate(
      condition = str_replace(condition, "_mean", ""),
      condition = factor(condition, levels = treatments),
      condition_numeric = as.numeric(condition),
      jitter_offset = (as.numeric(factor(clade_name)) %% 7 - 3) * 0.02,
      condition_jitter = condition_numeric + jitter_offset
    )
  
  # Asignar grupos de color
  plot_data_long <- plot_data_long %>%
    mutate(
      color_group = if(highlight_invaders) {
        case_when(
          is_invader ~ paste0("TotalRank_", invader_rank),
          is_resident ~ "Resident",
          TRUE ~ "Other"
        )
      } else {
        case_when(
          is_invader ~ "Invader",
          TRUE ~ paste0("Rank_", center_treatment_rank)
        )
      },
      line_color = if(highlight_invaders) color_group else ifelse(is_invader, NA, color_group)
    )
  
  # Calcular tamaño de puntos
  plot_data_long <- plot_data_long %>%
    mutate(
      point_size = case_when(
        abundance == 0 ~ 0.5,
        abundance < 0.01 ~ 1,
        abundance < 0.1 ~ 2,
        TRUE ~ 3
      )
    )
  
  # Identificar taxones para etiquetar
  if (show_labels) {
    if (highlight_invaders) {
      top_abundant <- abundance_means %>%
        filter(is_resident) %>%
        mutate(max_abundant_abundance = pmax(!!!syms(mean_cols))) %>%
        arrange(desc(max_abundant_abundance)) %>%
        slice_head(n = label_top_abundant)
    } else {
      top_abundant <- abundance_means %>%
        filter(!is_invader & !!sym(center_treatment_col) > 0) %>%
        arrange(desc(!!sym(center_treatment_col))) %>%
        slice_head(n = label_top_abundant)
    }
    
    top_invaders <- abundance_means %>%
      filter(is_invader) %>%
      mutate(max_invader_abundance = rowSums(select(., all_of(other_treatments_cols)))) %>%
      arrange(desc(max_invader_abundance)) %>%
      slice_head(n = label_invaders)
    
    top_labels <- bind_rows(
      top_abundant %>% mutate(label_type = "abundant"),
      top_invaders %>% mutate(label_type = "invader")
    ) %>%
      distinct(clade_name, .keep_all = TRUE)
    
    top_labels$label <- top_labels$clade_name
  }
  
  # Configurar esquema de colores
  rank_groups <- unique(plot_data_long$color_group)
  
  # Ordenar rangos y configurar esquema de colores
  if (highlight_invaders) {
    # Highlight invaders: Rainbow para invasores, gris para residentes
    ranks_to_color <- rank_groups[grepl("TotalRank_", rank_groups)]
    ranks_to_color <- ranks_to_color[order(as.numeric(sub("TotalRank_", "", ranks_to_color)))]
    n_ranks <- length(ranks_to_color)
    
    color_scheme <- c(
      "Resident" = "gray50",
      setNames(rev(rainbow(n_ranks)), ranks_to_color)
    )
  } else {
    # Estándar: Gris para invasores, rainbow para abundancia
    ranks_numeric <- rank_groups[rank_groups != "Invader"] %>%
      sub("Rank_", "", .) %>%
      as.numeric()
    ranks_to_color <- rank_groups[rank_groups != "Invader"][order(ranks_numeric)]
    if ("Invader" %in% rank_groups) {
      ranks_to_color <- c(ranks_to_color, "Invader")
    }
    
    n_ranks <- length(ranks_to_color)
    color_scheme <- c(
      "Invader" = "gray50",
      setNames(rev(rainbow(n_ranks - 1)), ranks_to_color[ranks_to_color != "Invader"])
    )
  }
  
  # Crear el plot
  p <- ggplot(plot_data_long, aes(x = condition_jitter, y = abundance)) +
    # Borde blanco para highlight (solo si highlight_invaders = TRUE)
    {if(highlight_invaders) geom_point(aes(color = color_group, size = point_size + 0.5), color = "white", alpha = 1, shape = 16)} +
    # Líneas y puntos principales
    geom_line(aes(group = clade_name, color = line_color), alpha = 0.3, linewidth = 0.5) +
    geom_point(aes(color = color_group, size = point_size), alpha = 0.7) +
    scale_color_manual(
      values = color_scheme,
      na.value = NA, na.translate = FALSE, guide = "none"
    ) +
    scale_size_continuous(range = c(0.5, 4), guide = "none") +
    scale_x_continuous(
      breaks = 1:length(treatments),
      labels = treatments
    ) +
    scale_y_log10(labels = scales::comma_format()) +
    labs(
      title = paste0("Comparación Rainbow", 
                     if(highlight_invaders) " - Highlight Invaders" else "",
                     if(!is.null(taxon_level)) paste0(" - ", taxon_level) else "",
                     " - ", paste(treatments, collapse = " vs ")),
      subtitle = if(highlight_invaders) "Invasores coloreados por abundancia (Rainbow), Residentes en gris" else paste0("Orden: ", paste(treatments, collapse = " -> ")),
      x = "Tratamiento",
      y = "Abundancia Relativa (log10)",
      caption = paste0("n = ", nrow(abundance_means), " taxa",
                      if(highlight_invaders) "\nColores: Rainbow = invasores (ROJO=más abundante), Gris = residentes" else "\nColores: Rainbow = abundancia, Gris = invasores",
                      if(show_labels) paste0("\nEtiquetas: Top ", label_top_abundant, " abundantes + Top ", label_invaders, " invasores") else "\nSin etiquetas")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 9, hjust = 0),
      legend.position = "right"
    )
  
  # Agregar etiquetas si se solicitan
  if (show_labels && exists("top_labels") && nrow(top_labels) > 0) {
    label_data <- top_labels %>%
      select(clade_name, label, all_of(mean_cols)) %>%
      pivot_longer(
        cols = all_of(mean_cols),
        names_to = "condition",
        values_to = "abundance"
      ) %>%
      mutate(
        condition = str_replace(condition, "_mean", ""),
        condition = factor(condition, levels = treatments),
        condition_numeric = as.numeric(condition),
        jitter_offset = (as.numeric(factor(clade_name)) %% 7 - 3) * 0.02,
        condition_jitter = condition_numeric + jitter_offset
      ) %>%
      filter(abundance > 0)
    
    if (nrow(label_data) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = label_data,
        aes(x = condition_jitter, y = abundance, label = label),
        size = 2.5, color = "black", alpha = 0.8,
        max.overlaps = 20, force = 2, force_pull = 0.5
      )
    }
  }
  
  # Guardar plot si se especifica directorio
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    filename <- paste0("rainbow_comparison", 
                       if(highlight_invaders) "_highlight_invaders" else "",
                       "_", paste(treatments, collapse = "_"),
                       if(!is.null(taxon_level)) paste0("_", taxon_level) else "",
                       ".png")
    
    output_file <- file.path(output_dir, filename)
    ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
    cat("Plot guardado:", basename(output_file), "\n")
  }
  
  return(p)
}


# --- EJEMPLOS DE USO SIMPLES ------------------------------------------------

# Ejemplo 1: Comparación básica Exp2 (estándar)
metadata <- read_tsv("metadata_FULL.tsv", show_col_types = FALSE)
abundance_data <- read_tsv("abundance_table_wide_SGB.tsv", show_col_types = FALSE)
metadata_exp2 <- metadata %>% filter(Experiment == 2)

p1 <- rainbow_comparison_plot(
  abundance_data = abundance_data,
  metadata = metadata_exp2,
  treatment_col = "Treatment",
  center_treatment = "Treatment9",
  output_dir = "rainbow_examples"
)

# Ejemplo 2: Highlight invaders
p2 <- rainbow_comparison_plot(
  abundance_data = abundance_data,
  metadata = metadata_exp2,
  treatment_col = "Treatment",
  center_treatment = "Treatment9",
  highlight_invaders = TRUE,
  output_dir = "rainbow_examples",
  show_labels = TRUE
)

metadata_exp3 <- metadata %>% filter(Experiment == 3)

# Ejemplo 3: Sample y sample 17
p3 <- rainbow_comparison_plot(
  abundance_data = abundance_data,
  metadata = metadata_exp3,
  treatment_col = "Sample",
  center_treatment = "42",
  highlight_invaders = FALSE,
  output_dir = "rainbow_examples",
  show_labels = F
)
