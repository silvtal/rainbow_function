# funcion_rainbow_generica.R
# Función genérica para crear plots de comparación rainbow/taxonomía entre tratamientos
# 
# === DESCRIPCIÓN ===
# Esta función crea plots de comparación entre múltiples tratamientos usando
# un esquema de colores rainbow (por abundancia) o un esquema de color taxonómico
# (por reino, clase, etc.), respetando siempre el filtro de taxones "grises" (invasores o residentes).
# 
# === PARÁMETROS ===
# - abundance_data: tabla de abundancias (taxones como filas, muestras como columnas)
# - metadata: tabla de metadatos con información de tratamientos
# - treatment_col: nombre de la columna en metadata que indica el tratamiento
# - sample_col: nombre de la columna en metadata que indica el ID de muestra (default: "Sample")
# - center_treatment: tratamiento que va en el centro del plot (para ordenamiento)
# - highlight_invaders: booleano para resaltar invasores (solo aplica con coloring rainbow)
# - taxonomy_table: Tabla de taxonomía con una primera columna indicando el taxón (por ej la secuencia, "otu1"...) y las columnas taxonómicas
# - color_by_taxon_col: Nombre de la columna taxonómica a usar para colorear (ej: "class")
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
# library(scales) # Necesario para hue_pal()

# --- FUNCIÓN PRINCIPAL DE COMPARACIÓN RAINBOW/TAXONOMÍA -----------------------
rainbow_comparison_plot <- function(
    abundance_data,
    metadata,
    treatment_col,
    sample_col = "Sample",
    center_treatment = NULL,
    highlight_invaders = FALSE,
    taxonomy_table = NULL,
    color_by_taxon_col = NULL,
    taxon_level = NULL,
    output_dir = NULL,
    show_labels = FALSE,
    label_top_abundant = 15,
    label_invaders = 10
){
  
  # Determinar si se usa coloreado taxonómico
  use_taxonomy_coloring <- !is.null(taxonomy_table) && 
    !is.null(color_by_taxon_col) && 
    color_by_taxon_col %in% colnames(taxonomy_table)
  
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
  # Asegurar que los nombres de las filas (taxones) sean una columna llamada 'clade_name'
  if (!"clade_name" %in% colnames(abundance_data)) {
    # Suponemos que si no está 'clade_name', el taxón ID es el nombre de la fila
    abundance_data <- abundance_data %>%
      rownames_to_column("clade_name")
  }
  
  # Obtener columnas de muestras para cada tratamiento
  treatment_samples <- list()
  
  for (treatment in treatments) {
    treatment_metadata <- metadata %>%
      filter(!!sym(treatment_col) == treatment)
    
    treatment_samples[[treatment]] <- treatment_metadata[[sample_col]]
    
    cat("  ", treatment, ":", length(treatment_samples[[treatment]]), "muestras\n")
  }
  
  # Verificar que haya muestras suficientes
  if (any(sapply(treatment_samples, length) == 0)) {
    stop("Algunos tratamientos no tienen muestras")
  }
  
  # Calcular abundancias medias para cada tratamiento
  all_sample_cols <- unlist(treatment_cols)
  abundance_means <- abundance_data %>%
    mutate(
      across(all_of(all_sample_cols[all_sample_cols %in% colnames(abundance_data)]), 
             ~as.numeric(replace_na(., 0)))
    )
  
  # Agregar columnas de medias por tratamiento
  for (treatment in treatments) {
    cols_to_use <- treatment_samples[[treatment]]
    cols_to_use <- cols_to_use[cols_to_use %in% colnames(abundance_means)]
    
    # Filtrar solo columnas numéricas
    numeric_cols <- cols_to_use[sapply(abundance_means[cols_to_use], is.numeric)]
    
    cat("  -> Tratamiento:", treatment, "| Columnas encontradas:", length(cols_to_use), 
        "| Columnas numéricas:", length(numeric_cols), "\n")
    
    if (length(numeric_cols) > 0) {
      abundance_means[[paste0(treatment)]] <- rowMeans(
        select(abundance_means, all_of(numeric_cols)), 
        na.rm = TRUE
      )
    } else {
      cat("  -> Warning: No hay columnas numéricas para", treatment, ", asignando 0\n")
      abundance_means[[paste0(treatment)]] <- 0
    }
  }
  
  # Filtrar datos válidos
  abundance_means <- abundance_means %>%
    filter(
      if_all(all_of(treatments), ~!is.na(.) & !is.infinite(.))
    )
  
  # === CLASIFICACIÓN (SIEMPRE NECESARIA) ===========================
  other_treatments <- treatments[treatments != center_treatment_col]
  
  abundance_means <- abundance_means %>%
    mutate(
      is_invader = !!sym(center_treatment) == 0 & 
        rowSums(select(., all_of(other_treatments)) > 0) > 0,
      is_resident = !!sym(center_treatment) > 0,
      center_treatment_rank = rank(-!!sym(center_treatment), ties.method = "min"),
      total_mean = rowSums(select(., all_of(treatments))),
      invader_rank = rank(-total_mean, ties.method = "min")
    )
  
  # === Eliminar taxones con Abundancia Total 0 (Informativamente inútiles) ===
  abundance_means <- abundance_means %>%
    filter(total_mean > 0)
  # ======================================================================================
  
  # --- CONFIGURACIÓN DE COLOREADO (TAXONOMÍA FILTRADA O RAINBOW) -----------------------
  if (use_taxonomy_coloring) {
    cat("Usando esquema: Taxonomía (Coloreando solo taxones que no están en el grupo gris).\n")
    color_type <- "Taxonomy_Filtered"
    
    # Renombrar la primera columna de la tabla de taxonomía (Taxón ID) a 'clade_name'
    names(taxonomy_table)[1] <- "clade_name"
    
    # Unir datos de abundancia con la taxonomía
    abundance_means <- abundance_means %>%
      left_join(
        select(taxonomy_table, all_of(c("clade_name", color_by_taxon_col))),
        by = "clade_name"
      )
    
  } else {
    cat("Usando coloreado por abundancia (Rainbow/Ranking)\n")
    color_type <- "Rainbow"
  }
  
  # Preparar datos para el plot (formato largo)
  plot_data_long <- abundance_means %>%
    pivot_longer(
      cols = all_of(treatments),
      names_to = "condition",
      values_to = "abundance"
    ) %>%
    mutate(
      condition = factor(condition, levels = treatments),
      condition_numeric = as.numeric(condition),
      # Usar clade_name para el jitter, asumiendo que es la ID única del taxón
      jitter_offset = (as.numeric(factor(clade_name)) %% 7 - 3) * 0.02,
      condition_jitter = condition_numeric + jitter_offset
    )
  
  # === ASIGNAR GRUPOS DE COLOR Y LINEA =======================================
  
  # Bloque de IF/ELSE base para asignar is_gray_taxon (soluciona el error de tamaño)
  if (highlight_invaders) {
    # Si highlight_invaders es TRUE: Los Residentes son Grises.
    plot_data_long$is_gray_taxon <- plot_data_long$is_resident
  } else {
    # Si highlight_invaders es FALSE: Los Invasores son Grises.
    plot_data_long$is_gray_taxon <- plot_data_long$is_invader
  }
  
  if (use_taxonomy_coloring) {
    # TAXONOMÍA FILTRADA: Usar taxón para colorear lo que NO es gris, y grupo gris para lo que SÍ lo es.
    
    plot_data_long <- plot_data_long %>%
      mutate(
        color_group_base = factor(!!sym(color_by_taxon_col)), # El color taxonómico
        
        # Definir qué grupo es el gris y qué taxones se mapean a él
        gray_group_name = if_else(highlight_invaders, "Resident", "Invader"),
        
        # color_group: Si debe ser gris, usa el nombre del grupo gris. Si no, usa el taxón.
        color_group = if_else(
          is_gray_taxon, 
          gray_group_name, 
          as.character(color_group_base)
        ),
        
        # line_color: Las líneas grises (Invader en modo estándar) NO se deben dibujar (NA).
        # Las líneas grises (Resident en highlight mode) SÍ se deben dibujar (en gris/color_group).
        line_color = if_else(
          !highlight_invaders & is_gray_taxon, # Caso: Standard mode + Invader (gris/sin línea)
          NA_character_, 
          color_group # Caso: Coloreado por Taxonomía o Residentes (en highlight mode)
        )
      )
    
  } else {
    # RAINBOW ESTÁNDAR
    
    # Bloque IF/ELSE base para calcular color_group_calc (soluciona el error de tamaño en modo Rainbow)
    if (highlight_invaders) {
      # Highlight Invaders: Colorear Invasores, Residentes en Gris
      color_group_calc <- with(plot_data_long, 
                               if_else(is_invader, paste0("TotalRank_", invader_rank), "Resident"))
      line_color_calc <- color_group_calc # Las líneas deben ser coloreadas
    } else {
      # Estándar: Colorear Residentes por Rank, Invasores en Gris
      color_group_calc <- with(plot_data_long, 
                               if_else(is_invader, "Invader", paste0("Rank_", center_treatment_rank)))
      line_color_calc <- with(plot_data_long, 
                              if_else(is_invader, NA_character_, color_group_calc)) # Las líneas de invasores son NA
    }
    
    plot_data_long <- plot_data_long %>%
      mutate(
        color_group = color_group_calc,
        line_color = line_color_calc
      )
  }
  
  plot_data_long <- plot_data_long %>%
    mutate(color_group = factor(color_group))
  
  # Calcular tamaño de puntos
  plot_data_long <- plot_data_long %>%
    mutate(
      point_size = case_when(
        abundance == 0 ~ 0.5,
        abundance < 0.001 ~ 0.8,
        abundance < 0.01 ~ 1.5,
        abundance < 0.1 ~ 2.5,
        TRUE ~ 3.5
      )
    )
  
  # Identificar taxones para etiquetar (se usa total_mean, independiente del coloring)
  if (show_labels) {
    top_labels <- abundance_means %>%
      arrange(desc(total_mean)) %>%
      slice_head(n = label_top_abundant + label_invaders)
    
    top_labels$label <- top_labels$clade_name
  }
  
  # === CONFIGURAR ESQUEMA DE COLORES =========================================
  
  # Obtener todos los grupos de color (taxonómicos y gris)
  rank_groups <- unique(plot_data_long$color_group)
  
  if (use_taxonomy_coloring) {
    # TAXONOMÍA FILTRADA: Colores taxonómicos + Gris
    
    # Identificar el grupo gris
    gray_group_name <- if (highlight_invaders) "Resident" else "Invader"
    
    # Filtrar solo los grupos NO grises para definir la paleta
    tax_groups <- rank_groups[rank_groups != gray_group_name]
    n_groups <- length(tax_groups)
    
    # === Usar scales::hue_pal() para colores perceptualmente uniformes ===
    # tax_colors <- scales::hue_pal()(n_groups) # ESCALA ALTERNATIVA
    tax_colors <- rainbow(n_groups)
    # ==================================================================================
    
    color_scheme <- c(
      setNames("gray50", gray_group_name), # El grupo gris primero
      setNames(tax_colors, tax_groups)
    )
    
    # Configuración de leyenda
    color_guide <- guide_legend(title = color_by_taxon_col, ncol = 1)
    
    plot_title <- paste0("Comparación: Taxonomía (Filtrada por Invasión)")
    plot_subtitle <- paste0("Coloreado por ", color_by_taxon_col, ". ", 
                            if(!is.null(taxon_level)) taxon_level else "")
    plot_caption <- paste0("Colores: Taxonomía (", color_by_taxon_col, ") vs ", gray_group_name, " (Gris)")
    
  } else {
    # Color por Abundancia (RAINBOW)
    
    if (highlight_invaders) {
      # Highlight invaders: Rainbow para invasores, gris para residentes
      ranks_to_color <- rank_groups[grepl("TotalRank_", rank_groups)]
      ranks_to_color <- ranks_to_color[order(as.numeric(sub("TotalRank_", "", ranks_to_color)))]
      n_ranks <- length(ranks_to_color)
      
      color_scheme <- c(
        "Resident" = "gray50",
        setNames(rev(rainbow(n_ranks)), ranks_to_color)
      )
      
      plot_title <- "Comparación Rainbow - Highlight Invaders"
      plot_subtitle <- "Invasores coloreados por abundancia (Rainbow), Residentes en gris"
      plot_caption <- paste0("Colores: Rainbow = invasores (ROJO=más abundante), Gris = residentes",
                             if(show_labels) paste0("\nEtiquetas: Top ", label_top_abundant, " abundantes + Top ", label_invaders, " invasores") else "\nSin etiquetas")
      
    } else {
      # Estándar: Gris para invasores, rainbow para abundancia
      ranks_numeric <- rank_groups[rank_groups != "Invader"] %>%
        sub("Rank_", "", .) %>%
        as.numeric()
      ranks_to_color <- rank_groups[rank_groups != "Invader"][order(ranks_numeric)]
      
      n_ranks <- length(ranks_to_color)
      color_scheme <- c(
        "Invader" = "gray50",
        setNames(rev(rainbow(n_ranks)), ranks_to_color)
      )
      
      plot_title <- "Comparación Rainbow (Ranking por Abundancia)"
      plot_subtitle <- paste0("Orden: ", paste(treatments, collapse = " -> "))
      plot_caption <- paste0("Colores: Rainbow = abundancia, Gris = invasores",
                             if(show_labels) paste0("\nEtiquetas: Top ", label_top_abundant, " abundantes + Top ", label_invaders, " invasores") else "\nSin etiquetas")
    }
    color_guide <- "none"
  }
  
  # Título común
  plot_title <- paste0(plot_title, 
                       if(!is.null(taxon_level)) paste0(" - ", taxon_level) else "",
                       " - ", paste(treatments, collapse = " vs "))
  
  # Crear el plot
  p <- ggplot(plot_data_long, aes(x = condition_jitter, y = abundance)) +
    # Borde blanco para highlight (si aplica)
    {if(highlight_invaders && !use_taxonomy_coloring) geom_point(aes(color = color_group, size = point_size + 0.5), color = "white", alpha = 1, shape = 16)} +
    # Líneas y puntos principales
    geom_line(aes(group = clade_name, color = line_color), alpha = 0.5, linewidth = 0.5) +
    geom_point(aes(color = color_group, size = point_size), alpha = 0.8) +
    scale_color_manual(
      values = color_scheme,
      na.value = "black", na.translate = FALSE, guide = color_guide
    ) +
    scale_size_continuous(range = c(0.5, 4), guide = "none") +
    scale_x_continuous(
      breaks = 1:length(treatments),
      labels = treatments
    ) +
    scale_y_log10(labels = scales::comma_format()) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Tratamiento",
      y = "Abundancia Relativa (log10)",
      caption = plot_caption
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
      legend.position = "right",
      legend.box.margin = margin(l = 10)
    )
  
  # Agregar etiquetas si se solicitan
  if (show_labels && exists("top_labels") && nrow(top_labels) > 0) {
    label_data <- top_labels %>%
      select(clade_name, label, all_of(treatments)) %>%
      pivot_longer(
        cols = all_of(treatments),
        names_to = "condition",
        values_to = "abundance"
      ) %>%
      mutate(
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
    
    filename <- paste0(
      if(use_taxonomy_coloring) paste0("taxonomy_filtered_", color_by_taxon_col) else paste0("rainbow", if(highlight_invaders) "_highlight_invaders" else ""),
      "_", paste(treatments, collapse = "_"),
      if(!is.null(taxon_level)) paste0("_", taxon_level) else "",
      ".png")
    
    output_file <- file.path(output_dir, filename)
    ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
    cat("Plot guardado:", basename(output_file), "\n")
  }
  
  return(p)
}