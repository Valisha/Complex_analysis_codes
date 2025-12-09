######## clonotype tracking plot

library(ggalluvial)
library(ggtext)
library(dplyr)
library(scales)

plot_clonotype_alluvial <- function(
    seurat_object,
    condition_col = "conditions",     # e.g. "origin_new"
    clonotype_col = "cdr3_TRA",       # cdr3_TRA or cdr3_TRB
    clone_bin_col = "clone_bin_tra",  # for ranking clonotypes
    cond_levels,                      # conditions to include
    reference_condition,              # DO NOT remove this — source condition
    topN = 50,
    legend_false = "None",
    title_text,
    show_sequences = FALSE   # NEW 🔥
) {
  
  meta <- seurat_object@meta.data
  
  # ---------------------------
  # 1. Filter metadata to relevant conditions
  # ---------------------------
  meta2 <- meta %>%
    filter(
      .data[[condition_col]] %in% cond_levels,
      !is.na(.data[[clonotype_col]]),
      .data[[clonotype_col]] != ""
    ) %>%
    mutate(
      !!condition_col := factor(.data[[condition_col]], levels = cond_levels)
    )
  
  # ---------------------------
  # 2. Select top N UNIQUE clonotypes from the reference condition
  # ---------------------------
  top_clonotypes <- meta2 %>%
    filter(.data[[condition_col]] == reference_condition) %>%
    group_by(.data[[clonotype_col]]) %>%  
    summarise(total_size = sum(.data[[clone_bin_col]]), .groups = "drop") %>%
    arrange(desc(total_size)) %>%
    slice_head(n = topN) %>%
    pull(.data[[clonotype_col]])
  top_core <- ifelse(nchar(top_clonotypes) >= 3,
                     str_sub(top_clonotypes, 2, -2),
                     NA_character_)
  
  hits <- top_core[top_core %in% gp33_final$`Chain 1 | CDR3 Curated`]
  
  # Map back original clonotypes whose core matches GP33 curated
  gp33_clonotypes <- top_clonotypes[
    !is.na(top_core) &
      top_core %in% gp33_final$`Chain 1 | CDR3 Curated`
  ]
  gp33_clonotypes_comp <- meta %>%
    filter(`Chain 1 | CDR3 Curated` %in% gp33_clonotypes) %>%
    pull(.data[[clonotype_col]])
  
  if (length(top_clonotypes) == 0) {
    stop("No clonotypes available in reference condition for ranking.")
  }
  
  # ---------------------------
  # 3. Keep only those clonotypes across ALL selected conditions
  # ---------------------------
  meta3 <- meta2 %>%
    filter(.data[[clonotype_col]] %in% top_clonotypes)
  
  if (nrow(meta3) == 0) {
    stop("The selected top clonotypes do not appear in other conditions.")
  }
  
  # ---------------------------
  # 4. Count frequencies
  # ---------------------------
  counts <- meta3 %>%
    dplyr::count(.data[[condition_col]], .data[[clonotype_col]], name = "Frequency")
  
  # ---------------------------
  # 5. Compute GLOBAL proportions (same denominator for all conditions)
  # ---------------------------
  total_cells <- sum(counts$Frequency)
  counts <- counts %>%
    mutate(Proportion = Frequency / total_cells)
  
  # ---------------------------
  # 6. Identify clonotypes shared across all selected conditions
  # ---------------------------
  shared_clonotypes <- counts %>%
    group_by(.data[[clonotype_col]]) %>%
    summarise(n_cond = n_distinct(.data[[condition_col]]), .groups = "drop") %>%
    filter(n_cond == length(cond_levels)) %>%
    pull(.data[[clonotype_col]])
  
  # ---------------------------
  # 7. Palette
  # ---------------------------
  all_ids <- sort(unique(counts[[clonotype_col]]))
  pal <- hue_pal()(length(all_ids))
  names(pal) <- all_ids
  
  # bold these clonotypes:
  # 1. clonotypes shared across all conditions
  # 2. clonotypes whose core matches GP33 curated sequences
  # bold_ids <- union(shared_clonotypes, gp33_clonotypes)
  bold_ids <- gp33_clonotypes
  print(bold_ids)
  legend_labels <- ifelse(
    all_ids %in% bold_ids,
    paste0("**", all_ids, "**"),   # markdown bold
    all_ids
  )
  
  # # ---------------------------
  # # 8. Title
  # # ---------------------------
  # if (is.null(title_text)) {
  #   title_text <- sprintf(
  #     "Top %d clonotypes from %s tracked across: %s",
  #     topN,
  #     reference_condition,
  #     paste(cond_levels, collapse = " → ")
  #   )
  # }
  
  
  # ---------------------------
  # 9. Final Alluvial Plot
  # ---------------------------
  p <- ggplot(
    counts,
    aes(
      x = .data[[condition_col]],
      stratum  = .data[[clonotype_col]],
      alluvium = .data[[clonotype_col]],
      y        = Proportion,
      fill     = .data[[clonotype_col]]
    )
  ) +
    geom_flow(stat = "alluvium", color = "grey25", size = 0.3) +
    geom_stratum(width = 0.28, color = "white") +
    scale_fill_manual(
      values = pal,
      breaks = all_ids,
      labels = legend_labels,
      guide = guide_legend(override.aes = list(size = 5))
    ) +
    labs(x = NULL, y = "Proportion", title = title_text) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 12, face = "bold", angle = 45, vjust = 1, hjust = 1),
      axis.text.y  = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = legend_false,
      legend.text = ggtext::element_markdown(size = 8)
    )
  if (show_sequences) {
    
    if (show_sequences) {
      p <- p +
        geom_text(
          stat  = "stratum",
          aes(label = after_stat(stratum)),  # stratum = clonotype
          size = 2.5,
          color = "black",
          fontface = "bold",
          check_overlap = TRUE
        )
    }
  }
  
  # return list
  return(list(
    plot = p,
    top_clonotypes = top_clonotypes,
    shared_clonotypes = shared_clonotypes,
    palette = pal
  ))
}





