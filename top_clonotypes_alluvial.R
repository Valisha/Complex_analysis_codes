plot_top_clonotypes_alluvial <- function(
    meta,
    origin_col      = "conditions_final",
    origin_levels   = c("Parenchyma", "Midline"),
    chain           = c("TRA", "TRB"),
    cdr3_col        = NULL,
    top_n           = 50,        # <--- now default 50
    gp33            = FALSE,
    gp_tra_seq      = NULL,
    gp_trb_seq      = NULL,
    title           = NULL,
    condition_palette = condition_palette
) {
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
  library(scales)
  library(ggVennDiagram)
  library(rlang)
  
  chain <- match.arg(chain)
  if (is.null(cdr3_col)) {
    cdr3_col <- paste0("cdr3_", chain)
  }
  
  # interpret gp33 flag
  gp33_flag <- FALSE
  if (is.logical(gp33)) gp33_flag <- isTRUE(gp33)
  else if (is.character(gp33)) gp33_flag <- tolower(gp33) %in% c("yes","true","gp33","gp-33")
  
  origin_sym <- sym(origin_col)
  cdr3_sym   <- sym(cdr3_col)
  
  # 1) subset to chosen origins + non-empty CDR3
  meta2 <- meta %>%
    filter(
      .data[[origin_col]] %in% origin_levels,
      !is.na(.data[[cdr3_col]]),
      .data[[cdr3_col]] != ""
    ) %>%
    mutate(!!origin_col := factor(.data[[origin_col]], levels = origin_levels))
  
  if (nrow(meta2) == 0) stop("No cells left after filtering.")
  
  # 2) compute *global* top-N clonotypes across all origins -------------------
  clone_totals_global <- meta2 %>%
    group_by(!!cdr3_sym) %>%
    summarise(total_cells = n(), .groups = "drop") %>%
    arrange(desc(total_cells))
  
  n_available <- nrow(clone_totals_global)
  top_ids <- clone_totals_global %>%
    slice_head(n = min(top_n, n_available)) %>%
    pull(!!cdr3_sym)
  
  # 3) counts per origin for those top clonotypes -----------------------------
  counts_top <- meta2 %>%
    filter(.data[[cdr3_col]] %in% top_ids) %>%
    dplyr::count(!!origin_sym, !!cdr3_sym, name = "n_cells")
  
  if (nrow(counts_top) == 0) stop("No counts to plot after top-N filtering.")
  
  # 4) gp33 tagging (no filtering)
  counts_top <- counts_top %>% mutate(is_gp33 = FALSE)
  
  # counts_top already has: origin, clonotype (cdr3), n_cells, is_gp33
  
  label_df <- counts_top %>%
    arrange(!!origin_sym, !!cdr3_sym) %>%
    group_by(!!origin_sym) %>%
    mutate(
      ymin = dplyr::lag(cumsum(n_cells), default = 0),
      ymax = cumsum(n_cells),
      y    = (ymin + ymax) / 2
    ) %>%
    ungroup()
  
  if (gp33_flag) {
    if (is.null(gp_tra_seq)) gp_tra_seq <- get0("gp_tra_seq", ifnotfound = NULL)
    if (is.null(gp_trb_seq)) gp_trb_seq <- get0("gp_trb_seq", ifnotfound = NULL)
    
    if (chain == "TRA") {
      if (is.null(gp_tra_seq)) stop("gp33=TRUE but gp_tra_seq is missing.")
      counts_top$is_gp33 <- counts_top[[cdr3_col]] %in% gp_tra_seq
      
    } else if (chain == "TRB") {
      if (is.null(gp_trb_seq)) stop("gp33=TRUE but gp_trb_seq is missing.")
      counts_top$is_gp33 <- counts_top[[cdr3_col]] %in% gp_trb_seq
    }
  } else {
    counts_top$is_gp33 <- TRUE  # show all
  }
  
  gp_ids <- sort(unique(counts_top[[cdr3_col]][counts_top$is_gp33]))
  if (gp33_flag && length(gp_ids) == 0) {
    warning("No gp33 clonotypes among top-N; showing regular plot.")
    counts_top$is_gp33 <- TRUE
  }
  
  all_ids <- sort(unique(counts_top[[cdr3_col]]))
  counts_top[[cdr3_col]] <- factor(counts_top[[cdr3_col]], levels = all_ids)
  
  palette_vec <- hue_pal()(length(all_ids))
  names(palette_vec) <- all_ids
  
  # title
  if (is.null(title)) {
    title <- paste0(
      "Top ", length(all_ids), " ", chain, " clonotypes across cohort (origins: ",
      paste(origin_levels, collapse = " vs "), ")",
      if (gp33_flag) " – gp33 highlighted" else ""
    )
  }
  # build legend labels, bold GP33 clonotypes
  legend_labels <- sapply(all_ids, function(id) {
    if (gp33_flag && id %in% gp_ids) {
      # bold label
      bquote(bold(.(id)))
    } else {
      id
    }
  })
  # PLOT ----------------------------------------------------------------------
  p <- ggplot(
    counts_top,
    aes(
      x        = !!origin_sym,
      stratum  = !!cdr3_sym,
      alluvium = !!cdr3_sym,
      y        = n_cells,
      fill     = !!cdr3_sym
    )
  ) +
    geom_flow(
      lode.guidance = "frontback",
      knot.pos      = 0.4,
      color         = "grey20",
      size          = 0.25,
      alpha         = 0.8      # <- constant, all flows visible
    ) +
    geom_stratum(width = 0.25, color = "white", alpha = 1) +
    geom_text(
      stat  = "stratum",
      aes(label = after_stat(stratum), alpha = is_gp33),
      size          = 2.5,
      color         = "black",
      check_overlap = TRUE
    ) +
    scale_fill_manual(
      values = palette_vec,
      breaks = all_ids,                       # show ALL sequences in legend
      labels = legend_labels[all_ids],
      name   = if (gp33_flag) paste0("gp33 ", chain, " clonotypes")
      else           paste0("Top ", chain, " clonotypes")
    ) +
    scale_alpha_manual(
      values = c(`FALSE` = 0.3, `TRUE` = 1),
      guide  = "none"
    ) +
    labs(x = NULL, y = "Number of cells", title = title) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(size = 11),
      axis.title.y = element_text(size = 12)
    )
  venn_list <- lapply(origin_levels, function(org) {
    meta2 %>%
      filter(.data[[origin_col]] == org) %>%
      pull(.data[[cdr3_col]]) %>%
      unique()
  })
  names(venn_list) <- origin_levels
  
  # vector of fill colors in correct order
  fill_cols <- condition_palette[origin_levels]
  
  # draw venn to a grid object
  venn_plot <- VennDiagram::venn.diagram(
    x = venn_list,
    filename = NULL,               # return a grob instead of saving
    fill = fill_cols,
    alpha = 0.6,
    col = "black",
    cat.col = "black",
    cat.cex = 1.3,
    cex = 1.2,
    margin = 0.08
  )
  
  # Convert to ggplot-friendly object if needed:
  library(grid)
  grid.newpage()
  grid.draw(venn_plot)
  
  return(list(
    alluvial = p,
    venn     = venn_plot,
    venn_data = venn_list
  ))
}








