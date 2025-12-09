### Function to create the dotplot with the DEGs for each cluster comparing one condition vs the rest 


deg_dotPlot_per_Cluster_bidir <- function(
    tcells,
    cluster,
    n_top = 25,
    pairwise_list,
    expr_assay = c("RNA", "SCT"),
    expr_slot  = c("scale.data", "data")
){
  expr_assay <- match.arg(expr_assay)
  expr_slot  <- match.arg(expr_slot)
  
  message("Processing cluster: ", cluster)
  
  # ================================
  # 0. BUILD PAIRWISE COMBINATIONS
  # ================================
  build_pairs_from_input <- function(x) {
    # Case 1: user passes a character vector directly
    if (is.character(x)) {
      conds <- x
      
      # Case 2: user passes list(c("A","B","C","D"))
    } else if (is.list(x) && length(x) == 1 && is.character(x[[1]]) && length(x[[1]]) > 2) {
      conds <- x[[1]]
      
      # Case 3: user already passes list of pairs
    } else if (is.list(x) && all(vapply(x, length, FUN.VALUE = 1L) == 2)) {
      return(x)
      
    } else {
      stop("pairwise_list must be either:
           (1) a character vector of conditions,
           (2) list(c('A','B','C',...)),
           or (3) a list of length-2 character vectors.")
    }
    
    # Generate all unordered pairs of conditions
    combn(conds, 2, simplify = FALSE)
  }
  
  pairwise_pairs <- build_pairs_from_input(pairwise_list)
  
  # ================================
  # 1. SUBSET CLUSTER
  # ================================
  sub_obj <- subset(tcells, subset = tcells_custom_annotations == cluster)
  DefaultAssay(sub_obj) <- expr_assay
  
  if (expr_assay == "RNA") {
    sub_obj <- FindVariableFeatures(sub_obj)
    sub_obj <- ScaleData(sub_obj, features = rownames(sub_obj))
  }
  if (expr_assay == "SCT") {
    sub_obj <- SCTransform(sub_obj)
    sub_obj <- PrepSCTFindMarkers(sub_obj)
  }
  
  # Use conditions_final as the identity
  Idents(sub_obj) <- sub_obj$conditions_final
  
  # ================================
  # 2. DEG IN BOTH DIRECTIONS
  # ================================
  all_markers <- list()
  
  for (pair in pairwise_pairs) {
    a <- pair[1]
    b <- pair[2]
    
    # ------- A up vs B -------
    comp1 <- paste0(a, " up vs ", b)
    deg1 <- FindMarkers(
      sub_obj,
      assay = expr_assay,
      slot  = expr_slot,
      ident.1 = a,   # up in A
      ident.2 = b,
      only.pos = TRUE,
      logfc.threshold = 0.25,
      min.pct = 0.1,
      return.thresh = 0.05
    )
    deg1$gene_symbol <- rownames(deg1)
    deg1$comparison  <- comp1
    all_markers[[comp1]] <- deg1
    
    # ------- B up vs A -------
    comp2 <- paste0(b, " up vs ", a)
    deg2 <- FindMarkers(
      sub_obj,
      assay = expr_assay,
      slot  = expr_slot,
      ident.1 = b,   # up in B
      ident.2 = a,
      only.pos = TRUE,
      logfc.threshold = 0.25,
      min.pct = 0.1,
      return.thresh = 0.05
    )
    deg2$gene_symbol <- rownames(deg2)
    deg2$comparison  <- comp2
    all_markers[[comp2]] <- deg2
  }
  
  markers_pw_df <- dplyr::bind_rows(all_markers)
  
  # ================================
  # 3. TOP GENES PER DIRECTION
  # ================================
  markers_pw_top <- markers_pw_df %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    dplyr::group_by(comparison) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = n_top) %>%
    dplyr::ungroup()
  # print(markers_pw_top)
  comparison_levels <- unique(markers_pw_top$comparison)
  
  gene_order <- markers_pw_top %>%
    dplyr::arrange(match(comparison, comparison_levels), dplyr::desc(avg_log2FC)) %>%
    dplyr::pull(gene_symbol)
  
  gene_set <- unique(gene_order)
  # ================================
  # 4. AVERAGE EXPRESSION IN THE “UP” CONDITION
  # ================================
  get_pairwise_avg_expr <- function(obj, genes, cond2, comp_name){
    Idents(obj) <- obj$conditions_final
    cells2 <- WhichCells(obj, idents = cond2)
    mat <- GetAssayData(obj, assay = expr_assay, slot = expr_slot)
    avg2 <- Matrix::rowMeans(mat[genes, cells2, drop = FALSE])
    
    tibble::tibble(
      gene_symbol = genes,
      comparison  = comp_name,
      avg_expr    = avg2
    )
  }
  
  expr_list <- list()
  
  for (pair in pairwise_pairs){
    a <- pair[1]
    b <- pair[2]
    
    comp1 <- paste0(a, " up vs ", b)
    comp2 <- paste0(b, " up vs ", a)
    
    # NEW: expression from the UP condition
    expr_list[[comp1]] <- get_pairwise_avg_expr(sub_obj, gene_set, cond2 = a, comp_name = comp1)
    expr_list[[comp2]] <- get_pairwise_avg_expr(sub_obj, gene_set, cond2 = b, comp_name = comp2)
  }
  
  avg_expr_long <- dplyr::bind_rows(expr_list)
  
  # ================================
  # 5. MERGE + PREP DATA
  # ================================
  full_df <- markers_pw_df %>%
    dplyr::select(gene_symbol, comparison, avg_log2FC) %>%
    dplyr::right_join(avg_expr_long, by = c("gene_symbol","comparison")) %>%
    dplyr::mutate(
      avg_log2FC = ifelse(is.na(avg_log2FC), 0, avg_log2FC),
      dot_size   = avg_log2FC,  # positive only (only.pos = TRUE)
      comparison = factor(comparison, levels = comparison_levels),
      gene_symbol = factor(gene_symbol, levels = gene_set),
      y_num = as.numeric(gene_symbol)
    ) %>%
    dplyr::filter(!is.na(comparison))
  
  expr_min <- min(full_df$avg_expr, na.rm = TRUE)
  expr_max <- max(full_df$avg_expr, na.rm = TRUE)
  
  # ================================
  # 6. DOT PLOT
  # ================================
  p <- ggplot2::ggplot(full_df, ggplot2::aes(x = comparison, y = gene_symbol)) +
    ggplot2::geom_point(ggplot2::aes(color = avg_expr, size = dot_size)) +
    ggplot2::scale_size_continuous(range = c(1, 10), name = "log2FC (up)") +
    ggplot2::scale_color_gradientn(
      colours = c("navy", "royalblue", "white", "orange", "red"),
      limits = c(expr_min, expr_max),
      name = paste0(expr_assay, " ", expr_slot)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle(paste("Bidirectional Pairwise DEG Dot Plot —", cluster))
  
  markersdf_return <- markers_pw_top %>%
    mutate(log10pvalue = -log10(p_val_adj)) %>%
    dplyr::select(avg_log2FC, log10pvalue, gene_symbol, comparison)
  
  return(list(plot = p, 
              markers_df = markers_pw_df))
}
############ HEATMAP
deg_heatmap_per_cluster <- function(tcells,
                                    cluster,
                                    n_top = 25, 
                                    pairwise_list) {
  
  message("Processing cluster: ", cluster)
  
  # -----------------------------
  # 1. Subset cluster
  # -----------------------------
  sub_obj <- subset(tcells, subset = tcells_custom_annotations == cluster)
  DefaultAssay(sub_obj) <- "RNA"
  
  sub_obj <- NormalizeData(sub_obj)
  sub_obj <- FindVariableFeatures(sub_obj, nfeatures = 2000)
  sub_obj <- ScaleData(sub_obj, features = rownames(sub_obj))
  
  # -----------------------------
  # 2. Compute DEGs
  # -----------------------------
  all_markers <- list()
  
  for (pair in pairwise_list) {
    a <- pair[1]; b <- pair[2]
    comp_name <- paste0(a, " vs ", b)
    
    Idents(sub_obj) <- sub_obj$conditions_final
    
    degs <- FindMarkers(
      sub_obj,
      assay = "RNA",
      slot  = "scale.data",
      ident.1 = a,
      ident.2 = b,
      only.pos = FALSE,
      logfc.threshold = 0.25,
      min.pct = 0.1,
      return.thresh = 0.05
    )
    
    degs$gene_symbol <- sub("\\.[0-9]+$", "", rownames(degs))
    degs$comparison  <- comp_name
    
    all_markers[[comp_name]] <- degs
  }
  
  markers_pw_df <- bind_rows(all_markers)
  
  # -----------------------------
  # 3. Select Top Genes + Order
  # -----------------------------
  markers_pw_top <- markers_pw_df %>%
    filter(p_val_adj < 0.05) %>%
    group_by(comparison) %>%
    slice_max(order_by = abs(avg_log2FC), n = n_top, with_ties = FALSE) %>%
    ungroup()
  
  comparison_levels <- unique(markers_pw_top$comparison)
  
  gene_order <- markers_pw_top %>%
    arrange(match(comparison, comparison_levels), desc(abs(avg_log2FC))) %>%
    pull(gene_symbol)
  
  gene_set <- unique(gene_order)
  
  # -----------------------------
  # 4. Compute Avg Expr (ident.2)
  # -----------------------------
  get_pairwise_avg_expr <- function(obj, genes, cond1, cond2, comp_name) {
    Idents(obj) <- obj$conditions_final
    cells2 <- WhichCells(obj, idents = cond2)
    
    mat <- GetAssayData(obj, assay="RNA", layer="scale.data")
    
    avg2 <- Matrix::rowMeans(mat[genes, cells2, drop = FALSE])
    
    tibble(
      gene_symbol = genes,
      comparison  = comp_name,
      avg_expr    = avg2,
      group_label = comp_name
    )
  }
  
  expr_list <- list()
  for (pair in pairwise_list) {
    a <- pair[1]; b <- pair[2]
    comp_name <- paste0(a, " vs ", b)
    expr_list[[comp_name]] <- get_pairwise_avg_expr(sub_obj, gene_set, a, b, comp_name)
  }
  
  avg_expr_long <- bind_rows(expr_list)
  
  # -----------------------------
  # 5. Merge DE + expression
  # -----------------------------
  full_df <- markers_pw_df %>%
    dplyr::select(gene_symbol, comparison, avg_log2FC) %>%
    right_join(avg_expr_long, by = c("gene_symbol", "comparison")) %>%
    mutate(
      avg_log2FC  = ifelse(is.na(avg_log2FC), 0, avg_log2FC),
      gene_symbol = factor(gene_symbol, levels = gene_set),
      comparison  = factor(comparison, levels = comparison_levels)
    )
  
  # -----------------------------
  # 6. Heatmap + Right Annotation Strip
  # -----------------------------
  # Create a categorical palette for groups
  # categorical palette for groups
  group_colors <- setNames(
    RColorBrewer::brewer.pal(length(comparison_levels), "Set2"),
    comparison_levels
  )
  
  p <- ggplot(full_df, aes(x = comparison, y = gene_symbol)) +
    
    # heatmap (continuous)
    geom_tile(aes(fill = avg_expr), color="grey80") +
    
    # log2FC numbers
    geom_text(aes(label = sprintf("%.2f", avg_log2FC)), size = 2.5) +
    
    # right-side annotation strip (discrete)
    geom_tile(
      aes(
        x = length(comparison_levels) + 1,
        y = gene_symbol,
        fill = group_label
      ),
      color = "white"
    ) +
    
    # continuous scale for heatmap
    scale_fill_gradientn(
      colours = c("navy", "royalblue", "white", "orange", "red"),
      name = "Scaled expr (RNA)",
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      guide = guide_colorbar(order = 1)
    ) +
    
    # discrete scale for grouping strip
    scale_fill_manual(
      values = group_colors,
      name = "Comparison group",
      aesthetics = "fill",
      guide = guide_legend(override.aes = list(size = 5), order = 2)
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=45, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    ) +
    ggtitle(paste("Pairwise DEG Heatmap —", cluster))
  return(p)
}



