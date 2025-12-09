############# DEG heatmap for seurat clusters ###############


sct_heatmap_one_resolution <- function(
    obj,
    resolution,                      # e.g. "0.25"
    assay            = "SCT",
    de_slot          = "data",
    test.use         = "wilcox",
    min.pct          = 0.25,
    logfc.threshold  = 0.25,
    top_n_per_cluster= 20,
    top_table_n      = 50,           # for Excel wide table
    add_gaps         = TRUE,         # insert NA rows between cluster blocks in heatmap
    group_colors     = NULL,         # named vector cluster->color; NULL = auto
    out_dir          = "~/Library/CloudStorage/Box-Box/Kaech Lab Folder/Valisha/SK_EXPBWC015/OUTPUTS/",
    basename_prefix  = "Tcells",
    prep_sct         = TRUE          # call PrepSCTFindMarkers inside
) {
  tcells$seurat_clusters <- factor(tcells@meta.data[[resolution]])
  
  # Lock Idents to seurat_clusters (no ambiguity)
  Idents(tcells) <- tcells$seurat_clusters
  table(Idents(tcells))  # should match table(seurat_obj$seurat_clusters)
  table(tcells@meta.data$seurat_clusters)
  
  # required for DE on SCT when multiple models exist
  tcells <- PrepSCTFindMarkers(tcells, assay = "SCT")
  
  markers_all <- FindAllMarkers(
    tcells,
    assay          = "SCT",
    slot           = "data",   # use 'data' for SCT-based DE
    only.pos       = TRUE,
    test.use       = "wilcox",
    min.pct        = 0.25,
    logfc.threshold= 0.25
  )
  
  # Handle Seurat version differences in logFC column name
  lfc_col <- if ("avg_log2FC" %in% colnames(markers_all)) "avg_log2FC" else "avg_logFC"
  
  # Top 50 genes per cluster (by |logFC|, breaking ties by p_val_adj)
  top50_by_cluster <- markers_all %>%
    group_by(cluster) %>%
    arrange(desc(abs(.data[[lfc_col]])), p_val_adj, .by_group = TRUE) %>%
    slice_head(n = 50) %>%
    ungroup()
  
  # Assume you already have top10_features grouped per cluster
  top10_per_cluster <- top50_by_cluster %>%
    group_by(cluster) %>%
    slice_head(n = 20) %>%
    pull(gene)
  
  # Ensure valid genes
  top10_per_cluster <- intersect(top10_per_cluster, rownames(tcells))
  
  # Order features by cluster, then add spacer rows (fake "genes")
  features_grouped <- top50_by_cluster %>%
    group_by(cluster) %>%
    slice_head(n = 20) %>%
    mutate(order = row_number()) %>%
    arrange(cluster, order) %>%
    pull(gene)
  
  # Insert dummy separators between each cluster’s block
  features_with_gaps <- unlist(lapply(split(features_grouped, 
                                            top50_by_cluster$cluster[
                                              match(features_grouped, top50_by_cluster$gene)
                                            ]),
                                      function(x) c(x, NA)))  # NA adds a gap
  
  # Draw heatmap with custom cluster colors
  p_heatmap <- DoHeatmap(
    tcells,
    features = features_with_gaps,
    group.by = "seurat_clusters",
    assay    = "SCT",
    size     = 3,
    raster   = FALSE,
    group.colors = ss   # <---- custom palette applied here
  ) + 
    scale_fill_gradientn(
      colors   = c("yellow", "red"), 
      na.value = "white"  # NA separators show as white gaps
    ) +
    ggtitle("Top 20 DEGs per cluster (with res=0.25)")
  
  print(p_heatmap)

}










