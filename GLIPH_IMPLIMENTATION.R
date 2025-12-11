## Run GLIPH2

### ----- Extract CDR3 TRA / TRB from your merged object -----
library(dplyr)
library(turboGliph)

nrow(cd8_tcells_tcr@meta.data)
# length(cd8_tcells_tcr@meta.data$barcode_core)
length(cd8_tcells_tcr@meta.data$tcr_barcode_core)
head(cd8_tcells_tcr@meta.data$tcr_barcode_core)
head(cd8_tcells_tcr@meta.data$barcode)
head(rownames(cd8_tcells_tcr@meta.data))

gliph_input <- cd8_tcells_tcr@meta.data %>%
  mutate(
    cell_id = tcr_barcode_core,
    CDR3a   = cdr3_TRA,
    CDR3b   = cdr3_TRB,
    # GLIPH2 works best with V/J genes; you can keep them NA if you don’t have them
    TRAV = str_split_fixed(v_gene, ";", 3)[,1],
    TRAJ = str_split_fixed(j_gene, ";", 3)[,1],
    TRBV = str_split_fixed(v_gene, ";", 3)[,2],
    TRBJ = str_split_fixed(j_gene, ";", 3)[,2]
  ) %>%
  dplyr::select(cell_id, CDR3b, TRBV, TRBJ, CDR3a, TRAV, TRAJ)

# Keep only productive T cells with both alpha and beta if you want:
# gliph_input <- gliph_input %>% filter(!is.na(CDR3a) | !is.na(CDR3b))
gliph_input

# 
# # Save to GLIPH2 input format:
# write.table(gliph_input,
#             file = paste0(out_dir, "gliph2_input.txt"),
#             sep = "\t",
#             row.names = FALSE,
#             quote = FALSE)

# remotes::install_github("HetzDra/turboGliph")
results <- turboGliph::gliph2(
  cdr3_sequences =gliph_input,
  # cdr3_columns = c("CDR3b", "CDR3a"),   # supports paired chains
  result_folder = paste0(out_dir, "gliph2_output_20251210-v1")
  # do_vdj = TRUE,     # uses V/J genes if provided
  # do_gene = TRUE,
  # do_length = TRUE
)

results$connections

motifs <- read.table("~/Library/CloudStorage/Box-Box/Kaech Lab Folder/Valisha/AD Serial Infection Project  (Irene & Brian W)/Exp1+Exp2 SCRIPTS and FIGURES/OUTPUTS/gliph2_output_20251210-v1/all_motifs.txt", header=TRUE, sep="\t")
clusters <- read.table("~/Library/CloudStorage/Box-Box/Kaech Lab Folder/Valisha/AD Serial Infection Project  (Irene & Brian W)/Exp1+Exp2 SCRIPTS and FIGURES/OUTPUTS/gliph2_output_20251210-v1/cluster_member_details.txt", header=TRUE, sep="\t")

head(clusters)
results$cluster_properties
clusters

gliph_clusters <- clusters %>%
  dplyr::select(cell_id, CDR3b, seq_ID, tag, ultCDR3b) %>%
  mutate(
    GLIPH_cluster = paste0("clust_", seq_ID),
    GLIPH_motif   = tag,
    GLIPH_ultimate = ultCDR3b
  )

gliph_collapsed <- gliph_clusters %>%
  group_by(cell_id) %>%
  summarise(
    GLIPH_cluster = dplyr::first(unlist(GLIPH_cluster)),
    GLIPH_motif   = dplyr::first(unlist(GLIPH_motif)),
    GLIPH_rep     = dplyr::first(unlist(GLIPH_ultimate)),
    CDR3b         = dplyr::first(unlist(CDR3b)),
    .groups = "drop"
  )
dim(gliph_collapsed)

gliph_clusters %>%
  filter(cell_id == "AAACCCCAGCAGGGGT")
gliph_collapsed %>%
  filter(cell_id == "AAACCCCAGCAGGGGT")

length(gliph_collapsed$cell_id[gliph_collapsed$cell_id %in% rownames(cd8_tcells_tcr@meta.data)])

gliph_cells <- intersect(
  rownames(cd8_tcells_tcr@meta.data),
  gliph_collapsed$cell_id
)

length(gliph_cells)
# should be 1796

cd8_tcells_gliph <- subset(cd8_tcells_tcr, cells = gliph_cells)

ncol(cd8_tcells_gliph)
# 1796

nrow(cd8_tcells_gliph@meta.data)
# 1796

identical(Cells(cd8_tcells_gliph),
          rownames(cd8_tcells_gliph@meta.data))
head(Cells(cd8_tcells_gliph))

meta_gliph <- cd8_tcells_gliph@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  left_join(gliph_collapsed, by = "cell_id") %>%
  tibble::column_to_rownames("cell_id")

cd8_tcells_gliph@meta.data <- meta_gliph

cd8_tcells_gliph$GLIPH_cluster <- meta_gliph$GLIPH_cluster
cd8_tcells_gliph$GLIPH_motif   <- meta_gliph$GLIPH_motif
cd8_tcells_gliph$GLIPH_rep     <- meta_gliph$GLIPH_rep

stopifnot(
  identical(Cells(cd8_tcells_gliph), rownames(cd8_tcells_gliph@meta.data)),
  identical(Cells(cd8_tcells_gliph), rownames(cd8_tcells_gliph@reductions$umap@cell.embeddings))
)

length(unique(cd8_tcells_gliph@meta.data$GLIPH_motif))
gliph_umap <-  DimPlot(cd8_tcells_gliph, group.by="GLIPH_motif", label=TRUE)

ggsave(
  filename = file.path(
    out_dir,
    sprintf("gliph_umap.png")
  ),
  plot   = gliph_umap,
  device = "png",
  width  = 40,
  height = 40,
  units  = "in"
)

# cd8_tcells_tcr$GLIPH_cluster <- cd8_tcells_tcr@meta.data$GLIPH_cluster
# cd8_tcells_tcr$GLIPH_motif   <- cd8_tcells_tcr@meta.data$GLIPH_motif
# cd8_tcells_tcr$GLIPH_rep     <- cd8_tcells_tcr@meta.data$GLIPH_ultimate
# # after all your joins:
# rownames(cd8_tcells_tcr@meta.data) <- colnames(cd8_tcells_tcr)
# rownames(cd8_tcells_tcr@meta.data)[grep("_", rownames(cd8_tcells_tcr@meta.data))]
# 
# DimPlot(cd8_tcells_tcr, group.by="GLIPH_cluster", label=TRUE, repel=TRUE)
# DimPlot(cd8_tcells_tcr, group.by="GLIPH_motif", label=TRUE)




