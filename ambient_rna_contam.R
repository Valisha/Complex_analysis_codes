
### Ambient RNA contamination

base_dir <- "/Users/valishashah/Library/CloudStorage/Box-Box/Kaech Lab Folder/Valisha/PDAC Diet CITESeq (Stephen)/DATA/250625 fish soy csf pdgf pdac scRNAseq2/counts_download_20251015/counts_download_20251015/"
pools <- c("Fish_WT1","Soy_WT1","Soy_iFat1_Pdgfra_cre1","Soy_iFat1_Csf1r_cre1",
           "Fish_WT2","Soy_WT2","Soy_iFat1_Pdgfra_cre2","Soy_iFat1_Csf1r_cre2")

## ---------- helpers ----------
# Ensure we always get the RNA matrix (not ADTs)
get_rna_matrix <- function(x) {
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) return(x[["Gene Expression"]])
    stop("Input returned a list but doesn’t contain 'Gene Expression'.")
  }
  x
}

read_raw_matrix <- function(sample_dir) {
  h5 <- file.path(sample_dir, "count", "sample_filtered_feature_bc_matrix.h5")
  if (file.exists(h5)) return(get_rna_matrix(Read10X_h5(h5, use.names=TRUE, unique.features=TRUE)))
  mtx <- file.path(sample_dir, "count", "sample_filtered_feature_bc_matrix")
  if (dir.exists(mtx)) return(get_rna_matrix(Read10X(mtx, gene.column=2)))
  alt <- list.files(file.path(sample_dir, "count"),
                    pattern="^sample_filtered_feature_bc_matrix", full.names=TRUE)
  if (length(alt)) {
    if (any(grepl("\\.h5$", alt)))
      return(get_rna_matrix(Read10X_h5(alt[grepl("\\.h5$", alt)][1], use.names=TRUE, unique.features=TRUE)))
    dirs <- alt[dir.exists(alt)]
    if (length(dirs)) return(get_rna_matrix(Read10X(dirs[1], gene.column=2)))
  }
  stop("Couldn’t find matrix (.h5 or folder).")
}

# Run EmptyNN on a cells×genes matrix, with guaranteed rownames
.run_emptynn_cellsxgenes <- function(cg, threshold=100, k=10, iteration=10, verbose=TRUE) {
  # Rownames must exist and correspond to cells
  if (is.null(rownames(cg))) {
    stop("cells×genes matrix has NULL rownames; cannot align results back to raw.")
  }
  fit <- emptynn(cg, threshold=threshold, k=k, iteration=iteration, verbose=verbose)
  fit
}

# Attempt both orientations, but ALSO create an index map from training rows -> raw cols
run_emptynn_orient <- function(raw_mat, threshold=100, k=10, iteration=10, verbose=TRUE) {
  # raw_mat is genes×cells (standard 10x). We want cells×genes.
  cellsxgenes <- if (nrow(raw_mat) > ncol(raw_mat)) t(raw_mat) else raw_mat
  
  # Make sure names are synced to raw cols (barcodes)
  if (is.null(rownames(cellsxgenes))) rownames(cellsxgenes) <- colnames(raw_mat)
  # And genes on columns:
  if (is.null(colnames(cellsxgenes))) colnames(cellsxgenes) <- rownames(raw_mat)
  
  # First try cells×genes
  try1 <- try(.run_emptynn_cellsxgenes(cellsxgenes, threshold, k, iteration, verbose), silent=TRUE)
  if (!inherits(try1, "try-error")) {
    return(list(nn=try1, mat=cellsxgenes, map_rows_to_rawcols=seq_len(ncol(raw_mat))))
  }
  
  # If that errored specifically asking to transpose, flip; otherwise surface the error
  msg <- as.character(attr(try1, "condition"))
  if (!grepl("Please transpose", msg, ignore.case=TRUE)) {
    stop("EmptyNN failed (cells×genes): ", msg)
  }
  
  message("EmptyNN requested transpose; retrying with flipped orientation ...")
  # Flipped: genes×cells (not what EmptyNN expects), so transpose back to cells×genes explicitly
  # We rebuild cells×genes to ensure rownames map to raw colnames
  cellsxgenes2 <- t(raw_mat)  # rows=cells (raw cols), cols=genes (raw rows)
  if (is.null(rownames(cellsxgenes2))) rownames(cellsxgenes2) <- colnames(raw_mat)
  if (is.null(colnames(cellsxgenes2))) colnames(cellsxgenes2) <- rownames(raw_mat)
  
  try2 <- try(.run_emptynn_cellsxgenes(cellsxgenes2, threshold, k, iteration, verbose), silent=TRUE)
  if (inherits(try2, "try-error")) {
    stop("EmptyNN failed on both orientations:\n", msg)
  }
  list(nn=try2, mat=cellsxgenes2, map_rows_to_rawcols=seq_len(ncol(raw_mat)))
}

# Convert nn.keep -> raw matrix COLUMN INDICES (not names)
nn_keep_to_raw_col_idx <- function(nn_fit, train_mat, map_rows_to_rawcols) {
  n_train <- nrow(train_mat)
  
  if (is.logical(nn_fit$nn.keep)) {
    stopifnot(length(nn_fit$nn.keep) == n_train)
    idx_in_train <- which(nn_fit$nn.keep)
  } else if (is.numeric(nn_fit$nn.keep)) {
    idx_in_train <- as.integer(nn_fit$nn.keep)
    stopifnot(all(idx_in_train >= 1 & idx_in_train <= n_train))
  } else if (!is.null(names(nn_fit$nn.keep))) {
    # Named vector: TRUE at kept names
    keep_names <- names(nn_fit$nn.keep)[nn_fit$nn.keep]
    # Match names to training rownames
    idx_in_train <- match(keep_names, rownames(train_mat))
    idx_in_train <- idx_in_train[!is.na(idx_in_train)]
  } else {
    stop("Unrecognized nn.keep format.")
  }
  
  # Map training-row indices -> original raw column indices
  raw_col_idx <- map_rows_to_rawcols[idx_in_train]
  raw_col_idx
}

## ---------- main ----------
base_dir <- "/Users/valishashah/Library/CloudStorage/Box-Box/Kaech Lab Folder/Valisha/PDAC Diet CITESeq (Stephen)/DATA/250625 fish soy csf pdgf pdac scRNAseq2/counts_download_20251015/counts_download_20251015/"
pools <- c("Fish_WT1","Soy_WT1","Soy_iFat1_Pdgfra_cre1","Soy_iFat1_Csf1r_cre1",
           "Fish_WT2","Soy_WT2","Soy_iFat1_Pdgfra_cre2","Soy_iFat1_Csf1r_cre2")

objs_emptynn <- list()

for (pool in pools) {
  sample_dir <- file.path(base_dir, pool)
  raw <- read_raw_matrix(sample_dir)   # genes×cells
  message("RAW dims (genes×cells): ", paste(dim(raw), collapse=" x"), " — ", pool)
  
  if (ncol(raw) < 500L || nrow(raw) < 1000L) {
    warning(sprintf("Skipping %s (too few cells or genes: %dx%d).", pool, nrow(raw), ncol(raw)))
    next
  }
  
  # Start with conservative params while debugging
  fit <- run_emptynn_orient(raw_mat=raw, threshold=100, k=5, iteration=5, verbose=TRUE)
  
  keep_idx <- nn_keep_to_raw_col_idx(fit$nn, fit$mat, fit$map_rows_to_rawcols)
  
  # If nothing kept, relax a bit and retry once
  if (length(keep_idx) == 0L) {
    message("No cells kept; retrying with a milder threshold...")
    fit <- run_emptynn_orient(raw_mat=raw, threshold=20, k=5, iteration=5, verbose=TRUE)
    keep_idx <- nn_keep_to_raw_col_idx(fit$nn, fit$mat, fit$map_rows_to_rawcols)
  }
  
  message(sprintf("Keeping %d / %d cells for %s", length(keep_idx), ncol(raw), pool))
  stopifnot(length(keep_idx) > 0L)
  
  counts_for_seurat <- raw[, keep_idx, drop=FALSE]
  retained <- runSeurat(counts=counts_for_seurat, resolution=0.2)
  print(DimPlot(retained, reduction="tsne") + ggtitle(sprintf("EmptyNN (raw): %s", pool)) + NoLegend())
  objs_emptynn[[pool]] <- retained
  
  # Hygiene on Apple Silicon
  try({ tensorflow::tf$keras$backend$clear_session() }, silent=TRUE)
  invisible(gc())
}
