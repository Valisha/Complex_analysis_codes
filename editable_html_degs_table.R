######## editable html DEGs ###########

merged_results_whole <- readRDS('/Users/valishashah/Library/CloudStorage/Box-Box/Kaech Lab Folder/Valisha/Arvinas (Tom)/LCMV_Exp2/SCRIPTS/merged_results_whole_DEG_CD4vsCD8.Rds')

editable_html_table <- function(deg_results, id = "de_genes", cluster) {
  # 1) Construct the exact column names
  log2fc_col  <- paste0("log2FC_",        cluster)
  padj_col    <- paste0("padj_",          cluster)
  neglog_col  <- paste0("negLog10_padj_", cluster)
  avg_col     <- paste0("avgExpr_",       cluster)
  
  # 2) Sanity‐check
  needed <- c(log2fc_col, padj_col, neglog_col, avg_col)
  if (!all(needed %in% colnames(deg_results))) {
    stop("Missing columns: ", paste(setdiff(needed, colnames(deg_results)), collapse = ", "))
  }
  
  # 3) Subset + round to 2 decimals
  df <- as.data.frame(deg_results) %>%
    filter(!is.na(.data[[padj_col]])) %>%
    mutate(across(all_of(needed), ~ round(as.numeric(.), 2)))
  
  # 4) Column definitions (styling only)
  cols <- list()
  cols[[log2fc_col]] <- colDef(cell = function(x) x, style = function(x) {
    col <- ifelse(x >  0.5, "#008000", ifelse(x < -0.5, "#e00000", "#777"))
    list(color = col)
  })
  cols[[padj_col]]   <- colDef(cell = function(x) x, style = function(x) list(fontWeight = ifelse(x < 0.05, "bold", "")))
  cols[[neglog_col]] <- colDef(cell = function(x) x, style = function(x) list(fontWeight = ifelse(x > 1.3, "bold", "")))
  cols[[avg_col]]    <- colDef(cell = function(x) x)
  
  # 5) Render an *editable* reactable
  htmltools::tagList(
    reactable(
      df,
      elementId     = id,
      columns       = cols,
      defaultSorted = neglog_col,
      filterable    = TRUE,
      striped       = TRUE,
      highlight     = TRUE),
    htmltools::tags$button(
      "Download as CSV",
      onclick = sprintf(
        "Reactable.downloadDataCSV('%s','%s_results.csv')",
        id, id
      )
    )
  )
}

for (tbl_name in names(merged_results_whole)) {
  cluster <- sub("^clust_", "", tbl_name)
  print(cluster)
  id      <- paste0("de_genes_",
                    gsub("[^A-Za-z0-9]", "_", tbl_name))
  
  print(editable_html_table(
    deg_results = merged_results_whole[[tbl_name]],
    id          = id,
    cluster     = cluster
  ))
}


