######## Violin plots with wilcoxon test statistics #############


plot_violin <- function(df, gene, pvals, cluster_palette){
  
  ymax <- max(df$expr, na.rm = TRUE)
  ystart <- ymax + 0.25    # more space
  step <- 0.25             # vertical spacing between comparisons
  bracket_height <- 0.07   # vertical offset for bracket ends
  
  g <- ggplot(df, aes(x = condition, y = expr, fill = condition)) +
    geom_violin(trim = TRUE, scale = "width", color = "black") +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    scale_fill_manual(values = cluster_palette) +
    theme_classic(base_size = 14) +
    ggtitle(paste(gene, "Expression by Cluster")) +
    ylab("Expression Level") +
    xlab("") +
    coord_cartesian(ylim = c(0, ystart + nrow(pvals)*step)) +
    theme(
      axis.text.x = element_text(angle = 0, face="bold"),
      plot.title = element_text(size=18, face="bold")
    )
  
  # --- Add p-value labels + brackets ---
  for(i in seq_len(nrow(pvals))){
    
    row <- pvals[i,]
    x1 <- which(levels(df$condition) == row$group1)
    x2 <- which(levels(df$condition) == row$group2)
    xmid <- mean(c(x1, x2))
    y <- ystart + (i-1)*step
    
    # Text label
    g <- g + annotate(
      "text",
      x = xmid,
      y = y + bracket_height + 0.1,
      label = row$label,
      size = 5,
      fontface = "bold"
    )
    
    # Bracket lines
    g <- g +
      annotate("segment", x = x1, y = y, xend = x2, yend = y, linewidth = 0.7) +  # horizontal
      annotate("segment", x = x1, y = y, xend = x1, yend = y - bracket_height, linewidth = 0.7) +  # left vertical
      annotate("segment", x = x2, y = y, xend = x2, yend = y - bracket_height, linewidth = 0.7)    # right vertical
  }
  
  return(g)
}





