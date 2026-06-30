# expression_matrix is samples x genes
filter_high_variance_genes <- function(expression_matrix,
                                       median_threshold = 4,
                                       mad_threshold = 0.25,
                                       plots_dir){
  
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  # Convert genes × samples input to the samples × genes matrix expected below.
  # If a gene-name column ("hgnc_symbol" or "gene") is present, use it as
  # rownames and keep only numeric columns before transposing.
  # If neither column exists, transpose as-is.
  gene_col <- intersect(c("hgnc_symbol", "gene"), colnames(expression_matrix))
  if (length(gene_col) > 0) {
    gene_col <- gene_col[1]   # prefer "hgnc_symbol" over "gene"
    expression_matrix <- expression_matrix %>%
      dplyr::select(all_of(gene_col), where(is.numeric)) %>%
      tibble::column_to_rownames(var = gene_col) %>%
      as.matrix() %>%
      t()
  } else {
    expression_matrix <- t(as.matrix(expression_matrix))
  }
  
  all_exp_vec <- as.vector(as.matrix(expression_matrix))
  
  plot_path <- file.path(plots_dir, "Primary_Distribution.pdf")
  
  grDevices::pdf(plot_path, height = 10, width = 10)
  
  graphics::hist(
    all_exp_vec,
    breaks = 100,
    xlim = c(3, 16),
    main = "Histogram of Distribution of Gene Expression Values",
    xlab = "Expression",
    col = "skyblue",
    border = "white"
  )
  
  grDevices::dev.off()
  
  median_genes = matrixStats::colMedians(as.matrix(expression_matrix))
  df_median = data.frame(median = median_genes)
  
  median_plot <- ggplot2::ggplot(df_median, aes(x = median)) +
    geom_histogram(bins = 100, fill = "lightblue", color = "white") +
    coord_cartesian(xlim = c(3, 14)) + # Optional zoom
    labs(title = "Distribution of the Median of Each Gene Expression",
         x = "Median (per gene)",
         y = "Count of Genes") +
    theme_minimal()
  
  plot_path <- file.path(plots_dir, "Median_Distribution_Before_Filtering.pdf")
  
  grDevices::pdf(plot_path, height = 10, width = 10)
  print(median_plot)
  grDevices::dev.off()
  
  
  all_exp_median_filtered <- expression_matrix[, median_genes >= median_threshold]
  
  mad_gene <- matrixStats::colMads(as.matrix(all_exp_median_filtered))
  
  df_mad <- data.frame(mad = mad_gene)
  
  mad_plot <- ggplot2::ggplot(df_mad, aes(x = mad)) +
    geom_histogram(bins = 100, fill = "lightblue", color = "white") +
    coord_cartesian(xlim = c(0, 2.5)) +  # Optional zoom
    labs(title = "Distribution of Gene Expression MAD",
         x = "MAD (per gene)",
         y = "Count of Genes") +
    theme_minimal()
  
  plot_path <- file.path(plots_dir, "MAD_Distribution_Before_MAD_After_Median_Filtering.pdf")
  
  grDevices::pdf(plot_path, height = 10, width = 10)
  print(mad_plot)
  grDevices::dev.off()
  
  all_exp_median_mad_filtered <- all_exp_median_filtered[, mad_gene >= mad_threshold]
  
  all_exp_median_mad_vec <- as.vector(as.matrix(all_exp_median_mad_filtered))
  
  
  plot_path <- file.path(plots_dir, "Final_Filtered_Distribution.pdf")
  
  grDevices::pdf(plot_path, height = 10, width = 10)
  
  graphics::hist(all_exp_median_mad_vec,
                 breaks = 50,
                 xlim = c(3, 16),
                 main = "Histogram of Distribution of Gene Expression Values",
                 xlab = "Expression",
                 col = "skyblue",
                 border = "white")
  
  grDevices::dev.off()
  
  return(all_exp_median_mad_filtered)
  
}