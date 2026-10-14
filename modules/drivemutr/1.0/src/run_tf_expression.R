log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(readr)
library(tibble)
library(ggplot2)

# The WGCNA normalized_expression.tsv is genes x samples with trailing `gene_id`
# and `hgnc_symbol` columns. Reshape it to a samples x genes matrix so per-TF means come straight from colMeans.
norm <- readr::read_tsv(snakemake@input[["normalized_expression"]], show_col_types = FALSE)
if (!"hgnc_symbol" %in% colnames(norm)) {
  stop("normalized_expression is missing the 'hgnc_symbol' column expected from the WGCNA module")
}
sample_cols    <- setdiff(colnames(norm), c("gene_id", "hgnc_symbol"))
gene_by_sample <- as.matrix(norm[, sample_cols, drop = FALSE])
rownames(gene_by_sample) <- norm[["hgnc_symbol"]]
expression_matrix <- t(gene_by_sample)  # samples x genes

tf_names             <- readr::read_lines(snakemake@input[["tf_names"]])
expression_threshold <- snakemake@params[["expression_threshold"]]

tf_exp   <- expression_matrix[, colnames(expression_matrix) %in% tf_names, drop = FALSE]
tf_means <- colMeans(tf_exp, na.rm = TRUE)

# TFs whose mean expression meets the threshold are considered "expressed"
tf_names_expressed <- names(tf_means)[tf_means >= expression_threshold]
readr::write_lines(tf_names_expressed, snakemake@output[["names"]])

message(sprintf("TFs in matrix: %d | expressed (>= %s): %d",
                length(tf_means), expression_threshold, length(tf_names_expressed)))

# Distribution of per-TF mean expression, with a dashed line at the threshold
plot_df <- tibble::tibble(tf_mean = as.numeric(tf_means))
p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = tf_mean)) +
  ggplot2::geom_histogram(bins = 40, fill = "#4ea8ff", colour = "white", alpha = 0.9) +
  ggplot2::geom_vline(xintercept = expression_threshold,
                      linetype = "dashed", colour = "#ef506a", linewidth = 0.8) +
  ggplot2::annotate("text", x = expression_threshold, y = Inf,
                    label = paste0("threshold = ", expression_threshold),
                    colour = "#ef506a", hjust = -0.05, vjust = 1.6, size = 3.6) +
  ggplot2::labs(
    x = "Mean normalized expression per TF",
    y = "Number of transcription factors",
    title = "Distribution of TF mean expression",
    subtitle = sprintf("%d / %d TFs above threshold",
                       length(tf_names_expressed), length(tf_means))
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(snakemake@output[["plot"]], p, width = 7, height = 5)