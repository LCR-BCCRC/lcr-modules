get_co_expressed_genes_modules <- function(expression_matrix, powers = 1:15, cor_method = "bicor", network_type = "signed hybrid", plots_dir){
  
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  sft_filtered_median_mad <- pickSoftThreshold(expression_matrix,
                                                        powerVector = powers,
                                                        corFnc = cor_method,
                                                        networkType = network_type)
  
  opt_power = sft_filtered_median_mad$fitIndices %>% filter(SFT.R.sq >= 0.8) %>% slice_min(SFT.R.sq) %>% pull(Power)
  
  soft_threshold_plot_path <- file.path(
    plots_dir,
    "Soft_Threshold_Diagnostic_Plots.pdf"
  )
  
  grDevices::pdf(
    soft_threshold_plot_path,
    height = 6,
    width = 12
  )
  
  graphics::par(
    mfrow = c(1, 2),
    mar = c(5, 5, 4, 2) + 0.1
  )
  
  # Plot 1: Scale independence
  graphics::plot(
    sft_filtered_median_mad$fitIndices[, 1],
    sft_filtered_median_mad$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale-Free Topology Model Fit, signed R²",
    type = "n",
    main = "Scale independence"
  )
  
  graphics::text(
    sft_filtered_median_mad$fitIndices[, 1],
    sft_filtered_median_mad$fitIndices[, 2],
    labels = sft_filtered_median_mad$fitIndices[, 1],
    col = "red"
  )
  
  graphics::abline(
    h = 0.80,
    col = "blue",
    lty = 2
  )
  
  # Plot 2: Mean connectivity
  graphics::plot(
    sft_filtered_median_mad$fitIndices[, 1],
    sft_filtered_median_mad$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = "Mean connectivity"
  )
  
  graphics::text(
    sft_filtered_median_mad$fitIndices[, 1],
    sft_filtered_median_mad$fitIndices[, 5],
    labels = sft_filtered_median_mad$fitIndices[, 1],
    col = "red"
  )
  
  grDevices::dev.off()
  
  
  
  net_cor_filtered <- blockwiseModules(
    expression_matrix,
    power            = opt_power,
    corType          = cor_method,
    networkType      = network_type,
    maxBlockSize     = 25000,             # single block if RAM permits
    minModuleSize    = 20,
    deepSplit        = 4,
    mergeCutHeight   = 0.1,
    reassignThreshold= 0.05,
    numericLabels    = TRUE,
    randomSeed       = 4321,
    saveTOMs = TRUE,
    saveTOMFileBase = "pwrTOM",
    verbose          = 3)
  mergedColors = labels2colors(net_cor_filtered$colors)
  names(mergedColors) = names(net_cor_filtered$colors)
  geneTree <- net_cor_filtered$dendrograms[[1]]
  blockIdx <- net_cor_filtered$blockGenes[[1]]
  cols4plot <- labels2colors(net_cor_filtered$colors)[ blockIdx ]
  
  dendrogram_plot_path <- file.path(
    plots_dir,
    "Gene_Dendrogram_Module_Colors.pdf"
  )
  
  grDevices::pdf(
    dendrogram_plot_path,
    height = 8,
    width = 14
  )
  
  WGCNA::plotDendroAndColors(geneTree,
                      cols4plot,
                      groupLabels = "Module colours",
                      dendroLabels = FALSE, hang = 0.03)
  
  grDevices::dev.off()
  
  return(list(
    module_colors = mergedColors,
    network       = net_cor_filtered
  ))
}