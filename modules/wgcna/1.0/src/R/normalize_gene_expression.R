normalize_gene_expression_by_pathology <- function(
      raw_gene_exp_path,
      failed_qc_path,
      samples_metadata_path,
      pathology_arg = NULL,
      plots_dir,
      batches    = c("protocol", "ffpe_or_frozen"),
      cohort_var = "cohort",
      bio_var    = "pathology",
      min_reads  = 3e6,
      threads   = 1
  ) {
  
  # Register parallel backend — MulticoreParam uses forking (Linux/macOS only).
  # Falls back to SnowParam on Windows.
  backend <- if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(threads)
  } else {
    BiocParallel::MulticoreParam(threads)
  }
  BiocParallel::register(backend)
  message(paste0("Parallelisation: ", threads, " thread(s) via ", class(backend)[1]))

  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_raw_gene_exp = read_tsv(raw_gene_exp_path)
  
  samples_metadata = read_tsv(samples_metadata_path)
  
  if (!is.null(pathology_arg)) {
    pathology_metadata <- samples_metadata %>% filter(pathology == pathology_arg)
  } else {
    message("No pathology filter specified — using all samples.")
    pathology_metadata <- samples_metadata
  } 
  
  pathology_sample_ids = pathology_metadata %>% pull(sample_id)

  # Drop failed QC samples
  if (file.exists(failed_qc_path)) {
    dropped_rna_seq <- readr::read_tsv(
      failed_qc_path,
      show_col_types = FALSE
    )
    
    dropped_rna_seq_sample_id <- dropped_rna_seq %>%
      dplyr::pull(sample_id)
    
    pathology_sample_ids <- setdiff(
      pathology_sample_ids,
      dropped_rna_seq_sample_id
    )
    
    pathology_metadata = pathology_metadata %>% filter(sample_id %in% pathology_sample_ids)
  }
  
  # Get raw expression data for the selected samples
  all_pathology_raw_exp <- all_raw_gene_exp %>%
    dplyr::select(gene, any_of(pathology_sample_ids))
  
  # Set gene names as rownames to make it as a matrix
  all_pathology_raw_exp <- all_pathology_raw_exp %>% column_to_rownames(var = "gene")
  
  
  ### DESeq2 ###
  
  # It takes a raw counts matrix and outputs a variance-
  # stabilizing transformed matrix. Batch effect correction
  # is performed based on specified variables in the input
  # metadata.
  
  matrix <- all_pathology_raw_exp
  meta = pathology_metadata %>% column_to_rownames(var = "sample_id")
  
  
  # Filter the matrix to samples meeting the minimum read threshold
  matrix_filt <- matrix[, colSums(matrix) > min_reads]
  
  
  dropped <- colnames(matrix)[!colnames(matrix) %in% colnames(matrix_filt)]
  
  dropped <- colSums(matrix[, dropped, drop = FALSE])
  
  if (length(dropped) > 0) {
    print(paste0(
      "The following samples have total reads below the threshold of ",
      min_reads,
      ":\n"
    ))
    print(data.frame(total_reads = dropped))
  }
  
  
  
  meta <- meta[colnames(matrix_filt), ]
  
  meta_samples <- length(rownames(meta))
  matrix_samples <- length(colnames(matrix_filt))
  message(paste0(
    "There are ",
    matrix_samples,
    " samples in the matrix and ",
    meta_samples,
    " samples in the metadata."
  ))
  if (!identical(rownames(meta), colnames(matrix_filt))) {
    stop("The sample_ids in the matrix and metadata are not identical")
  }
  
  # Remove any batches with only one factor level
  final_batches <- c()
  
  for (i in batches) {
    if (length(unique(meta[[i]])) <= 1) {
      print(paste0(
        "The batch '", i, "' has only one level: '", unique(meta[[i]]), "'. Dropping."
      ))
    } else {
      final_batches <- c(final_batches, i)
      meta[[i]] <- factor(meta[[i]])
    }
  }
  
  print(paste0("Final batches for batch effect correction: ", paste0(final_batches, collapse = ", ")))
  
  # Run DESeq2
  
  message("Running DESeq2...")
  matrix_filt <- as.matrix(matrix_filt)
  
  dds <- DESeqDataSetFromMatrix(matrix_filt,
                                colData = meta,
                                design = ~1
  )
  
  
  dds <- DESeq(dds,
               parallel = TRUE
  )
  
  # Perform variance stabilzing transformation
  
  message("Performing variance stabilizing transformation...")
  vsd <- vst(dds, blind = FALSE)
  
  
  # Plot PCA for each uncorrected batch
  for (batch in c(final_batches, bio_var, cohort_var)) {
    plot_path <- file.path(plots_dir, paste0("PCA_uncorrected_", batch, ".pdf"))
    pdf(plot_path, height = 10, width = 10)
    print(plotPCA(vsd, intgroup = batch))
    dev.off()
  }
  
  
  mat <- assay(vsd)
  
  
  # Perform batch effect correction on each batch
  
  if (length(final_batches) >= 1) {
    message("Performing batch effect correction...")
    vsd_corr <- vsd
    
    batch_matrix <- model.matrix(
      as.formula(paste("~", paste(final_batches, collapse = " + "))),
      meta
    )
    print(paste0("there are ", nrow(batch_matrix), " rows in the batch matrix"))
    print(paste0("there are ", ncol(mat), " columns in the expression matrix"))
    
    bio_vars_present <- bio_var[bio_var %in% colnames(meta)]
    bio_vars_ok      <- bio_vars_present[
      sapply(bio_vars_present, function(v) length(unique(meta[[v]])) > 1)
    ]
    
    if (length(bio_vars_ok) > 0) {
      design_matrix <- model.matrix(
        as.formula(paste("~", paste(bio_vars_ok, collapse = " + "))),
        meta
      )
      message(paste0(
        "Using design matrix for biological variable(s): ",
        paste(bio_vars_ok, collapse = ", ")
      ))
      mat <- removeBatchEffect(mat, covariates = batch_matrix, design = design_matrix)
    } else {
      mat <- removeBatchEffect(mat, covariates = batch_matrix)
    }
    
    
    assay(vsd_corr) <- mat
    
    # Plot PCA for each corrected batch plus cohort
    for (batch in c(final_batches, cohort_var, bio_var)) {
      plot_path <- file.path(plots_dir, paste0("PCA_corrected_", batch, ".pdf"))
      pdf(plot_path, height = 10, width = 10)
      print(plotPCA(vsd_corr, intgroup = batch))
      dev.off()
    }
  }
  
  mat <- data.frame(mat, check.names = FALSE)
  
  
  mat_renamed <- mat
  mat_renamed$gene_id <- rownames(mat)
  
  mat_renamed = GAMBLR.helpers::translate_ensg_ids(df = mat_renamed, ensg_column_in = "gene_id")
  mat_renamed = mat_renamed %>% tidyr::drop_na(hgnc_symbol)
  
  message("Complete.")
  
  genes_dropped = c("AKAP17A", "AMD1P2", "ASMT", "ASMTL", "ASMTL-AS1", "CD99", "CD99P1", "CRLF2", "CSF2RA", "DDX11L16", "DHRSX", "DHRSX-IT1", "DPH3P2", "DUXAP8", "ELOCP24", "FABP5P13", "GOLGA8M", "GPR84-AS1", "GTPBP6", "HERC3", "IL3RA", "IL9R", "ITFG2-AS1", "KRT18P53", "LINC00102", "LINC00106", "LINC00685", "LINC01238", "LINC02203", "LINC02968", "MIR3690", "MIR6089", "P2RY8", "PINX1", "PLCXD1", "POLR2J3", "PPP2R3B", "RNA5SP498", "RPL14P5", "SCARNA4", "SHOX", "SIGLEC5", "SLC25A6", "SNORA16A", "SNORA17B", "SNORA50A", "SNORD38B", "SNORD3D", "SPRY3", "TBCE", "TNFRSF10A-DT", "TRPC6P1", "VAMP7", "WASH6P", "WASIR1", "ZBED1")
  mat_renamed <- mat_renamed %>% 
    filter(!(hgnc_symbol %in% genes_dropped))
  
  return(mat_renamed)
}

