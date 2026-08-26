# Arguments:
#   target_regions_df   one-row-per-gene tibble (cadd_adjusted.rds)
#   mergedColors        named character vector produced by WGCNA:
#                         names  = gene names (the WGCNA gene set)
#                         values = colour strings ("turquoise", "blue", ...)
#   expression_matrix   samples x genes expression matrix. May contain MORE genes
#                         than mergedColors (e.g. the full normalized matrix);
#                         module membership is resolved by gene NAME below, so a
#                         superset matrix is fine and the two need not be aligned.
#   sample_id_aliases   long-format linkage data frame with columns:
#                         sample_id, id_type, all_ids
#                         Pass NULL (default) when all modalities share the same IDs.

assign_module_color_and_expression <- function(target_regions_df,
                                                mergedColors,
                                                expression_matrix,
                                                sample_id_aliases = NULL) {

  target_gene_regions_data <- target_regions_df %>%
    mutate(
      # Look up the module colour of each target gene by name (NA if the gene is
      # not part of any WGCNA module, e.g. it was dropped by the variance filter).
      Module_Color = unname(mergedColors[Hugo_Symbol])
    ) %>%
    mutate(
      Module = purrr::map(Module_Color, ~
        tryCatch({
          if (is.na(.x)) return(NULL)
          # Genes belonging to this module (by name), restricted to those present
          # in the expression matrix.
          module_genes <- names(mergedColors)[mergedColors == .x]
          module_genes <- intersect(module_genes, colnames(expression_matrix))
          as.data.frame(expression_matrix) %>%
            dplyr::select(all_of(module_genes)) %>%
            rownames_to_column("Tumor_Sample_Barcode") %>%
            annotate_ids(sample_id_aliases = sample_id_aliases)
        }, error = function(e) NULL)
      )
    )

  return(target_gene_regions_data)
}