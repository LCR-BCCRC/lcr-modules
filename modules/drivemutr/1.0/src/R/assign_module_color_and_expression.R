# Arguments:
# target_regions_df      one-row-per-gene tibble (cadd_adjusted.rds)
#   mergedColors           named character vector produced by WGCNA:
#                            names  = gene names (== colnames of expression matrix)
#                            values = colour strings ("turquoise", "blue", …)
#   high_var_expression_matrix
#                          samples × genes expression matrix (filtered_expression.rds)
#   sample_id_aliases      long-format linkage data frame with columns:
#                            sample_id, id_type, all_ids
#                          Pass NULL (default) when all modalities share the same IDs.

assign_module_color_and_expression <- function(target_regions_df,
                                                mergedColors,
                                                high_var_expression_matrix,
                                                sample_id_aliases = NULL) {


  target_gene_regions_data <- target_regions_df %>%
    mutate(
      Module_Color = mergedColors[
        match(Hugo_Symbol, colnames(high_var_expression_matrix))
      ]
    ) %>%
    mutate(
      Module = purrr::map(Module_Color, ~
        tryCatch(
          as.data.frame(high_var_expression_matrix) %>%
            dplyr::select(all_of(
              colnames(high_var_expression_matrix)[mergedColors == .x]
            )) %>%
            rownames_to_column("Tumor_Sample_Barcode") %>%
            annotate_ids(sample_id_aliases = sample_id_aliases),
          error = function(e) NULL
        )
      )
    )
  
  return(target_gene_regions_data)
}