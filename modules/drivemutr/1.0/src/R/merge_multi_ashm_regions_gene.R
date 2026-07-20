merge_multi_ashm_regions_gene <- function(target_regions_df){
  
  if ("CN" %in% colnames(target_regions_df)) {
    
    target_regions_df <- target_regions_df %>%
      dplyr::group_by(Hugo_Symbol) %>%
      dplyr::summarise(
        Pathology = dplyr::first(Pathology),
        NCBI_Build = dplyr::first(NCBI_Build),
        Gene_Region = list(
          tibble::tibble(
            Gene_Region = Gene_Region,
            Region = Region
          )
        ),
        CN = CN[1],
        ssm = list(
          dplyr::bind_rows(ssm) %>%
            dplyr::distinct()
        ),
        .groups = "drop"
      )
    
  } else {
    
    target_regions_df <- target_regions_df %>%
      dplyr::group_by(Hugo_Symbol) %>%
      dplyr::summarise(
        Pathology = dplyr::first(Pathology),
        NCBI_Build = dplyr::first(NCBI_Build),
        Gene_Region = list(
          tibble::tibble(
            Gene_Region = Gene_Region,
            Region = Region
          )
        ),
        ssm = list(
          dplyr::bind_rows(ssm) %>%
            dplyr::distinct()
        ),
        .groups = "drop"
      )
  }
  
  return(target_regions_df)
}