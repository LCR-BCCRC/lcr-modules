get_ids <- function(ssm_df, sample_col = "Tumor_Sample_Barcode") {
  if (is.null(ssm_df) || nrow(ssm_df) == 0) return(character())
  unique(na.omit(as.character(ssm_df[[sample_col]])))
}

filter_matched_DNA_RNA_samples <- function(target_regions_df){
  
  # filter to those samples we have matched DNA and RNA seq data
  target_gene_regions_data <- target_regions_df %>%
    mutate(
      Module = map2(Module, CN, ~ {
        tryCatch(
          .x %>%
            mutate(IDs_All = IDs) %>%
            unnest(IDs, keep_empty = TRUE) %>%
            left_join(.y, by = c("all_ids" = "sample_id")) %>%
            group_by(Tumor_Sample_Barcode) %>%
            arrange((CN)) %>%
            slice_head() %>%
            dplyr::select(-id_type, -all_ids) %>%
            filter(!is.na(CN)) %>%
            mutate(
              CN = case_when(CN < 2 ~ 1, 
                             CN == 2 ~ 2,
                             CN > 2 ~ 3)
            )
          ,
          error = function(e) NULL
        )
      } 
      ),
      ssm_matched = map(Module, ~ {
        if (is.null(.x)) return(NULL)
        .x %>%
          dplyr::select(Tumor_Sample_Barcode, IDs_All, CN) 
      }),
      ssm_matched = map2(ssm_matched, ssm, ~ {
        mutated_ids = get_ids(.y)
        tryCatch(
          .x %>%
            mutate(IDs = IDs_All) %>%                 # keep the original list-column
            tidyr::unnest(IDs_All, keep_empty = TRUE) %>%  # creates id_type, all_ids
            full_join(.y %>% dplyr::select(-IDs), by = c("all_ids" = "Tumor_Sample_Barcode"), relationship = "many-to-many") %>%
            group_by(Tumor_Sample_Barcode, Start_Position) %>%
            arrange((CADD_score)) %>%
            slice_head() %>%
            dplyr::select(-id_type, -all_ids) %>% 
            filter(!is.na(Tumor_Sample_Barcode)) %>%
            mutate(
              Mutated = map_int(IDs, ~{
                
                any_mut  <- any(.x$all_ids %in% mutated_ids,  na.rm = TRUE)
                
                dplyr::case_when(
                  any_mut ~ 1L,
                  TRUE    ~ 0L
                )
              })
            )
          ,
          error = function(e) tibble()
        )
      } 
      )
    )
  
  return(target_gene_regions_data)
  
}