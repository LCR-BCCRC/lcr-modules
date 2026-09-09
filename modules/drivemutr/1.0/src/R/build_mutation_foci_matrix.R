build_mutation_foci_matrix = function(ssm_df, lambda = 1){
  
  if (length(lambda) != 1) {
    stop("Expects a single lambda value.")
  }
  
  suffix <- lambda_suffix(lambda)
  
  base_group_col <- paste0("group_", suffix)
  custom_group_col <- paste0(base_group_col, "_custom")
  
  base_value_col <- paste0("max_CADD_", suffix)
  custom_value_col <- paste0("max_CADD_", suffix, "_custom")
  
  has_base <- base_group_col %in% names(ssm_df)
  has_custom <- custom_group_col %in% names(ssm_df)
  
  out <- list()
  if (has_base) {
    out[[suffix]] <- tryCatch({
      ssm_df %>%
        group_by(!!rlang::sym(base_group_col)) %>%
        filter(n()>5) %>%
        ungroup() %>%
        dplyr::select(Tumor_Sample_Barcode, IDs, CN,
                      !!rlang::sym(base_group_col),
                      !!rlang::sym(base_value_col),
                      Mutated) %>%
        unique() %>%
        pivot_wider(names_from = !!rlang::sym(base_group_col),
                    values_from = !!rlang::sym(base_value_col)) %>%
        dplyr::select(-any_of("NA")) %>%
        mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))
    },  error = function(e) tibble::tibble())
  }
  
  if (has_custom) {
    out[[paste0(suffix, "_custom")]] <- tryCatch({
      ssm_df %>%
        group_by(!!rlang::sym(custom_group_col)) %>%
        filter(n()>5) %>%
        ungroup() %>%
        dplyr::select(Tumor_Sample_Barcode, IDs, CN,
                      !!rlang::sym(custom_group_col),
                      !!rlang::sym(custom_value_col),
                      Mutated) %>%
        unique() %>%
        pivot_wider(names_from = !!rlang::sym(custom_group_col),
                    values_from = !!rlang::sym(custom_value_col)) %>%
        dplyr::select(-any_of("NA")) %>%
        mutate(across(where(is.numeric), ~ replace(., is.na(.), 0)))
    }, error = function(e) tibble::tibble())
  }
  
  return(out)
}