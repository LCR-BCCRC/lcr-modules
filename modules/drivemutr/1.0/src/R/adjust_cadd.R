adjust_cadd = function(ssm_df, lambda = 1){
  
  if (length(lambda) != 1) {
    stop("Expects a single lambda value.")
  }
  
  base_group_col <- paste0("group_", lambda_suffix(lambda))
  custom_group_col <- paste0(base_group_col, "_custom")
  
  
  has_base <- base_group_col %in% names(ssm_df)
  has_custom <- custom_group_col %in% names(ssm_df)
  
  if (!has_base && !has_custom) {
    return(ssm_df)
  }
  
  tryCatch({
    
    out <- ssm_df
    
    if (has_base) {
      out <- out %>%
        dplyr::group_by(Tumor_Sample_Barcode, !!rlang::sym(base_group_col)) %>%
        dplyr::mutate(max_CADD = max(CADD_score, na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
    
    if (has_custom) {
      out <- out %>%
        dplyr::group_by(Tumor_Sample_Barcode, !!rlang::sym(custom_group_col)) %>%
        dplyr::mutate(max_CADD_custom = max(CADD_score, na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
    
    out
    
  }, error = function(e) tibble::tibble())
}