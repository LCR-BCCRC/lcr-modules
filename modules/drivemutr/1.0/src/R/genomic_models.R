calc_genomic_model <- function(df) {
  tryCatch({
    if (is.null(df) || nrow(df) == 0) return(tibble::tibble())
    
    df2 <- df %>% ungroup() %>% dplyr::select(-(1:5))
    
    stats::lm(residual ~ ., data = df2) %>% broom::tidy()
    
  }, error = function(e) tibble::tibble())
}

calc_significant_foci <- function(pred, df) {
  tryCatch({
    if (is.null(pred) || nrow(pred) == 0 || is.null(df) || nrow(df) == 0) {
      return(tibble::tibble())
    }
    
    sig_pred <- pred %>%
      dplyr::filter(term != "(Intercept)", p.value < 0.05) %>%
      dplyr::pull(term) %>%
      gsub("`", "", x = .)
    
    df2 <- df %>% ungroup() %>% dplyr::select(dplyr::any_of(c("residual", sig_pred)))
    
    stats::lm(residual ~ ., data = df2) %>% broom::tidy()
    
  }, error = function(e) tibble::tibble())
}

annotate_genomic_models <- function(target_regions_df) {

  # Compute Genomic_Model and Significant_Foci inside matched_mut_foci_data
  target_regions_df <- target_regions_df %>%
    dplyr::mutate(
      matched_mut_foci_data = purrr::map(
        matched_mut_foci_data,
        ~ {
          if (is.null(.x) || nrow(.x) == 0) return(.x)
          .x %>%
            dplyr::mutate(
              Genomic_Model    = purrr::map(Module_Model, calc_genomic_model),
              Significant_Foci = purrr::map2(Genomic_Model, Module_Model, calc_significant_foci)
            )
        }
      )
    )
  # Significant_Foci with final coefficients and genomic region coords
  annotate_foci_coefficients(target_regions_df)
}

annotate_foci_coefficients <- function(target_results) {
  target_results <- tryCatch({
    target_results %>%
      dplyr::mutate(
        matched_mut_foci_data = purrr::map2(
          ssm_matched,
          matched_mut_foci_data,
          ~ {
            ssm_matched_df <- .x
            nested_df <- .y
            
            if (is.null(ssm_matched_df) || is.null(nested_df) || nrow(nested_df) == 0) {
              return(nested_df)
            }
            
            nested_df %>%
              dplyr::mutate(
                Significant_Foci = purrr::pmap(
                  list(Genomic_Model, Significant_Foci, lambda_name),
                  function(Genomic_Model, Significant_Foci, lambda_name) {
                    tryCatch({
                      if (is.null(Genomic_Model) || nrow(Genomic_Model) == 0 ||
                          is.null(Significant_Foci) || nrow(Significant_Foci) == 0) {
                        return(tibble::tibble())
                      }
                      
                      genomic_model_tbl <- Genomic_Model %>%
                        dplyr::select(-dplyr::any_of("statistic")) %>%
                        dplyr::filter(term != "(Intercept)", term != "NA") %>%
                        dplyr::rename(
                          foci = term,
                          primary_coef = estimate,
                          primary_p.value = p.value,
                          primary_std.error = std.error
                        )
                      
                      sig_foci_tbl <- Significant_Foci %>%
                        dplyr::select(-dplyr::any_of("statistic")) %>%
                        dplyr::filter(term != "(Intercept)", term != "NA") %>%
                        dplyr::rename(
                          foci = term,
                          final_coef = estimate,
                          final_p.value = p.value,
                          final_std.error = std.error
                        )
                      
                      lambda_val_col <- paste0("group_", lambda_name)
                      
                      if (!lambda_val_col %in% names(ssm_matched_df)) {
                        return(tibble::tibble())
                      }
                      
                      dplyr::left_join(sig_foci_tbl, genomic_model_tbl, by = "foci") %>%
                        dplyr::mutate(
                          foci = gsub("`", "", foci),
                          ssm_filtered = purrr::map(
                            foci,
                            ~ ssm_matched_df %>%
                              dplyr::ungroup() %>%
                              dplyr::filter(.[[lambda_val_col]] == .x)
                          ),
                          region = purrr::map_chr(
                            ssm_filtered,
                            ~ {
                              if (is.null(.x) || nrow(.x) == 0) return(NA_character_)
                              chr <- unique(.x$Chromosome)[1]
                              paste0(
                                chr, ":",
                                min(.x$Start_Position, na.rm = TRUE),
                                "-",
                                max(.x$End_Position, na.rm = TRUE)
                              )
                            }
                          )
                        )
                      
                    }, error = function(e) {
                      message(conditionMessage(e))
                      tibble::tibble()
                    })
                  }
                )
              )
          }
        )
      )
  }, error = function(e) {
    message(conditionMessage(e))
    target_results
  })
}