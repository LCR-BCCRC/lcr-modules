extract_interaction_terms <- function(fit) {
  tryCatch({
    
    if (is.null(fit)) {
      return(tibble::tibble())
    }
    
    coef_df <- as.data.frame(coef(fit)) %>%
      tibble::as_tibble()
    
    # Case 1: coef(fit) already contains rule descriptions
    if (all(c("coefficient", "rule", "description") %in% colnames(coef_df))) {
      
      interaction_terms <- coef_df %>%
        dplyr::filter(
          coefficient != 0,
          !is.na(description),
          stringr::str_detect(description, stringr::fixed("&"))
        ) %>%
        dplyr::mutate(
          variables = stringr::str_extract_all(description, "X[0-9]+"),
          n_variables = purrr::map_int(variables, ~ length(unique(.x)))
        ) %>%
        dplyr::filter(n_variables >= 2)
      
      return(interaction_terms)
    }
    
    # Case 2: fallback — extract directly from fit$rules
    if (!is.null(fit[["rules"]])) {
      
      interaction_terms <- fit[["rules"]] %>%
        tibble::as_tibble() %>%
        dplyr::filter(
          !is.na(description),
          stringr::str_detect(description, stringr::fixed("&"))
        ) %>%
        dplyr::mutate(
          variables = stringr::str_extract_all(description, "X[0-9]+"),
          n_variables = purrr::map_int(variables, ~ length(unique(.x)))
        ) %>%
        dplyr::filter(n_variables >= 2)
      
      return(interaction_terms)
    }
    
    tibble::tibble()
    
  }, error = function(e) {
    tibble::tibble()
  })
}

get_rulefit_shapley_results <- function(target_results){
  target_results <- target_results %>%
    dplyr::mutate(
      matched_mut_foci_data = purrr::map(
        matched_mut_foci_data,
        function(matched_df) {
          
          if (is.null(matched_df) || nrow(matched_df) == 0) {
            return(matched_df)
          }
          
          matched_df %>%
            dplyr::mutate(
              Interaction_Results = purrr::map2(
                Significant_Foci,
                Module_Model,
                function(Significant_Foci, Module_Model) {
                  
                  foci_cols <- Significant_Foci %>%
                    dplyr::filter(foci != "CN") %>%
                    dplyr::pull(foci)
                  
                  rulefit_df <- Module_Model %>%
                    dplyr::select(
                      1,
                      5,
                      dplyr::any_of(foci_cols),
                      residual
                    )
                  
                  # Rename columns whose names are only numbers: "1" -> "X1"
                  numeric_col_names <- names(rulefit_df)[
                    grepl("^[0-9]+$", names(rulefit_df))
                  ]
                  
                  names(rulefit_df)[names(rulefit_df) %in% numeric_col_names] <-
                    paste0("X", numeric_col_names)
                  
                  # Remove the first two columns before fitting RuleFit
                  # but keep residual as the response
                  rulefit_model_data <- rulefit_df %>% ungroup() %>%
                    dplyr::select(-c(1, 2))
                  
                  rulefit_model <- tryCatch({
                    set.seed(5678)
                    pre::pre(
                      residual ~ .,
                      data = rulefit_model_data,
                      family = "gaussian",
                      maxdepth = 2L,
                      ntrees = 1000L
                    )
                  }, error = function(e) {
                    NULL
                  })
                  
                  interaction_terms <- extract_interaction_terms(rulefit_model)
                  
                  phi_mean <- tibble::tibble(
                    feature = character(),
                    phi_mean = double()
                  )
                  shap_mean_plot <- NULL
                  shap_beeswarm_plot <- NULL
                  
                  if (nrow(interaction_terms) > 0) {
                    
                    vars <- unique(unlist(interaction_terms$variables))
                    rows_ids = which(rowSums(rulefit_model_data[vars] != 0, na.rm = TRUE) == length(vars))
                    
                    X <- rulefit_model_data[, setdiff(names(rulefit_model_data),"residual"), drop = FALSE]
                    y <- rulefit_model_data$residual
                    
                    predictor_object <- tryCatch({
                      iml::Predictor$new(
                        model = rulefit_model,
                        data = X,
                        y = y,
                        predict.function = function(model, newdata) {
                          as.numeric(stats::predict(model, newdata))
                        }
                      )
                    }, error = function(e) {
                      NULL
                    })
                    
                    if (!is.null(predictor_object) && length(rows_ids) > 0) {
                      
                      phi_res = lapply(rows_ids, function(i){
                        sh <- Shapley$new(predictor_object, x.interest = X[i, , drop = FALSE])
                        res <- as.data.frame(sh$results) %>% filter(phi != 0)
                        data.frame(
                          row_index = i,
                          feature = res$feature,
                          phi = res$phi,
                          feature.value = res$feature.value
                        )
                      })
                      phi_res_df <- do.call(rbind, phi_res)
                      
                      phi_res_df <- phi_res_df %>%
                        mutate(
                          feature_value = as.numeric(sub(".*=", "", feature.value))
                        )
                      
                      phi_mean = phi_res_df %>% group_by(feature) %>% summarise(phi_mean = mean(phi))
                      
                      shap_mean_plot <- ggplot2::ggplot(
                        phi_mean,
                        ggplot2::aes(x = .data$phi_mean, y = .data$feature)
                      ) +
                        ggplot2::geom_col(
                          width = 0.75,
                          fill = "grey40"
                        ) +
                        ggplot2::geom_vline(
                          xintercept = 0,
                          linetype = "dashed",
                          linewidth = 0.7
                        ) +
                        ggplot2::theme_minimal(base_size = 14) +
                        ggplot2::theme(
                          panel.grid.major.y = ggplot2::element_blank(),
                          panel.border = ggplot2::element_rect(
                            color = "black",
                            fill = NA,
                            linewidth = 0.8
                          )
                        ) +
                        ggplot2::labs(
                          x = "Mean Shapley value (phi)",
                          y = NULL
                        )
                      
                      shap_beeswarm_plot <- ggplot2::ggplot(
                        phi_res_df,
                        ggplot2::aes(
                          x = .data$phi,
                          y = .data$feature,
                          color = .data$feature_value
                        )
                      ) +
                        ggplot2::geom_vline(
                          xintercept = 0,
                          linetype = "dashed",
                          linewidth = 0.6
                        ) +
                        ggbeeswarm::geom_quasirandom(
                          width = 0.25,
                          alpha = 0.9,
                          size = 2
                        ) +
                        ggplot2::scale_color_gradient(
                          low = "blue",
                          high = "red",
                          na.value = "grey70"
                        ) +
                        ggplot2::labs(
                          x = "SHAP value (impact on model output)",
                          y = NULL,
                          color = "Feature value"
                        ) +
                        ggplot2::theme_minimal(base_size = 13) +
                        ggplot2::theme(
                          panel.grid.major.y = ggplot2::element_blank()
                        )
                      
                    }
                  }
                  
                  tibble::tibble(
                    RuleFit_data = list(rulefit_df),
                    RuleFit_result = list(rulefit_model),
                    Interaction_terms = list(interaction_terms),
                    Phi_mean = list(phi_mean),
                    SHAP_Mean_Plot = list(shap_mean_plot),
                    SHAP_Beeswarm_Plot = list(shap_beeswarm_plot)
                  )
                }
              )
            )
        }
      )
    )
}
