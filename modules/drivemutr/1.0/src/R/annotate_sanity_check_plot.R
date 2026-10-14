annotate_sanity_check_plot <- function(target_results) {
  target_results = tryCatch({
    target_results %>%
      dplyr::mutate(
        matched_mut_foci_data = purrr::map2(
          matched_mut_foci_data, Hugo_Symbol,
          function(matched_mut_foci_data, Hugo_Symbol){
            
            if (is.null(matched_mut_foci_data) || nrow(matched_mut_foci_data) == 0) return(matched_mut_foci_data)
            
            matched_mut_foci_data %>%
              dplyr::mutate(
                Sanity_Check_Plot = purrr::pmap(
                  list(Module_Model, Significant_Foci, lambda_name),
                  function(Module_Model, Significant_Foci, lambda_name){
                  tryCatch({
                    sig_terms <- Significant_Foci %>%
                      dplyr::filter(foci != "CN") %>%
                      dplyr::pull(foci)
                    
                    if (length(sig_terms) == 0) {
                      return(NULL)
                    }
                    
                    mut_cols <- names(Module_Model)[grepl("^\\d+(?:_\\d+)?$", names(Module_Model))]

                    sanity_df_non_mut = Module_Model %>%
                      filter(rowSums(across(all_of(mut_cols))) == 0) %>%
                      dplyr::mutate(group_name = "Not_Mutated") %>%
                      dplyr::select(-dplyr::all_of(mut_cols))
                    
                    sanity_df_mut = Module_Model %>%
                      filter(rowSums(across(all_of(mut_cols))) > 0) %>%
                      dplyr::mutate(group_name = "Mutated") %>%
                      dplyr::select(-dplyr::all_of(mut_cols))
                    
                    sanity_df_mut_groups <- Module_Model %>%
                      tidyr::pivot_longer(dplyr::all_of(mut_cols), names_to = "group_name", values_to = "mut_value") %>%
                      dplyr::filter(mut_value > 0) %>%
                      dplyr::select(-mut_value)
                    
                    sanity_df_curated <- dplyr::bind_rows(sanity_df_non_mut, sanity_df_mut, sanity_df_mut_groups)
                    
                    sanity_df_curated_sig = sanity_df_curated %>%
                      dplyr::filter(group_name %in% c("Not_Mutated", "Mutated", sig_terms))
                    
                    group_vals = unique(sanity_df_curated_sig$group_name)
                    
                    group_vals_num = setdiff(group_vals, c("Not_Mutated", "Mutated"))
                    
                    group_order <- group_vals_num[order(
                      as.numeric(sub("_.*$", "", group_vals_num)),
                      as.numeric(ifelse(grepl("_", group_vals_num), sub("^.*_", "", group_vals_num), 0))
                    )]
                    
                    levs = c("Mutated", "Not_Mutated", group_order)
                    
                    sanity_df_curated_sig$group_name <- factor(sanity_df_curated_sig$group_name, levels = levs, ordered = TRUE)
                    
                    gene_col <- Hugo_Symbol
                    
                    sanity_df_curated_sig$y_value <- sanity_df_curated_sig[[gene_col]]
                    
                    ggplot(sanity_df_curated_sig, aes(x = group_name, y = y_value)) +
                      geom_boxplot(outlier.shape = NA, width = 0.6) +
                      ggbeeswarm::geom_quasirandom(
                        aes(color = as.factor(CN)),
                        groupOnX = TRUE,
                        width = 0.15,
                        size = 2,
                        alpha = 0.9
                      ) +
                      labs(
                        y = "Normalized expression",
                        x = "Mutation status",
                        color = "CN",
                        title = paste0("Expression by foci group and mutation status: ", gene_col,
                                       " (lambda = ", gsub("_", ".", sub("^lambda_", "", lambda_name)), ")")
                      ) +
                      theme(
                        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)
                      )
                    
                  }, error = function(e) NULL)
                })
              )
          }
        )
      )
  }, error = function(e) {
    message(conditionMessage(e))
    target_results
  })          
}