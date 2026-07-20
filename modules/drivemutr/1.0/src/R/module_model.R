calc_ME = function(genomic_matrix, Module, Hugo_Symbol){
  tryCatch({
    if (is.null(genomic_matrix) || is.null(Module)) return(NULL)
    qrW <- genomic_matrix %>%
      dplyr::select(CN, matches("^\\d+$")) %>%
      # dplyr::rename_with(~ paste0("G", .x), matches("^\\d+$")) %>%
      model.matrix(~ ., data = .) %>%
      qr()
    module_mx <- Module %>%
      tibble::column_to_rownames("Tumor_Sample_Barcode") %>%
      dplyr::select(-dplyr::all_of(c("CN", "IDs_All", Hugo_Symbol))) %>%
      as.matrix()
    R    <- qr.resid(qrW, module_mx)
    mu_r <- colMeans(R, na.rm = TRUE)
    sd_r <- apply(R, 2, sd, na.rm = TRUE)
    
    Z    <- sweep(sweep(R, 2, mu_r, "-"), 2, sd_r, "/")
    
    pca <- prcomp(Z, center = FALSE, scale. = FALSE)
    v1  <- pca$rotation[, 1]               # PC1 loadings (genes)
    scores <- as.numeric(Z %*% v1)         # ME_purged (length n)
    
    cors <- cor(scores, Z, use = "pairwise.complete.obs")
    if (mean(cors, na.rm = TRUE) < 0) { v1 <- -v1; scores <- -scores }
    names(scores) <- Module$Tumor_Sample_Barcode
    Module %>%
      dplyr::mutate(
        ME = unname(scores[as.character(Tumor_Sample_Barcode)])
      ) %>%
      left_join(genomic_matrix %>% dplyr::select(Tumor_Sample_Barcode, Mutated), by = "Tumor_Sample_Barcode") %>%
      dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol , ME , Mutated)
  }, error = function(e) tibble()
  )
}




calc_module_model <- function(genomic_matrix, ME_data, Hugo_Symbol){
  tryCatch({
    module_model_train = ME_data %>% filter(Mutated == 0) %>%
      dplyr::select(-Mutated) %>%
      column_to_rownames(var = "Tumor_Sample_Barcode")
    
    module_model_test = ME_data %>% filter(Mutated == 1) %>%
      dplyr::select(-Mutated) %>%
      column_to_rownames(var = "Tumor_Sample_Barcode")
    
    
    module_model = lm(module_model_train[[Hugo_Symbol]] ~ ME , data = module_model_train)
    module_model_nonmut_resid = augment(module_model) %>%
      dplyr::select(.rownames ,.resid)
    
    module_model_mut_resid = data.frame(
      .rownames = rownames(module_model_test),
      .resid = module_model_test[[Hugo_Symbol]] - predict(module_model, newdata = module_model_test)
    )
    
    module_model_all_resid = bind_rows(module_model_nonmut_resid, module_model_mut_resid)
    colnames(module_model_all_resid) = c("Tumor_Sample_Barcode", "residual")
    genomic_matrix %>% left_join(module_model_all_resid, by = "Tumor_Sample_Barcode")
  },error = function(e) tibble()
  )
}





annotate_module_model_results <- function(target_regions_df){
  
  target_gene_regions_data <- target_regions_df %>%
    mutate(
      matched_mut_foci_data = purrr::pmap(
        list(matched_mut_foci_data, Module, Hugo_Symbol),
        function(matched_mut_foci_data, Module, Hugo_Symbol) {
          
          if (is.null(matched_mut_foci_data) || nrow(matched_mut_foci_data) == 0) {
            return(matched_mut_foci_data)
          }
          
          matched_mut_foci_data %>%
            mutate(
              ME_data = purrr::map(
                genomic_matrix, 
                ~ calc_ME(.x, Module, Hugo_Symbol)
              ),
              Module_Model = purrr::map2(
                genomic_matrix,
                ME_data,
                ~ calc_module_model(.x, .y, Hugo_Symbol)
              ),
              Module_Model = purrr::map2(
                ME_data,
                Module_Model,
                ~ tryCatch(
                  .x %>%
                    dplyr::left_join(
                      .y %>% dplyr::select(-Mutated),
                      by = "Tumor_Sample_Barcode"
                    ),
                  error = function(e) tibble::tibble()
                )
              )
            )
        }
      )
    )
  
  return(target_gene_regions_data)
}
