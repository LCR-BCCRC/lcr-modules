annotate_ids <- function(ssm_df, sample_id_aliases = NULL) {
    tryCatch(
        ssm_df %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                IDs = list({
                if (is.na(Tumor_Sample_Barcode)) {
                    tibble::tibble(id_type = NA_character_, all_ids = NA_character_)
                } else if (is.null(sample_id_aliases)) {
                    tibble::tibble(
                    id_type = c("mrna_sample_id", "genome_sample_id", "capture_sample_id"),
                    all_ids = Tumor_Sample_Barcode
                    )
                } else {
                    id_match <- which(sample_id_aliases$all_ids %in% Tumor_Sample_Barcode)
                    if (length(id_match) > 0) {
                    sid <- sample_id_aliases$sample_id[id_match[1]]
                    sample_id_aliases %>%
                    dplyr::filter(sample_id == sid) %>%
                        dplyr::select(id_type, all_ids)
                    } else {
                        tibble::tibble(id_type = NA_character_, all_ids = Tumor_Sample_Barcode)
                    }
                }
                })
            ) %>%
            dplyr::ungroup()
            , error = function(e) {
            message("annotate_ids failed: ", conditionMessage(e))
            ssm_df
            }
        )
}
