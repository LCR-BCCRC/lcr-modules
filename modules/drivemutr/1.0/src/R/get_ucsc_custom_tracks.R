get_ucsc_custom_tracks <- function(target_results, outdir) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Color: blue = RNA-matched sample (ME present),
  # red = DNA-only sample. Score = genomic-model p-value for the foci.

  target_ucsc_points <- tryCatch({
    target_results %>%
      dplyr::mutate(
        matched_mut_foci_data = purrr::map2(
          ssm,
          matched_mut_foci_data,
          ~ {
            mut_df    <- .x
            nested_df <- .y

            if (is.null(mut_df) || is.null(nested_df) || nrow(nested_df) == 0) {
              return(nested_df)
            }

            nested_df %>%
              dplyr::mutate(
                ssm = purrr::pmap(
                  list(ME_data, Genomic_Model, lambda_name),
                  function(ME_data, Genomic_Model, lambda_name) {
                    tryCatch({
                      if (is.null(ME_data)       || nrow(ME_data) == 0 ||
                          is.null(Genomic_Model)  || nrow(Genomic_Model) == 0) {
                        return(tibble::tibble())
                      }

                      group_lambda_col <- paste0("group_", lambda_name)

                      ssm_joined <- mut_df %>%
                        dplyr::mutate(IDs_All = IDs) %>%
                        tidyr::unnest(IDs, keep_empty = TRUE) %>%
                        dplyr::left_join(
                          ME_data,
                          by = c("all_ids" = "Tumor_Sample_Barcode")
                        ) %>%
                        dplyr::group_by(Tumor_Sample_Barcode, Start_Position) %>%
                        dplyr::arrange(CADD_score) %>%
                        dplyr::slice_head(n = 1) %>%
                        dplyr::ungroup() %>%
                        dplyr::select(-id_type, -all_ids)

                      coef_df <- Genomic_Model %>%
                        dplyr::filter(term != "(Intercept)", term != "NA", term != "CN") %>%
                        dplyr::select(term, estimate, std.error, p.value) %>%
                        dplyr::mutate(
                          term = gsub("`", "", term),
                          term = as.character(term)
                        )

                      ssm_joined %>%
                        dplyr::mutate(foci_value = as.character(.[[group_lambda_col]])) %>%
                        dplyr::left_join(coef_df, by = c("foci_value" = "term")) %>%
                        dplyr::mutate(
                          p.value = dplyr::if_else(is.na(p.value), 1, p.value),
                          itemRgb = dplyr::if_else(
                            is.na(ME),
                            "255,0,0",   # red  = no RNA match
                            "0,0,255"    # blue = RNA-matched
                          )
                        ) %>%
                        dplyr::transmute(
                          Chromosome,
                          Start_Position,
                          End_Position,
                          name    = foci_value,
                          score   = p.value,
                          strand  = ".",
                          thickStart = Start_Position,
                          thickEnd   = End_Position,
                          itemRgb
                        )

                    }, error = function(e) tibble::tibble())
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

  ucsc_points_df <- target_ucsc_points %>%
    dplyr::select(Hugo_Symbol, matched_mut_foci_data) %>%
    tidyr::unnest(matched_mut_foci_data) %>%
    dplyr::select(Hugo_Symbol, lambda_name, ssm) %>%
    tidyr::unnest(ssm)

  ucsc_points_df %>%
    dplyr::group_by(lambda_name) %>%
    dplyr::group_walk(~ {
      lam <- .y$lambda_name[[1]]

      track_points <- paste0(
        'track name="Mutations_', lam,
        '" description="Mutations_', lam,
        '" visibility=2 itemRgb="On"'
      )

      out_df <- .x %>%
        dplyr::select(
          Chromosome, Start_Position, End_Position,
          name, score, strand, thickStart, thickEnd, itemRgb
        )

      readr::write_lines(
        c(track_points, readr::format_tsv(out_df, col_names = FALSE)),
        file.path(outdir, paste0("Mutation_Points_", lam, ".tsv"))
      )
    })


  # Color encodes genomic-model p-value significance.

  ucsc_blocks_df <- ucsc_points_df %>%
    dplyr::group_by(Hugo_Symbol, lambda_name, name) %>%
    dplyr::summarise(
      Chromosome     = dplyr::first(Chromosome),
      End_Position   = max(Start_Position, na.rm = TRUE) + 1L,
      Start_Position = min(Start_Position, na.rm = TRUE),
      score          = dplyr::first(score),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      score         = dplyr::if_else(is.na(score), 1, score),
      min_nonsig_gray = 100,
      max_gray        = 220,
      gamma           = 1.2,
      p_clamped       = pmin(pmax(score, 0.1), 1),
      t               = (p_clamped - 0.1) / 0.9,
      gray            = round(min_nonsig_gray + (max_gray - min_nonsig_gray) * (t ^ gamma)),
      itemRgb = dplyr::case_when(
        score < 0.05  ~ "237,113,23",                                # orange
        score <= 0.1  ~ "255,215,0",                                 # yellow
        score >  0.1  ~ sprintf("%d,%d,%d", gray, gray, gray)        # gray gradient
      ),
      strand     = ".",
      thickStart = Start_Position,
      thickEnd   = End_Position
    ) %>%
    dplyr::select(
      lambda_name, Chromosome, Start_Position, End_Position,
      name, score, strand, thickStart, thickEnd, itemRgb
    )

  ucsc_blocks_df %>%
    dplyr::group_by(lambda_name) %>%
    dplyr::group_walk(~ {
      lam <- .y$lambda_name[[1]]

      track_blocks <- paste0(
        'track name="Blocks_', lam,
        '" description="Blocks_', lam,
        '" visibility=2 itemRgb="On"'
      )

      out_df <- .x %>%
        dplyr::select(
          Chromosome, Start_Position, End_Position,
          name, score, strand, thickStart, thickEnd, itemRgb
        )

      readr::write_lines(
        c(track_blocks, readr::format_tsv(out_df, col_names = FALSE)),
        file.path(outdir, paste0("Mutation_Blocks_", lam, ".tsv"))
      )
    })


  # Transcription Factors 
  # One BED9 row per TF binding site disrupted by a significant-focus mutation.
  # Dark green colour for all rows.

  target_ucsc_tfs <- tryCatch({
    target_results %>%
      dplyr::mutate(
        matched_mut_foci_data = purrr::map(
          matched_mut_foci_data,
          ~ {
            matched_df <- .x

            if (is.null(matched_df) || nrow(matched_df) == 0) {
              return(matched_df)
            }

            matched_df %>%
              dplyr::mutate(
                Significant_Foci = purrr::map2(
                  Significant_Foci,
                  lambda_name,
                  ~ {
                    sig_foci   <- .x
                    lambda_val <- .y

                    if (is.null(sig_foci) || nrow(sig_foci) == 0) {
                      return(sig_foci)
                    }

                    sig_foci %>%
                      dplyr::mutate(
                        TF_Results = purrr::map(
                          TF_Results,
                          ~ {
                            tf_df <- .x

                            if (is.null(tf_df) || nrow(tf_df) == 0) {
                              return(tibble::tibble())
                            }

                            tf_df %>%
                              dplyr::mutate(
                                Chromosome     = seqnames,
                                Start_Position = start + purrr::map_dbl(motifPos, ~ sort(as.numeric(.x[1:2]))[1]),
                                End_Position   = end   + purrr::map_dbl(motifPos, ~ sort(as.numeric(.x[1:2]))[2]),
                                p.value        = 0,
                                itemRgb        = "9,79,41"   # dark green
                              ) %>%
                              dplyr::transmute(
                                Chromosome,
                                Start_Position,
                                End_Position,
                                name       = geneSymbol,
                                score      = p.value,
                                strand     = ".",
                                thickStart = Start_Position,
                                thickEnd   = End_Position,
                                itemRgb
                              )
                          }
                        )
                      )
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

  ucsc_tfs_df <- target_ucsc_tfs %>%
    dplyr::select(Hugo_Symbol, matched_mut_foci_data) %>%
    tidyr::unnest(matched_mut_foci_data) %>%
    dplyr::select(Hugo_Symbol, lambda_name, Significant_Foci) %>%
    tidyr::unnest(Significant_Foci) %>%
    dplyr::select(Hugo_Symbol, lambda_name, TF_Results) %>%
    tidyr::unnest(TF_Results)

  ucsc_tfs_df %>%
    dplyr::group_by(lambda_name) %>%
    dplyr::group_walk(~ {
      lam <- .y$lambda_name[[1]]

      track_tfs <- paste0(
        'track name="Transcription_Factors_', lam,
        '" description="Transcription_Factors_', lam,
        '" visibility=2 itemRgb="On"'
      )

      out_df <- .x %>%
        dplyr::select(
          Chromosome, Start_Position, End_Position,
          name, score, strand, thickStart, thickEnd, itemRgb
        )

      readr::write_lines(
        c(track_tfs, readr::format_tsv(out_df, col_names = FALSE)),
        file.path(outdir, paste0("Transcription_Factors_", lam, ".tsv"))
      )
    })

  invisible(NULL)
}