get_nucleotides_prop = function(region){
  m <- stringr::str_match(region, "^(?:chr)?([^:]+):(\\d+)-(\\d+)$")
  if (anyNA(m)) return(tibble(A=NA_real_, C=NA_real_, G=NA_real_, T=NA_real_))
  
  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", m[1,2]),
    ranges   = IRanges::IRanges(start = as.integer(m[1,3]), end = as.integer(m[1,4]))
  )
  
  seq    <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, gr)
  counts <- Biostrings::letterFrequency(seq, letters = c("A","C","G","T"))
  pct    <- as.numeric(counts / sum(counts))
  
  tibble(A=pct[1], C=pct[2], G=pct[3], T=pct[4])
}


annotate_nucleotides_prop = function(target_results){
  target_results <- target_results %>%
    dplyr::mutate(
      Gene_Region = purrr::map(
        Gene_Region,
        ~ {
          .x %>%
            dplyr::mutate(
              ATCG_prop = purrr::map(Region, get_nucleotides_prop)
            ) %>%
            tidyr::unnest_wider(ATCG_prop)
        }
      )
    ) %>%
    dplyr::mutate(
      aSHM_avg_ATCG_prop = purrr::map(
        Gene_Region,
        ~ {
          .x %>%
            dplyr::summarise(
              A = mean(A, na.rm = TRUE),
              C = mean(C, na.rm = TRUE),
              G = mean(G, na.rm = TRUE),
              T = mean(T, na.rm = TRUE)
            )
        }
      )
    )
}



bed_maker = function(maf){
  maf_to_bed_fixed = maf %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      chrom = ifelse(grepl("^chr", Chromosome), Chromosome, paste0("chr", Chromosome)),
      
      Reference_Allele = dplyr::case_when(
        Variant_Type == "DEL" ~ paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, GRanges(chrom, IRanges(start = Start_Position - 1, end = Start_Position - 1)))), Reference_Allele),
        TRUE ~ as.character(Reference_Allele)
      ),
      Tumor_Seq_Allele2 = dplyr::case_when(
        Variant_Type == "DEL" ~ as.character(getSeq(
          BSgenome.Hsapiens.UCSC.hg19,
          GRanges(chrom, IRanges(start = Start_Position - 1, end = Start_Position - 1)))),
        Variant_Type == "INS" ~ paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, GRanges(chrom, IRanges(start = Start_Position, end = Start_Position)))), Tumor_Seq_Allele2, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, GRanges(chrom, IRanges(start = End_Position, end = End_Position))))),
        TRUE ~ as.character(Tumor_Seq_Allele2)
      ),
      
      chromStart = dplyr::if_else(Variant_Type == "DEL", Start_Position - 2L, Start_Position - 1L),
      chromEnd   = End_Position,
      name       = paste0(chrom, ":", chromStart, ":", Reference_Allele, ":", Tumor_Seq_Allele2)
    ) %>%
    transmute(chrom, chromStart, chromEnd, name)
  return(maf_to_bed_fixed)
}


# Don't need to run seprately as it is part of get_motifbreakr()
bed_file_snps = function(bed){
  bed = bed %>% unique()
  
  #tf_temp_file <- tempfile(fileext = ".bed")
  tf_temp_file <- file.path(getwd(), "my_tf_temp_file.bed")
  
  write.table(bed, file = tf_temp_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  snps_tf <- motifbreakR::snps.from.file(
    file = tf_temp_file,
    indels = TRUE,
    search.genome = BSgenome.Hsapiens.UCSC.hg19,
    format = "bed"
  )
  return(snps_tf)
}


MotifDB_human <- MotifDb::query(MotifDb::MotifDb, "Hsapiens")

get_motifbreakr = function(bed, A, C, G, T){
  
  snps_tf = bed_file_snps(bed)
  
  motif_result_temp <- motifbreakR::motifbreakR(
    snps_tf,
    filterp   = TRUE,
    pwmList   = MotifDB_human,
    threshold = 1e-4,
    method    = "ic",
    bkg       = c(A=A, C=C, G=G, T=T),
    BPPARAM   = BiocParallel::bpparam()
  )
  
  #TF_Results = as.data.frame(motif_result_temp, row.names = NULL)
  TF_Results = motif_result_temp
  return(TF_Results)
}



### Do not use ssm_matched here as we want to have all the mutations there (ssm needed)
fisher_aid_p = function(num_tf_mut_target, num_mut_target, maf_out_foci, motifbr){
  
  start_mut = as.character(maf_out_foci$Start_Position)
  
  num_tf_mut_out_foci = motifbr %>%
    dplyr::mutate(start = as.character(start)) %>%
    dplyr::group_by(geneSymbol) %>%
    dplyr::summarise(n_mut = sum(start_mut %in% start), # We do it as in bed we just have the unique exp (no rep)
                     .groups = "drop"
    )
  
  num_target_tf_mut_out = num_tf_mut_out_foci %>%
    filter(geneSymbol == num_tf_mut_target$TF_name) %>%
    pull(n_mut)
  num_target_tf_mut_out <- if (length(num_target_tf_mut_out) == 0) 0 else num_target_tf_mut_out 
  
  
  p_fisher <- fisher.test(matrix(c(num_tf_mut_target$n_mut,
                                   (num_mut_target - num_tf_mut_target$n_mut),
                                   num_target_tf_mut_out,
                                   (nrow(maf_out_foci)- num_target_tf_mut_out)
  ), nrow=2),alternative="greater")$p.value
  return(p_fisher)
}

annotate_tf_results <- function(target_regions_df, tf_names_expressed) {

  # 1. Annotate each aSHM region with A/C/G/T proportions and compute
  #    the per-gene average (aSHM_avg_ATCG_prop) used as motifbreakR background
  target_regions_df <- annotate_nucleotides_prop(target_regions_df)

  # 2. Convert each significant focus's SSM subset to BED format
  target_regions_df <- target_regions_df %>%
    dplyr::mutate(
      matched_mut_foci_data = purrr::map(
        matched_mut_foci_data, ~{
          if (is.null(.x) || nrow(.x) == 0) return(.x)
          .x %>%
            dplyr::mutate(
              Significant_Foci = purrr::map(
                Significant_Foci, ~{
                  if (is.null(.x) || nrow(.x) == 0) return(.x)
                  .x %>%
                    dplyr::mutate(
                      bed = purrr::map(ssm_filtered, bed_maker)
                    )
                }
              )
            )
        }
      )
    )

  # 3. Run motifbreakR for each focus using the gene-level ATCG background
  target_regions_df <- target_regions_df %>%
    dplyr::mutate(
      matched_mut_foci_data = purrr::map2(
        matched_mut_foci_data, aSHM_avg_ATCG_prop, ~{
          if (is.null(.x) || nrow(.x) == 0) return(.x)
          A_val <- .y$A[[1]]
          C_val <- .y$C[[1]]
          G_val <- .y$G[[1]]
          T_val <- .y$T[[1]]
          .x %>%
            dplyr::mutate(
              Significant_Foci = purrr::map(
                Significant_Foci, ~{
                  if (is.null(.x) || nrow(.x) == 0) return(.x)
                  .x %>%
                    dplyr::mutate(
                      motif_plotting_data = purrr::map(
                        bed,
                        ~ tryCatch(
                          get_motifbreakr(.x, A = A_val, C = C_val, G = G_val, T = T_val),
                          error = function(e) NULL
                        )
                      ),
                      TF_Results = purrr::map(
                        motif_plotting_data, ~tryCatch(
                          as.data.frame(.x, row.names = NULL),
                          error = function(e) tibble::tibble()
                        )
                      )
                    )
                }
              )
            )
        }
      )
    )

  # 4. Clean TF results: drop NA gene symbols and deduplicate per SNP/TF pair
  target_regions_df <- target_regions_df %>%
    dplyr::mutate(
      matched_mut_foci_data = purrr::map(
        matched_mut_foci_data, ~{
          if (is.null(.x) || nrow(.x) == 0) return(.x)
          .x %>%
            dplyr::mutate(
              Significant_Foci = purrr::map(
                Significant_Foci, ~{
                  if (is.null(.x) || nrow(.x) == 0) return(.x)
                  .x %>%
                    dplyr::mutate(
                      TF_Results = purrr::map(
                        TF_Results, ~{
                          tryCatch({
                            if (is.null(.x) || nrow(.x) == 0) return(.x)
                            .x %>%
                              dplyr::filter(geneSymbol != "NA") %>%
                              dplyr::group_by(SNP_id, geneSymbol) %>%
                              dplyr::slice_head(n = 1) %>%
                              dplyr::ungroup()
                          }, error = function(e) tibble::tibble())
                        }
                      )
                    )
                }
              )
            )
        }
      )
    )

  # 5. Retain only TFs that are expressed above threshold in this pathology
  target_regions_df <- target_regions_df %>%
    dplyr::mutate(
      matched_mut_foci_data = purrr::map(
        matched_mut_foci_data, ~{
          if (is.null(.x) || nrow(.x) == 0) return(.x)
          .x %>%
            dplyr::mutate(
              Significant_Foci = purrr::map(
                Significant_Foci, ~{
                  if (is.null(.x) || nrow(.x) == 0) return(.x)
                  .x %>%
                    dplyr::mutate(
                      TF_Results = purrr::map(
                        TF_Results, ~{
                          tryCatch({
                            if (is.null(.x) || nrow(.x) == 0) return(.x)
                            .x %>%
                              dplyr::mutate(
                                TF_name = geneSymbol,
                                TF_status = ifelse(
                                  TF_name %in% tf_names_expressed,
                                  "presence",
                                  "absence"
                                )
                              ) %>%
                              dplyr::filter(TF_status == "presence")
                          }, error = function(e) tibble::tibble())
                        }
                      )
                    )
                }
              )
            )
        }
      )
    )

  # 6. Fisher test: compare TF binding-site mutation rate inside vs. outside
  #    each significant focus; BH-adjust and filter to adjusted p < 0.05
  target_regions_df <- target_regions_df %>%
    dplyr::mutate(
      matched_mut_foci_data = purrr::pmap(
        list(matched_mut_foci_data, ssm, aSHM_avg_ATCG_prop),
        function(matched_mut_foci_data, ssm, aSHM_avg_ATCG_prop) {

          if (is.null(matched_mut_foci_data) || nrow(matched_mut_foci_data) == 0) {
            return(matched_mut_foci_data)
          }

          A_val <- aSHM_avg_ATCG_prop$A[[1]]
          C_val <- aSHM_avg_ATCG_prop$C[[1]]
          G_val <- aSHM_avg_ATCG_prop$G[[1]]
          T_val <- aSHM_avg_ATCG_prop$T[[1]]

          bed_all <- tryCatch(bed_maker(ssm), error = function(e) NULL)
          full_motifbr <- if (is.null(bed_all) || nrow(bed_all) == 0) {
            tibble::tibble()
          } else {
            tryCatch(
              as.data.frame(
                get_motifbreakr(bed_all, A = A_val, C = C_val, G = G_val, T = T_val),
                row.names = NULL
              ),
              error = function(e) tibble::tibble()
            )
          }

          mutation_positions <- ssm$Start_Position

          matched_mut_foci_data %>%
            dplyr::mutate(
              Significant_Foci = purrr::map2(
                Significant_Foci, lambda_name,
                ~{
                  sig_foci_df <- .x
                  lambda_name <- .y

                  if (is.null(sig_foci_df) || nrow(sig_foci_df) == 0 ||
                      !"TF_Results" %in% names(sig_foci_df)) {
                    return(sig_foci_df)
                  }

                  sig_foci_df %>%
                    dplyr::mutate(
                      TF_Results = purrr::map2(
                        TF_Results, foci,
                        ~{
                          tf_df    <- .x
                          foci_id  <- .y

                          if (is.null(tf_df) || nrow(tf_df) == 0) return(tf_df)
                          if (foci_id == "CN") return(tf_df)

                          group_col <- paste0("group_", lambda_name)
                          if (!group_col %in% names(ssm)) return(tf_df)

                          num_mut_foci <- ssm %>%
                            dplyr::filter(.[[group_col]] == foci_id) %>%
                            nrow()

                          out_foci_mutations <- ssm %>%
                            dplyr::filter(.[[group_col]] != foci_id) %>%
                            dplyr::ungroup()

                          foci_positions <- ssm %>%
                            dplyr::filter(.[[group_col]] == foci_id) %>%
                            dplyr::pull(Start_Position)

                          res_motifbr <- tryCatch(
                            full_motifbr %>%
                              dplyr::filter(
                                !as.character(start) %in% as.character(foci_positions)
                              ),
                              
                            error = function(e) tibble::tibble()
                          )

                          tf_df_counts <- tf_df %>%
                            dplyr::mutate(
                              start_num = as.numeric(start),
                              end_num   = as.numeric(end),
                              num_mutations = purrr::map2_int(
                                start_num, end_num,
                                ~ sum(mutation_positions >= .x & mutation_positions <= .y,
                                      na.rm = TRUE)
                              )
                            ) %>%
                            dplyr::select(-start_num, -end_num)

                          tf_summary <- tf_df_counts %>%
                            dplyr::select(start, end, TF_name, num_mutations) %>%
                            unique() %>%
                            dplyr::group_by(TF_name) %>%
                            dplyr::summarise(n_mut = sum(num_mutations), .groups = "drop") %>%
                            dplyr::rowwise() %>%
                            dplyr::mutate(
                              raw_p_value = tryCatch(
                                fisher_aid_p(
                                  num_tf_mut_target = tibble::tibble(TF_name = TF_name,
                                                                     n_mut   = n_mut),
                                  num_mut_target    = num_mut_foci,
                                  maf_out_foci      = out_foci_mutations,
                                  motifbr           = res_motifbr
                                ),
                                error = function(e) NA_real_
                              )
                            ) %>%
                            dplyr::ungroup() %>%
                            dplyr::mutate(
                              adjusted_p_value = p.adjust(raw_p_value, method = "BH")
                            )

                          tf_df_counts %>%
                            dplyr::left_join(
                              tf_summary %>%
                                dplyr::select(TF_name, n_mut, raw_p_value, adjusted_p_value),
                              by = "TF_name"
                            ) %>%
                            dplyr::filter(adjusted_p_value < 0.05)
                        }
                      )
                    )
                }
              )
            )
        }
      )
    )

  return(target_regions_df)
}