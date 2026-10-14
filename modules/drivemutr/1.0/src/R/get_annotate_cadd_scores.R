get_CADD_scores = function(target_regions_df, cadd_path = NULL){
  
  if(target_regions_df$NCBI_Build[1] == "GRCh37"){
    genome_build = "hg19"
  } else {
    genome_build = "hg38"
  }
  if (!is.null(cadd_path) && nzchar(cadd_path)) {
    bw_a <- file.path(cadd_path, "a.bw")
    bw_c <- file.path(cadd_path, "c.bw")
    bw_g <- file.path(cadd_path, "g.bw")
    bw_t <- file.path(cadd_path, "t.bw")
    missing <- Filter(function(f) !file.exists(f), c(bw_a, bw_c, bw_g, bw_t))
    if (length(missing) > 0) stop("CADD bigwig file(s) not found: ", paste(missing, collapse = ", "))
  } else {
    bw_a <- paste0("https://hgdownload.soe.ucsc.edu/gbdb/", genome_build, "/cadd1.7/a.bw")
    bw_c <- paste0("https://hgdownload.soe.ucsc.edu/gbdb/", genome_build, "/cadd1.7/c.bw")
    bw_g <- paste0("https://hgdownload.soe.ucsc.edu/gbdb/", genome_build, "/cadd1.7/g.bw")
    bw_t <- paste0("https://hgdownload.soe.ucsc.edu/gbdb/", genome_build, "/cadd1.7/t.bw")
  }
  
  all_regions <- GenomicRanges::GRanges(target_regions_df$Region)

  # BigWig files use chr-prefixed seqnames; add prefix if missing
  sl <- GenomeInfoDb::seqlevels(all_regions)
  GenomeInfoDb::seqlevels(all_regions) <- ifelse(
    grepl("^chr", sl, ignore.case = TRUE), sl, paste0("chr", sl)
  )
  
  a <- rtracklayer::import(bw_a, which = all_regions)
  c <- rtracklayer::import(bw_c, which = all_regions)
  g <- rtracklayer::import(bw_g, which = all_regions)
  t <- rtracklayer::import(bw_t, which = all_regions)
  
  cadd_a = as.data.frame(a)[, c("seqnames", "start", "score")]
  names(cadd_a) = c("chr", "position", "CADD_A")
  
  cadd_c = as.data.frame(c)[, c("seqnames", "start", "score")]
  names(cadd_c) = c("chr", "position", "CADD_C")
  
  cadd_g = as.data.frame(g)[, c("seqnames", "start", "score")]
  names(cadd_g) = c("chr", "position", "CADD_G")
  
  cadd_t = as.data.frame(t)[, c("seqnames", "start", "score")]
  names(cadd_t) = c("chr", "position", "CADD_T")
  
  cadd_a$chr <- sub("^chr", "", as.character(cadd_a$chr), ignore.case = TRUE)
  cadd_c$chr <- sub("^chr", "", as.character(cadd_c$chr), ignore.case = TRUE)
  cadd_g$chr <- sub("^chr", "", as.character(cadd_g$chr), ignore.case = TRUE)
  cadd_t$chr <- sub("^chr", "", as.character(cadd_t$chr), ignore.case = TRUE)
  
  cadd_data_all <- dplyr::full_join(cadd_a, cadd_c, by = c("chr", "position")) %>%
    dplyr::full_join(cadd_g, by = c("chr", "position")) %>%
    dplyr::full_join(cadd_t, by = c("chr", "position"))
  
  return(cadd_data_all)
}



annotate_CADD_scores = function(cadd_data, maf_data){
  
  if (is.null(cadd_data)){
    stop("CADD data for your region should be provided to this function! You can check `get_CADD_scores()`")
  }
  
  if (is.null(maf_data)) {
    stop("MAF file!!!")
  }
  
  # Capture IDs lookup from original input BEFORE CADD processing, which drops
  # list-columns via summarise for DNP/TNP/ONP rows.  Building it here ensures
  # all-MNP genes (where every post-summarise row has NULL IDs) are also covered.
  id_lookup_original <- if ("IDs" %in% names(maf_data)) {
    maf_data %>%
      dplyr::filter(!purrr::map_lgl(IDs, is.null)) %>%
      dplyr::select(Tumor_Sample_Barcode, IDs) %>%
      dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)
  } else {
    NULL
  }

  tryCatch(
    {
      maf_chr <- maf_data$Chromosome %>%
        as.character() %>%
        sub("^chr", "", ., ignore.case = TRUE) %>%
        unique()
      
      cadd_data <- cadd_data %>%
        dplyr::mutate(
          chr = sub("^chr", "", as.character(chr), ignore.case = TRUE)
        ) %>%
        dplyr::filter(chr %in% maf_chr)
      
      
      # Make MNPs into single nucleotide to assign their score first
      mnp_maf = maf_data %>% filter(Variant_Type %in% c("DNP", "TNP", "ONP"))
      mnp_maf <- mnp_maf %>%
        separate_rows(Reference_Allele, Tumor_Seq_Allele2, sep = "") %>%
        filter(Reference_Allele != "") %>%
        group_by(Tumor_Sample_Barcode,Start_Position) %>%
        mutate(
          Start_Position = Start_Position + row_number() - 1,
          End_Position = Start_Position
        ) %>%
        ungroup()
      
      rest_maf = maf_data %>% filter(!(Variant_Type %in% c("DNP", "TNP", "ONP")))
      maf = bind_rows(rest_maf, mnp_maf)
      
      # Assign CADD score
      maf <- maf %>%
        rowwise() %>%
        mutate(
          CADD_score = case_when(
            Variant_Type %in% c("SNP","DNP", "TNP", "ONP") ~ {
              rc <- filter(cadd_data, position == Start_Position)
              case_when(
                Tumor_Seq_Allele2 == "A" ~ rc$CADD_A,
                Tumor_Seq_Allele2 == "C" ~ rc$CADD_C,
                Tumor_Seq_Allele2 == "G" ~ rc$CADD_G,
                Tumor_Seq_Allele2 == "T" ~ rc$CADD_T,
                TRUE ~ NA_real_  # Ensures a value is always returned
              )
            },
            Variant_Type %in% c("DEL", "INS") ~ {
              rc <- filter(cadd_data, position >= Start_Position, position <= End_Position)
              max(c(rc$CADD_A, rc$CADD_C, rc$CADD_G, rc$CADD_T), na.rm = TRUE)
            },
            TRUE ~ NA_real_
          )
        ) %>%
        ungroup()
      
      
      # Fix the separated maf by joining together "DNP" , "TNP", "ONP" rows
      maf_fix <- maf %>%
        filter(Variant_Type %in% c("DNP", "TNP", "ONP")) %>%
        group_by(Tumor_Sample_Barcode, Variant_Type) %>%
        arrange(Start_Position, .by_group = TRUE) %>%
        
        mutate(
          block_id = cumsum(c(TRUE, diff(Start_Position) != 1))
        ) %>%
        
        group_by(block_id, .add = TRUE) %>%
        mutate(
          row_in_block = row_number(),
          chunk_size = case_when(
            Variant_Type == "DNP" ~ 2L,
            Variant_Type == "TNP" ~ 3L,
            Variant_Type == "ONP" ~ 4L,  # adjust 
            TRUE ~ 1L
          ),
          chunk_id = ceiling(row_in_block / chunk_size)
        ) %>%
        
        group_by(Tumor_Sample_Barcode, Variant_Type, block_id, chunk_id) %>%
        summarise(
          Hugo_Symbol = dplyr::first(Hugo_Symbol),
          NCBI_Build = dplyr::first(NCBI_Build),
          Chromosome = dplyr::first(Chromosome),
          Start_Position     = min(Start_Position),
          End_Position       = max(End_Position),
          Reference_Allele   = paste(Reference_Allele, collapse = ""),
          Tumor_Seq_Allele2  = paste(Tumor_Seq_Allele2, collapse = ""),
          Matched_Norm_Sample_Barcode = dplyr::first(Matched_Norm_Sample_Barcode),
          Variant_Classification = dplyr::first(Variant_Classification),
          CADD_score         = max(CADD_score, na.rm = TRUE),
          .groups            = "drop"
        ) %>% dplyr::select(-block_id, -chunk_id) 
      
      maf_other <- maf %>%
        filter(!Variant_Type %in% c("DNP", "TNP", "ONP"))
      
      maf_original <- bind_rows(maf_fix, maf_other)
      maf = maf_original[,names(maf)]
      
      
      # Scale CADD score to 0-1
      maf$CADD_score = as.numeric(scale(maf$CADD_score))
      sigmoid = function(x) {
        1 / (1 + exp(-x))
      }
      maf$CADD_score = sigmoid(maf$CADD_score)
      maf_data = maf
    }, error = function(e) {
      message("Error in ", maf_data$Hugo_Symbol[[1]])
    }
  )
  
  # Re-attach IDs using the pre-processing lookup (id_lookup_original).
  if (!is.null(id_lookup_original) && nrow(id_lookup_original) > 0) {
    maf_data <- maf_data %>%
      dplyr::select(-dplyr::any_of("IDs")) %>%
      dplyr::left_join(id_lookup_original, by = "Tumor_Sample_Barcode")
  }

  return(maf_data)
}