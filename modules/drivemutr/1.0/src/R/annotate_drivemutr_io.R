normalize_build <- function(x) {
  x_clean <- tolower(gsub("[^a-z0-9]", "", as.character(x)))
  if (x_clean %in% c("grch37", "hg19", "19", "37")) {
    return("GRCh37")
  } else if (x_clean %in% c("grch38", "hg38", "38")) {
    return("GRCh38")
  } else {
    stop("NCBI_Build must be one of: GRCh37/hg19/19/37 or GRCh38/hg38/38")
  }
}

# Parse "chrN:start-end" or "N:start-end" into a list(chr, start, end).
parse_region <- function(region_str) {
  parts <- strsplit(
    gsub("^(?:chr)?([^:]+):(\\d+)-(\\d+)$", "\\1|\\2|\\3", region_str, perl = TRUE),
    "\\|"
  )[[1]]
  if (length(parts) != 3) stop("Cannot parse region: ", region_str)
  list(
    chr   = parts[1],
    start = as.integer(parts[2]),
    end   = as.integer(parts[3])
  )
}



#   pathology          - string label stored in the Pathology column
#   genome_build       - "grch37"/"hg19"/"37" or "grch38"/"hg38"/"38"
#   genes_regions_list - list of c(Hugo_Symbol, "chrN:start-end") pairs
#   ssm                - data.frame in MAF format (all samples, all regions)
#   cnv                - data.frame with columns: sample_id, <gene1>, <gene2>, …
#                        (one column per gene, values = copy number).
#                        
#   ashm_regions       - logical; if TRUE, appends aSHM regions from
#                        GAMBLR.data (the data package does NOT need GAMBL
#                        infrastructure, only the R package to be installed)
#
# Returns:
#   A tibble matching the get_target_ssm_and_cnv() output:
#   Pathology | NCBI_Build | Hugo_Symbol | Gene_Region | Region | CN (list-col) | ssm (list-col)

annotate_drivemutr_io <- function(pathology,
                                  genome_build       = "grch37",
                                  genes_regions_list = list(),
                                  ssm,
                                  cnv,
                                  ashm_regions       = TRUE,
                                  sample_id_map      = NULL) {
  
  NCBI_Build <- normalize_build(genome_build)
  
  
  target_regions <- purrr::map_dfr(genes_regions_list, ~ tibble::tibble(
    Hugo_Symbol = .x[1],
    Gene_Region = .x[1],
    Region      = .x[2]
  ))
  
  if (isTRUE(ashm_regions)) {
    if (!requireNamespace("GAMBLR.data", quietly = TRUE)) {
      warning("GAMBLR.data is not installed; skipping aSHM regions.")
    } else {
      ashm <- if (NCBI_Build == "GRCh37") {
      GAMBLR.data::grch37_ashm_regions %>%
        dplyr::mutate(chr_name = stringr::str_remove(chr_name, "^chr")) %>%
        dplyr::bind_rows(
          tibble::tibble(
            chr_name           = "3",
            hg19_start         = 101569274L,
            hg19_end           = 101579866L,
            gene               = "NFKBIZ",
            region             = "Rest",
            regulatory_comment = NA_character_
          )
        )
      } else {
        GAMBLR.data::hg38_ashm_regions
      }

      
      regions_bed <- if (NCBI_Build == "GRCh37") {
        dplyr::mutate(ashm,
                      name         = paste(gene, region, sep = "_"),
                      region_coord = paste0(chr_name, ":", hg19_start, "-", hg19_end)
        )
      } else {
        dplyr::mutate(ashm,
                      name         = paste(gene, region, sep = "_"),
                      region_coord = paste0(chr_name, ":", hg38_start, "-", hg38_end)
        )
      }
      
      regions_bed <- regions_bed %>%
        dplyr::mutate(Hugo_Symbol = stringr::str_extract(name, "^[^_]+")) %>%
        dplyr::rename(Gene_Region = name, Region = region_coord)
      
      target_regions <- dplyr::bind_rows(target_regions, regions_bed)
    }
  }
  
  target_regions_data <- data.frame(
    Pathology   = pathology,
    NCBI_Build  = NCBI_Build,
    Hugo_Symbol = target_regions$Hugo_Symbol,
    Gene_Region = target_regions$Gene_Region,
    Region      = target_regions$Region,
    stringsAsFactors = FALSE
  )
  
  
  if (!"sample_id" %in% colnames(cnv)) {
    stop("`cnv` must have a 'sample_id' column.")
  }
  target_regions_data <- target_regions_data %>%
    dplyr::mutate(
      CN = purrr::map(Hugo_Symbol, ~ {
        if (.x %in% colnames(cnv)) {
          tibble::tibble(sample_id = cnv$sample_id, CN = cnv[[.x]])
        } else {
          NULL
        }
      })
    )
  
  
  ssm_norm <- ssm %>%
    dplyr::mutate(
      .chr_norm = sub("^chr", "", as.character(Chromosome), ignore.case = TRUE)
    )
  
  target_regions_data <- target_regions_data %>%
    dplyr::mutate(
      ssm = purrr::map(Region, function(region_str) {
        parsed   <- tryCatch(parse_region(region_str), error = function(e) NULL)
        if (is.null(parsed)) return(NULL)
        
        chr_norm <- sub("^chr", "", parsed$chr, ignore.case = TRUE)
        
        result <- ssm_norm %>%
          dplyr::filter(
            .chr_norm    == chr_norm,
            Start_Position >= parsed$start,
            Start_Position <= parsed$end
          ) %>%
          dplyr::select(-.chr_norm)
        
        if (nrow(result) == 0) NULL else result
      })
    )
  
  # Remove regions without any mutations (mirrors get_target_ssm_and_cnv)
  target_regions_data <- target_regions_data %>%
    dplyr::filter(purrr::map_lgl(ssm, ~ !is.null(.x)))
  
  # Annotate each SSM data frame with an IDs list-column so downstream steps
  # can cross-reference DNA barcodes to RNA/expression sample IDs.
  sample_id_aliases <- if (!is.null(sample_id_map)) {
    sample_id_map %>%
      tidyr::pivot_longer(
        cols      = c(mrna_sample_id, genome_sample_id, capture_sample_id),
        names_to  = "id_type",
        values_to = "all_ids"
      )
  } else {
    NULL
  }
  target_regions_data <- target_regions_data %>%
    dplyr::mutate(
      ssm = purrr::map(ssm, annotate_ids, sample_id_aliases = sample_id_aliases)
    )


  return(target_regions_data)
}