#!/usr/bin/env Rscript
#
# Notes: 
#   This script is intended for use with the dlbclone-1.0 module in LCR-modules. 
#   ...
#
#   ...

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressPackageStartupMessages({
  library(GAMBLR.utils)
  library(optparse)
  library(tibble)
  library(dplyr)
  library(readr)
})

# Determine arguments -----------------------------------------------------

# Define command line options
option_list <- list(
    make_option(c("-m", "--metadata"), type = "character", help = "path to test metadata file"),
    make_option(c("-a", "--metadata_sample_id_colname"), type = "character", help = "Column name in metadata file for sample IDs"),
    make_option(c("-c", "--sv_from_metadata"), type = "character", help = "A named vector that specifies the columns containing the oncogene translocation status for any SV that is annotated in the metadata"),
    make_option(c("-d", "--translocation_status"), type = "character", help = "A named vector that specifies the values corresponding to positive translocation, negative translocation, and NA for any SV that is annotated in the metadata"),    
    make_option(c("-f", "--maf"), type = "character", help = "path to test maf file"),
    make_option(c("-e", "--maf_sample_id_colname"), type = "character", help = "Column name in MAF file for sample IDs"),

    make_option(c("-o", "--output_matrix_dir"), type = "character", help = "output mutation matrix path"),
    make_option(c("-i", "--genes_noncoding"), type = "character", help = "Vector of gene symbols for synonymous mutations"),
    make_option(c("-j", "--genes_hotspot"), type="character", default="MYD88", help = "Vector specifying genes for which hotspot mutations should be separately annotated"),
    make_option(c("-z", "--genes_driver"), type = "character", default=c("CD79B","EZH2","NOTCH1","NOTCH2"), help = "Vector of gene symbols for known driver mutations"),
    make_option(c("-g", "--genes_coding"), type = "character", help = "Vector of gene symbols to include"),  
    make_option(c("-k", "--sv_value"), type = "integer", help = "structural variant value default: 2"),
    make_option(c("-l", "--noncoding_value"), type = "integer", help = "synonymous mutation value default: 1"),
    make_option(c("-n", "--coding_value"), type = "integer", help = "coding mutation value default: 2"),
    make_option(c("-r", "--driver_value"), type = "integer", help = "driver mutation value default: 2")
)

opt <- parse_args(OptionParser(option_list=option_list))

parse_nullable <- function(x) {
    if (x == "None") {
        return(NULL)
    }
    eval(parse(text = x))
}

# command line options; formatting meta/maf files to standard column names -------
metadata <- opt$metadata
metadata_sample_id_colname <- opt$metadata_sample_id_colname
sv_from_metadata <- parse_nullable(opt$sv_from_metadata)

trans_raw <- parse_nullable(opt$translocation_status)
if ("NA_STATUS" %in% names(trans_raw)) { # Convert special placeholder to actual NA
    names(trans_raw)[names(trans_raw) == "NA_STATUS"] <- NA
}
value_map <- setNames(names(trans_raw), unlist(trans_raw)) # Build reverse map: e.g. "yes" -> "POS", "no" -> "NEG", "Not Available" -> NA
maf <- opt$maf
maf_sample_id_colname <- opt$maf_sample_id_colname

# formatting meta/maf files to standard column names -----------------------------
sv_cols <- unlist(sv_from_metadata, use.names = FALSE)
test_metadata <- readr::read_tsv(
        metadata
    ) %>%
    rename(
        sample_id = metadata_sample_id_colname
    ) %>%
    mutate(across( # Assign POS, NEG, NA values
        all_of(sv_cols),
        ~ dplyr::recode(.x, !!!value_map, .default = .x)
    ))  

test_maf <- readr::read_tsv(
        maf
    ) %>%  
    mutate(
        sample_id = .data[[maf_sample_id_colname]]
    )

# command line options; assemble mutation matrix ---------------------------------
output_matrix_dir <- opt$output_matrix_dir
genes_noncoding <- parse_nullable(opt$genes_noncoding)
genes_hotspot <- parse_nullable(opt$genes_hotspot)
genes_driver <- parse_nullable(opt$genes_driver)
if (is.null(genes_driver)) {
     stop("Error: genes_driver cannot be NULL. Please provide a vector of known driver genes to include in the feature matrix.")
}
genes_coding <- parse_nullable(opt$genes_coding)
if (is.null(genes_coding)) {
    genes_coding <- c(
        "TNFRSF14","TP73","CCDC27","NOL9","MTOR","SPEN","ID3","ARID1A","ZC3H12A","RRAGC",
        "BCL10","CD58","NOTCH2","HIST2H2BE","CTSS","SEMA4A","FCGR2B","FCRLA","PRRC2C","SMG7",
        "BTG2","ITPKB","IRF2BP2","HNRNPU","SEC24C","PTEN","FAS","LCOR","EDRF1","WEE1","MPEG1",
        "MS4A1","CCND1","FADD","SESN3","BIRC3","ATM","POU2AF1","ETS1","WNK1","PTPN6","ETV6",
        "KRAS","KMT2D","STAT6","BTG1","ATP2A2","HVCN1","DTX1","SETD1B","BCL7A","FOXO1","RB1",
        "CHD8","NFKBIA" ,"ZFP36L1","TCL1A","YY1","CRIP1","TUBGCP5","MGA","B2M","RFX7","MAP2K1",
        "SIN3A","IDH2","CHD2","CREBBP","TFAP4","USP7","CIITA","SOCS1","IL4R","CD19","PLCG2",
        "MBTPS1","IRF8","TP53","IKZF3","STAT5B","STAT3","CD79B","GNA13","GRB2","ACTG1","KLHL14",
        "BCL2","MED16","TCF3","CD70","S1PR2","SMARCA4","JUNB","KLF2","JAK3","MEF2B","CD79A",
        "POU2F2","DNMT3A","BIRC6","FBXO11","BCL11A","XPO1","PCBP1","DYSF","DUSP2","CXCR4","ZEB2",
        "METAP1D","SF3B1","BCL2L1","PTPN1","BMP7","PRDM15","IGLL5","XBP1","RAC2","EP300","RFTN1",
        "OSBPL10","MYD88","SETD2","RHOA","GNAI2","NFKBIZ","TBL1XR1","KLHL6","EIF4A2","BCL6",
        "HNRNPD","TET2","FBXW7","MEF2C","ZNF608","TNIP1","EBF1","IRF4","FOXC1","CD83","HIST1H3B",
        "HIST1H1C","HIST1H2BC","HIST1H2AC","HIST1H1E","HIST1H2BG","HIST1H1D","HIST1H2AG","HIST1H2BK",
        "HIST1H1B","HIST1H2AM","HLA-A","HLA-C","HLA-B","LTB","PRRC2A","HLA-DMB","PIM1","CCND3",
        "TBCC","NFKBIE","EEF1A1","TMEM30A","ZNF292","BACH2","PRDM1","SGK1","TNFAIP3","INTS1","CARD11",
        "ACTB","PIK3CG","LRRN3","POT1","BRAF","EZH2","KMT2C","ATP6V1B2","PRKDC","LYN","TOX","PABPC1",
        "UBR5","MYC","CD274","CDKN2A","GRHPR","NOTCH1","TRAF2","P2RY8","TMSB4X","DDX3X","PIM2","TAF1",
        "BTK","UBE2A","PHF6","VMA21","ATP6AP1"
    )
}
sv_value <- opt$sv_value
noncoding_value <- opt$noncoding_value
driver_value <- opt$driver_value
coding_value <- opt$coding_value

# assemble mutation matrix -------------------------------------------------------

test_maf <- create_maf_data( #convert to S3 object for compatability
        maf_df = test_maf,
        genome_build = "grch37"
    ) %>%
    mutate(Chromosome = as.character(Chromosome))

test_matrix = GAMBLR.utils::assemble_genetic_features(
    these_samples_metadata = test_metadata,
    unannotated_maf = test_maf,
    sv_from_metadata = sv_from_metadata,
    genes_coding = genes_coding,
    genes_noncoding = genes_noncoding,
    genes_hotspot = genes_hotspot,
    genes_driver = genes_driver, 
    encoding_policy = list(coding=coding_value,driver=driver_value,noncoding=noncoding_value,sv=sv_value) 
    ) %>%
    rownames_to_column(var = "sample_id")

write_tsv(test_matrix, output_matrix_dir)
