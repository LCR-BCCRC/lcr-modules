library(Battenberg)
library(optparse)
#source("./src/R/battenberg/R/clonal_ascat.R")
#source("./src/R/battenberg/R/impute.R")
#devtools::load_all(path="/projects/rmorin/projects/gambl-repos/gambl-rmorin/src/R/battenberg")
#source("/projects/rmorin/projects/gambl-repos/gambl-rmorin/src/R/prepare_wgs.R")
option_list = list(
  make_option(c("-t", "--tumourname"), type="character", default=NULL, help="Samplename of the tumour", metavar="character"),
  make_option(c("-n", "--normalname"), type="character", default=NULL, help="Samplename of the normal", metavar="character"),
  make_option(c("--tb"), type="character", default=NULL, help="Tumour BAM file", metavar="character"),
  make_option(c("--nb"), type="character", default=NULL, help="Normal BAM file", metavar="character"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Directory where output will be written", metavar="character"),
  make_option(c("--skip_allelecount"), type="logical", default=FALSE, action="store_true", help="Provide when alleles don't have to be counted. This expects allelecount files on disk", metavar="character"),
  make_option(c("--skip_preprocessing"), type="logical", default=FALSE, action="store_true", help="Provide when pre-processing has previously completed. This expects the files on disk", metavar="character"),
  make_option(c("--skip_phasing"), type="logical", default=FALSE, action="store_true", help="Provide when phasing has previously completed. This expects the files on disk", metavar="character"),
  make_option(c("--cpu"), type="numeric", default=8, help="The number of CPU cores to be used by the pipeline (Default: 8)", metavar="character"),
  make_option(c("--bp"), type="character", default=NULL, help="Optional two column file (chromosome and position) specifying prior breakpoints to be used during segmentation", metavar="character"),
  make_option(c("-f","--reference_fasta"), type="character", default=NULL, help="Path to indexed genome fasta file (needed for CRAM compatability)", metavar="character"),
  make_option(c("--chr_prefixed_genome"), type="logical", default=FALSE, action="store_true", help="Flag to specify if the genome has chr prefixes in chromosome names", metavar="character"),
  make_option(c("--impute_log"), type="character", default="./", help="Full path for where to store impute logs. If blank, these will be written to the main output directory and cleared."),
  make_option(c("--ref"), type="character", default=NULL, help="Path to reference file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

REFERENCE_BASE = paste0(normalizePath(original_dir,"\\"), "/",opt$ref)
TUMOURNAME = opt$tumourname
NORMALNAME = opt$normalname

REFERENCE_FASTA = opt$reference_fasta
print(paste("using fasta:",REFERENCE_FASTA))
IS.MALE = opt$sex=="male" | opt$sex=="Male"
RUN_DIR = opt$o
CHR_PREFIXED = opt$chr_prefixed_genome
print(paste("chr prefix present?",CHR_PREFIXED))
SKIP_ALLELECOUNTING = opt$skip_allelecount
SKIP_PREPROCESSING = opt$skip_preprocessing
SKIP_PHASING = opt$skip_phasing
NTHREADS = opt$cpu
PRIOR_BREAKPOINTS_FILE = opt$bp
IMPUTE_LOG = opt$impute_log
verbose = TRUE
###############################################################################
# 2018-11-01
# A pure R Battenberg v2.2.9 WGS pipeline implementation.
# sd11 [at] sanger.ac.uk
###############################################################################

# General static
IMPUTEINFOFILE = paste0(REFERENCE_BASE,"/impute_info.txt")
print(IMPUTEINFOFILE)
REPLICCORRECTPREFIX = paste0(REFERENCE_BASE, "/battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_")
G1000PREFIX = paste0(REFERENCE_BASE, "/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr")
G1000PREFIX_AC = paste0(REFERENCE_BASE, "/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr")
GCCORRECTPREFIX = paste0(REFERENCE_BASE, "/battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_")
IMPUTE_EXE = "impute2"  #install using conda

PLATFORM_GAMMA = 1
PHASING_GAMMA = 1
SEGMENTATION_GAMMA = 10
SEGMENTATIIN_KMIN = 3
PHASING_KMIN = 1
CLONALITY_DIST_METRIC = 0
ASCAT_DIST_METRIC = 1
MIN_PLOIDY = 1.6
MAX_PLOIDY = 4.8
MIN_RHO = 0.1
MIN_GOODNESS_OF_FIT = 0.63
BALANCED_THRESHOLD = 0.51
MIN_NORMAL_DEPTH = 10
MIN_BASE_QUAL = 20
MIN_MAP_QUAL = 35
CALC_SEG_BAF_OPTION = 3

# WGS specific static
ALLELECOUNTER = "alleleCounter" #conda package that should have this: cancerit-allelecount
PROBLEMLOCI = paste0(REFERENCE_BASE, "/probloci.txt.gz")

print(PROBLEMLOCI);

# Change to work directory and load the chromosome information
original_dir = getwd()
setwd(RUN_DIR)
NORMALBAM = paste0(normalizePath(original_dir,"\\"), "/",opt$nb)
TUMOURBAM = paste0(normalizePath(original_dir,"\\"), "/",opt$tb)


#this should be the full path to the files after changing directories

#debugging lines added here:
#SKIP_ALLELECOUNTING = TRUE
#SKIP_PREPROCESSING = TRUE
#SKIP_PHASING = FALSE

battenberg(tumourname=TUMOURNAME, 
           normalname=NORMALNAME, 
           tumour_data_file=TUMOURBAM, 
           normal_data_file=NORMALBAM, 
           ismale=IS.MALE, 
           imputeinfofile=IMPUTEINFOFILE, 
           g1000prefix=G1000PREFIX, 
           g1000allelesprefix=G1000PREFIX_AC, 
           gccorrectprefix=GCCORRECTPREFIX, 
           repliccorrectprefix=REPLICCORRECTPREFIX, 
           problemloci=PROBLEMLOCI, 
           data_type="wgs",
           impute_exe=IMPUTE_EXE,
           allelecounter_exe=ALLELECOUNTER,
           nthreads=NTHREADS,
           platform_gamma=PLATFORM_GAMMA,
           phasing_gamma=PHASING_GAMMA,
           segmentation_gamma=SEGMENTATION_GAMMA,
           segmentation_kmin=SEGMENTATIIN_KMIN,
           phasing_kmin=PHASING_KMIN,
           clonality_dist_metric=CLONALITY_DIST_METRIC,
           ascat_dist_metric=ASCAT_DIST_METRIC,
           min_ploidy=MIN_PLOIDY,
           max_ploidy=MAX_PLOIDY,
           min_rho=MIN_RHO,
           min_goodness=MIN_GOODNESS_OF_FIT,
           uninformative_BAF_threshold=BALANCED_THRESHOLD,
           min_normal_depth=MIN_NORMAL_DEPTH,
           min_base_qual=MIN_BASE_QUAL,
           min_map_qual=MIN_MAP_QUAL,
           calc_seg_baf_option=CALC_SEG_BAF_OPTION,
           skip_allele_counting=SKIP_ALLELECOUNTING,
           skip_preprocessing=SKIP_PREPROCESSING,
           skip_phasing=SKIP_PHASING,
           prior_breakpoints_file=PRIOR_BREAKPOINTS_FILE,
           chr_prefixed=CHR_PREFIXED,
           verbose=verbose,
	         logfile_prefix=IMPUTE_LOG,
           ref_fasta=REFERENCE_FASTA)
