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
  make_option(c("--reference"), type="character", default=NULL, help="Path to reference file", metavar="character"),
  make_option(c("-f","--reference_fasta"), type="character", default=NULL, help="Path to indexed genome fasta file (needed for CRAM compatability)", metavar="character"),
  make_option(c("--chr_prefixed_genome"), type="logical", default=FALSE, action="store_true", help="Flag to specify if the genome has chr prefixes in chromosome names", metavar="character"),
  make_option(c("--impute_log"), type="character", default="./", help="Full path for where to store impute logs. If blank, these will be written to the main output directory and cleared.")
  ,make_option(c("--min_ploidy"), type="character", default=NULL, help="Minimum ploidy to pass explicitly to Battenberg (numeric)")
  ,make_option(c("--max_ploidy"), type="character", default=NULL, help="Maximum ploidy to pass explicitly to Battenberg (numeric)")
  ,make_option(c("--ploidy_constraint"), type="character", default=NULL, help="Optional ploidy constraint in the form MIN-MAX (e.g. '1.6-4.8'). This will be parsed unless explicit --min_ploidy/--max_ploidy are provided.")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
original_dir = getwd()

REFERENCE_BASE = opt$reference
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
IMPUTEINFOFILE = paste0(REFERENCE_BASE,"/battenberg_impute_v3/impute_info_fix.txt")
print(IMPUTEINFOFILE)
G1000PREFIX = paste0(REFERENCE_BASE,"/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr")
G1000PREFIX_AC = paste0(REFERENCE_BASE,"/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr")
GCCORRECTPREFIX = paste0(REFERENCE_BASE,"/battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_")
REPLICCORRECTPREFIX = paste0(REFERENCE_BASE,"/battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_")
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

# Resolve input BAM/CRAM paths to absolute paths BEFORE changing to the run directory.
# This avoids relative-paths being interpreted relative to the run directory, which
# can cause alleleCounter to be unable to locate the files when Snakemake passes
# project-relative paths.
if(!is.null(opt$nb) && nzchar(opt$nb)){
  NORMALBAM = tryCatch(normalizePath(opt$nb, winslash="/", mustWork=FALSE), error=function(e) opt$nb)
} else {
  NORMALBAM = opt$nb
}
if(!is.null(opt$tb) && nzchar(opt$tb)){
  TUMOURBAM = tryCatch(normalizePath(opt$tb, winslash="/", mustWork=FALSE), error=function(e) opt$tb)
} else {
  TUMOURBAM = opt$tb
}

# Change to work directory and load the chromosome information
# Make RUN_DIR if it doesn't exist and fail early with a clear message if
# setwd() cannot be performed. This avoids later `getcwd: cannot access parent`
# errors when the process' current directory has been removed or is unavailable.
if(is.null(RUN_DIR) || !nzchar(RUN_DIR)){
  stop("RUN_DIR (--output) was not provided to the runner. Please pass -o/--output.")
}
if(!dir.exists(RUN_DIR)){
  warning(paste0("RUN_DIR does not exist; attempting to create it: ", RUN_DIR))
  dir.create(RUN_DIR, recursive = TRUE, showWarnings = FALSE)
}



#this should be the full path to the files after changing directories

#debugging lines added here:
#SKIP_ALLELECOUNTING = TRUE
#SKIP_PREPROCESSING = TRUE
#SKIP_PHASING = FALSE
# If explicit min/max ploidy values were passed on the command line, use them (overwriting defaults above)
# If a ploidy constraint string was passed (MIN-MAX), parse it first. Explicit --min_ploidy/--max_ploidy override this.
if(!is.null(opt$ploidy_constraint) && nzchar(opt$ploidy_constraint)){
  parts = unlist(strsplit(opt$ploidy_constraint, "-"))
  if(length(parts) == 2){
    pmin = as.numeric(parts[1]); pmax = as.numeric(parts[2])
    if(!is.na(pmin) && !is.na(pmax)){
      MIN_PLOIDY = pmin
      MAX_PLOIDY = pmax
    } else {
      warning(paste0("Invalid numeric ploidy_constraint provided: ", opt$ploidy_constraint))
    }
  } else {
    warning(paste0("Invalid ploidy_constraint format (expected MIN-MAX): ", opt$ploidy_constraint))
  }
}
if(!is.null(opt$min_ploidy) && nzchar(opt$min_ploidy)){
  t = as.numeric(opt$min_ploidy)
  if(!is.na(t)) MIN_PLOIDY = t
}
if(!is.null(opt$max_ploidy) && nzchar(opt$max_ploidy)){
  t = as.numeric(opt$max_ploidy)
  if(!is.na(t)) MAX_PLOIDY = t
}

# Diagnostic pre-flight checks: print resolved paths and verify BAMs and indexes exist
cat("[battenberg pre-flight] Resolved paths:\n")
cat(sprintf("  TUMOURBAM: %s\n", ifelse(is.null(TUMOURBAM), "<NULL>", TUMOURBAM)))
cat(sprintf("  NORMALBAM: %s\n", ifelse(is.null(NORMALBAM), "<NULL>", NORMALBAM)))


check_file = function(p){
  if(is.null(p) || !nzchar(p)) return(FALSE)
  return(file.exists(p) && file.info(p)$size > 0)
}

tum_bam_ok = check_file(TUMOURBAM)
norm_bam_ok = check_file(NORMALBAM)

cat(sprintf("  TUMOURBAM exists and non-empty: %s\n", tum_bam_ok))
cat(sprintf("  NORMALBAM exists and non-empty: %s\n", norm_bam_ok))

# Check for BAM index files (.bai or .crai)
tum_bai = paste0(TUMOURBAM, ".bai")
tum_bai2 = sub("\\.bam$", ".bai", TUMOURBAM)
tum_crai = sub("\\.bam$", ".crai", TUMOURBAM)
norm_bai = paste0(NORMALBAM, ".bai")
norm_bai2 = sub("\\.bam$", ".bai", NORMALBAM)
norm_crai = sub("\\.bam$", ".crai", NORMALBAM)

cat(sprintf("  TUMOUR index candidates exist: %s, %s, %s\n", file.exists(tum_bai), file.exists(tum_bai2), file.exists(tum_crai)))
cat(sprintf("  NORMAL index candidates exist: %s, %s, %s\n", file.exists(norm_bai), file.exists(norm_bai2), file.exists(norm_crai)))

if(!tum_bam_ok || !norm_bam_ok){
  stop(paste0("Pre-flight error: Missing or empty BAMs. Tried:\n",
              "  Tumour: ", TUMOURBAM, " (exists? ", file.exists(TUMOURBAM), ", size: ", ifelse(file.exists(TUMOURBAM), file.info(TUMOURBAM)$size, NA), ")\n",
              "  Normal: ", NORMALBAM, " (exists? ", file.exists(NORMALBAM), ", size: ", ifelse(file.exists(NORMALBAM), file.info(NORMALBAM)$size, NA), ")\n",
              "If these files are symlinks, ensure the symlink targets are visible from the compute node."))
}

# Re-run with explicit bounds (overwrite previous call above by running again with explicit args)
# Ensure a `preprocess/` directory exists and run the allele-counting steps from there
# so per-chromosome alleleFrequencies files are created inside `preprocess/` as expected.
PREPROCESS_DIR = file.path(RUN_DIR, "preprocess")

# Decide where to run battenberg depending on whether preprocessing is being
# performed. If preprocessing is enabled (SKIP_PREPROCESSING == FALSE) we run
# inside `RUN_DIR/preprocess` so alleleCounter output lands there. If
# preprocessing is skipped (fit stage), we run from `RUN_DIR` where the
# Snakemake fit-rule should have symlinked the preserved preprocess files.
if(!SKIP_PREPROCESSING){
  if(!dir.exists(PREPROCESS_DIR)) dir.create(PREPROCESS_DIR, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(PREPROCESS_DIR)
  on.exit({ try(setwd(old_wd), silent=TRUE) }, add=TRUE)
} else {
  # Safety check for fit-stage: ensure preserved preprocess outputs are present
  # in the current RUN_DIR (the Snakemake rule should have symlinked them here).
  expected_files = c(
    file.path(RUN_DIR, sprintf("%s_alleleCounts.tab", TUMOURNAME)),
    file.path(RUN_DIR, sprintf("%s_mutantBAF.tab", TUMOURNAME)),
    file.path(RUN_DIR, sprintf("%s_mutantLogR.tab", TUMOURNAME)),
    file.path(RUN_DIR, sprintf("%s_mutantLogR_gcCorrected.tab", TUMOURNAME)),
    file.path(RUN_DIR, sprintf("%s_normalBAF.tab", TUMOURNAME)),
    file.path(RUN_DIR, sprintf("%s_normalLogR.tab", TUMOURNAME))
  )
  missing = expected_files[!sapply(expected_files, function(p) file.exists(p) && file.info(p)$size > 0)]
  if(length(missing) > 0){
    stop(paste0("Missing required preprocess files for fit stage. Expected the following files to be present (non-empty) in ", RUN_DIR, ":\n", paste(missing, collapse="\n"),
                "\n\nThe Snakemake fit rule should symlink or copy preprocess outputs into the fit directory before invoking this script."))
  }
  setwd_success <- tryCatch({
    setwd(RUN_DIR)
    TRUE
  }, error = function(e){
    warning(paste0("Failed to change working directory to RUN_DIR (", RUN_DIR, "): ", e$message))
    FALSE
  })
  cat(paste0("[battenberg] current working directory: ", tryCatch(getwd(), error=function(e) '<unavailable>'), "\n"))
  if(!setwd_success){
    stop(paste0("Cannot change to RUN_DIR: ", RUN_DIR, ". Aborting to avoid running with an invalid working directory."))
  }
}

battenberg(analysis = "paired",
           samplename = TUMOURNAME,
           normalname = NORMALNAME,
           sample_data_file = TUMOURBAM,
           normal_data_file = NORMALBAM,
           imputeinfofile = IMPUTEINFOFILE,
           g1000prefix = G1000PREFIX,
           problemloci = PROBLEMLOCI,
           gccorrectprefix = GCCORRECTPREFIX,
           repliccorrectprefix = REPLICCORRECTPREFIX,
           g1000allelesprefix = G1000PREFIX_AC,
           ismale = IS.MALE,
           data_type = "wgs",
           impute_exe = IMPUTE_EXE,
           allelecounter_exe = ALLELECOUNTER,
           nthreads = NTHREADS,
           platform_gamma = PLATFORM_GAMMA,
           phasing_gamma = PHASING_GAMMA,
           segmentation_gamma = SEGMENTATION_GAMMA,
           segmentation_kmin = SEGMENTATIIN_KMIN,
           phasing_kmin = PHASING_KMIN,
           clonality_dist_metric = CLONALITY_DIST_METRIC,
           ascat_dist_metric = ASCAT_DIST_METRIC,
           min_ploidy = MIN_PLOIDY,
           max_ploidy = MAX_PLOIDY,
           min_rho = MIN_RHO,
           min_goodness = MIN_GOODNESS_OF_FIT,
           uninformative_BAF_threshold = BALANCED_THRESHOLD,
           min_normal_depth = MIN_NORMAL_DEPTH,
           min_base_qual = MIN_BASE_QUAL,
           min_map_qual = MIN_MAP_QUAL,
           calc_seg_baf_option = CALC_SEG_BAF_OPTION,
           skip_allele_counting = SKIP_ALLELECOUNTING,
           skip_preprocessing = SKIP_PREPROCESSING,
           skip_phasing = SKIP_PHASING,
           prior_breakpoints_file = PRIOR_BREAKPOINTS_FILE)
