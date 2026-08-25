#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Tom Jones & Clare Puttick (McGranahanLab/mhc-hammer)
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import os
import re
import time
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version = "1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section


# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mhc_hammer"]`
CFG = op.setup_module(
    name = "mhc_hammer",
    version = "1.0",
    subdirectories = [
        "inputs", "mhc_reference", "preprocess", "hlahd", "patient_reference",
        "novoalign", "allele_bams", "dna_analysis", "mutations", "cohort_tables", "outputs"
    ],
)

# oncopipe's own CFG["paired_runs"] only excludes pair_status == "no_normal" -- it still
# includes pair_status == "unmatched" (oncopipe's `run_unpaired_tumours_with: "unmatched_normal"`
# substitute-normal mechanism, which pairs an otherwise-unpaired tumour with some OTHER patient's
# designated stand-in normal, e.g. `_shared.unmatched_normal_ids`). That mechanism is a reasonable
# approximation for copy-number callers, but HLA typing/CN-AIB/mutation-calling all depend on the
# patient's OWN real germline sample to establish their baseline HLA genotype -- a substitute
# normal from an unrelated patient is meaningless here, and every downstream rule that looks up
# "the patient's germline sample" filters CFG["samples"] directly by patient_id, so it always
# correctly finds zero real germline samples for a genuinely-unpaired patient regardless of
# pairing_config. Without this extra filter, an "unmatched" row survives into CFG["paired_runs"],
# gets requested as a real target by _mhc_hammer_all, and then hard-fails with a germline-sample
# assertion error. Narrow CFG["paired_runs"] in place (CFG is a live reference into
# config["lcr-modules"]["mhc_hammer"], so this propagates to every one of the ~8 call sites below,
# including the ones inside functions that re-fetch CFG from `config` directly) so it means
# exactly "runs where the tumour has its own real matched germline sample."
#
# NOTE: this makes `pairing_config.<seq_type>.run_unpaired_tumours_with` effectively a no-op for
# this module specifically -- "no_normal" and "unmatched_normal" now behave identically (both
# excluded), even though they mean very different things for other modules (e.g. battenberg/
# sequenza actually process "unmatched_normal" tumours against the substitute normal). Leave it
# set to "no_normal" (the shipped default) since that's the semantically-correct choice; setting
# it to "unmatched_normal" changes nothing here. Warn (rather than silently drop) so a tumour
# that some OTHER module happily processes via a substitute normal doesn't just vanish from
# mhc_hammer's cohort table with no visible explanation.
_mhc_hammer_dropped = CFG["paired_runs"][CFG["paired_runs"]["pair_status"] != "matched"]
if len(_mhc_hammer_dropped) > 0:
    print(
        f"WARNING [mhc_hammer]: excluding {len(_mhc_hammer_dropped)} tumour sample(s) with no "
        f"real matched germline sample (pair_status != 'matched', likely paired via oncopipe's "
        f"substitute-normal mechanism elsewhere in the pipeline) -- mhc_hammer cannot type HLA "
        f"alleles without the patient's own germline DNA"
    )
CFG["paired_runs"] = CFG["paired_runs"][CFG["paired_runs"]["pair_status"] == "matched"]

# Unconditional, always-on scope summary -- this module has been repeatedly hard to debug when
# CFG["samples"]/CFG["paired_runs"] turn out empty or unexpectedly scoped for a given invocation
# (e.g. GAMBL's own subset=/config-merge layer, upstream of this file, silently producing a
# different sample table between two runs of the nominally "same" subset). Print what this
# specific invocation actually sees so that's the first thing visible in the log, rather than
# something that has to be re-derived after the fact from a downstream rule's empty input list.
print(
    f"INFO [mhc_hammer]: {len(CFG['samples'])} sample(s) across "
    f"{CFG['samples']['patient_id'].nunique() if len(CFG['samples']) else 0} patient(s), "
    f"{len(CFG['paired_runs'])} matched tumour/normal pair(s) in scope for this invocation. "
)

# v1 is DNA-only (WES): HLA typing -> personalised reference -> Novoalign alignment ->
# allele-specific BAM splitting -> LOH/allelic-imbalance detection -> Mutect2+VEP mutation
# calling. RNA allelic expression/imbalance/repression and alt-splicing are out of scope --
# see modules/mhc_hammer/README.md.
#
# Design note on dynamic HLA alleles: the number/names of HLA alleles for a patient are only
# known after HLA-HD typing runs, and every upstream process that needs to iterate "per allele"
# does so with an internal bash/R loop that parses allele names from the personalised reference
# FASTA headers at runtime, rather than spawning one task per allele. This module preserves that
# design: rules that produce allele-named files (rules whose output filenames depend on the
# runtime HLA genotype) always also write at least one fixed-name manifest file in the same
# output directory, and every downstream consumer locates the allele-named siblings via
# `dirname()` + shell-side allele discovery -- never via individual Snakemake output/input
# declarations. This avoids needing a `checkpoint` anywhere in the DAG.

SCRIPTS_DIR = CFG["options"]["mhc_hammer_scripts_dir"]
NOVOALIGN_DIR = CFG["options"]["novoalign_dir"]
HLAHD_DIR = CFG["options"]["hlahd_dir"]
VEP_PATH = CFG["options"]["vep_path"]
VEP_CACHE = CFG["inputs"]["vep_cache"]
# Precomputed once here (not read from CFG inside the lambdas/functions below) for the same
# reason -- see the NOTE above _mhc_hammer_get_patient_id_for_sample: any lambda/function
# evaluated lazily (params functions included) must never close over the module-level `CFG` name
# itself, since op.cleanup_module(CFG) deletes it before those lambdas ever run. A plain module-
# level constant assigned here, before cleanup, is a completely different (never-deleted) name and
# is safe to reference from a lambda.
MIN_DEPTH = CFG["options"]["min_depth"]
CONTIG_READS = CFG["options"]["contig_reads"]
MHC_BUILD_COORDINATE_SYSTEM_PATTERNS = CFG["options"]["mhc_build_coordinate_system_patterns"]
MHC_COORDS_BY_SYSTEM = CFG["options"]["mhc_coords_by_system"]

# The chr-prefix half of this lives in _mhc_hammer_subset_bam's shell block (read directly from
# each sample's own BAM header at runtime -- see the long comment on options.mhc_coords_by_system
# in default.yaml for why that's not done here too). This only resolves the coordinate-SYSTEM
# half: which of the two canonical (hg19-style, hg38-style) HLA class I regions applies, via
# regex patterns matched against the raw genome_build string, first match wins.
def _mhc_hammer_get_mhc_region_coords(genome_build):
    system = None
    for pattern, candidate_system in MHC_BUILD_COORDINATE_SYSTEM_PATTERNS.items():
        if re.search(pattern, genome_build):
            system = candidate_system
            break
    if system is None:
        raise ValueError(
            f"mhc_hammer: genome_build '{genome_build}' didn't match any pattern in "
            f"options.mhc_build_coordinate_system_patterns -- add one in default.yaml (or a "
            f"project config override) that matches it."
        )
    if system not in MHC_COORDS_BY_SYSTEM:
        raise ValueError(
            f"mhc_hammer: genome_build '{genome_build}' matched coordinate system '{system}', "
            f"but options.mhc_coords_by_system has no entry for '{system}'."
        )
    return MHC_COORDS_BY_SYSTEM[system]

# HLA class II germline typing (separate, parallel path -- see the long comment on
# options.hla2_coords_by_system in default.yaml). Reuses MHC_BUILD_COORDINATE_SYSTEM_PATTERNS
# (that hg19-vs-hg38 classification is generic to any chr6 locus) against a different canonical
# coordinate table -- same structure as _mhc_hammer_get_mhc_region_coords above, just resolved
# against HLA2_COORDS_BY_SYSTEM instead of MHC_COORDS_BY_SYSTEM.
HLA2_COORDS_BY_SYSTEM = CFG["options"]["hla2_coords_by_system"]
HLA2_GENES = CFG["options"]["hla2_genes"]
# os.path.abspath(), not just the raw config value: unlike SCRIPTS_DIR (mhc_hammer_scripts_dir,
# always a real absolute path since the user supplies it directly), oncopipe resolves {MODSDIR}
# via os.path.normpath() -- not os.path.abspath() -- so it stays relative if the user's own
# _shared.lcr-modules config value is relative (as it typically is, e.g. "../"). Resolved to
# absolute here, once, at parse time (while the CWD is still guaranteed to be wherever Snakemake
# itself was invoked from -- the same directory _shared.lcr-modules is itself relative to), so
# it survives _mhc_hammer_hla2_hlahd's `cd {params.workdir}` later. Confirmed on a real run
# (2026-08-11): "Fatal error: cannot open file 'src/lcr-modules/modules/mhc_hammer/1.0/src/
# parse_hlahd_output.R': No such file or directory" -- a relative path, valid only from the
# original invocation directory, being read after the shell had already cd'd elsewhere.
HLA2_PARSE_SCRIPT = os.path.abspath(CFG["options"]["hla2_parse_script"])

# Same os.path.abspath() treatment as HLA2_PARSE_SCRIPT above, for the same reason -- this is
# another module-owned {MODSDIR}-templated script path.
MUTATIONS_TO_MAF_SCRIPT = os.path.abspath(CFG["options"]["mutations_to_maf_script"])

def _mhc_hammer_get_hla2_region_coords(genome_build):
    system = None
    for pattern, candidate_system in MHC_BUILD_COORDINATE_SYSTEM_PATTERNS.items():
        if re.search(pattern, genome_build):
            system = candidate_system
            break
    if system is None:
        raise ValueError(
            f"mhc_hammer: genome_build '{genome_build}' didn't match any pattern in "
            f"options.mhc_build_coordinate_system_patterns -- add one in default.yaml (or a "
            f"project config override) that matches it."
        )
    if system not in HLA2_COORDS_BY_SYSTEM:
        raise ValueError(
            f"mhc_hammer: genome_build '{genome_build}' matched coordinate system '{system}', "
            f"but options.hla2_coords_by_system has no entry for '{system}'."
        )
    return HLA2_COORDS_BY_SYSTEM[system]

# Upstream's own scripts hardcode the literal tokens "wxs"/"rnaseq" (its internal DNA-vs-RNA
# vocabulary) into filenames and into the regexes several cohort-level R scripts use to find and
# parse those filenames (e.g. `make_cohort_overview_table.R` matches `_wxs.library_size.txt$` and
# `_wxs_novoalign.hla_bam_read_count.csv$` literally) -- this is independent of and NOT the same
# as lcr-modules' own `seq_type` wildcard (`capture`/`genome`). Confirmed by reading the actual
# scripts in a local clone. Every filename segment that flows into one of those scripts' pattern
# matching must use this literal constant, never `{wildcards.seq_type}` directly.
MHC_SEQ = "wxs"
IMGT_URLS = CFG["options"]["imgt_release_urls"][CFG["options"]["imgt_release"]]

# Define rules to be run locally when using a compute cluster. Beyond the trivial symlink/target
# rules, this also includes rules whose real work is pure small-CSV/text bookkeeping with no raw
# BAM/FASTQ/FASTA IO at all (_mhc_hammer_generate_inventory, _mhc_hammer_mutations_to_maf,
# _mhc_hammer_cohort_table) -- fast enough that submitting them as separate cluster jobs is pure
# scheduling overhead. Rules that touch BAM/FASTQ/FASTA content directly (even small,
# already-subsetted ones, e.g. _mhc_hammer_parse_mutations's per-mutation BAM read-counting) are
# deliberately left off this list rather than assumed fast.
localrules:
    _mhc_hammer_input_bam,
    _mhc_hammer_download_reference,
    _mhc_hammer_generate_inventory,
    _mhc_hammer_mutations_to_maf,
    _mhc_hammer_cohort_table,
    _mhc_hammer_output_dna_analysis,
    _mhc_hammer_output_mutations,
    _mhc_hammer_output_mutations_maf,
    _mhc_hammer_output_cohort_table,
    _mhc_hammer_output_hla_final_result,
    _mhc_hammer_output_hla2_alleles,
    _mhc_hammer_all,


##### HELPER FUNCTIONS #####


# NOTE: every helper function below re-fetches `CFG` from the persistent `config` global
# (rather than closing over the module-level `CFG` variable) because Snakemake evaluates input/
# param functions lazily, after this whole file -- including the `op.cleanup_module(CFG)` call
# at the very end, which deletes the module-level `CFG` shortcut -- has already been parsed.
# This mirrors the identical pattern in modules/gridss/2.0/gridss.smk's `get_input_per_patient`/
# `get_split_vcf` helpers. Confirmed necessary by dry-run (a closure-based version fails with
# `NameError: name 'CFG' is not defined` as soon as any rule beyond input-symlinking is planned).

def _mhc_hammer_get_patient_id_for_sample(sample_id):
    # Every sample_id is unique within a seq_type, so a plain lookup against CFG["samples"]
    # (not CFG["runs"]) is sufficient regardless of whether the sample is tumour or normal.
    CFG = config["lcr-modules"]["mhc_hammer"]
    hits = op.filter_samples(CFG["samples"], sample_id = sample_id)
    return hits["patient_id"].tolist()[0]

def _mhc_hammer_get_patient_id_for_tumour(tumour_id, seq_type):
    # Mirrors the gridss.smk `get_split_vcf` pattern: recover patient_id for a pair-level rule
    # (wildcarded on tumour_id/normal_id, not patient_id) via CFG["runs"].
    CFG = config["lcr-modules"]["mhc_hammer"]
    hits = op.filter_samples(CFG["runs"], tumour_sample_id = tumour_id, tumour_seq_type = seq_type)
    return hits["tumour_patient_id"].tolist()[0]

def _mhc_hammer_get_germline_sample_id(wildcards):
    # Shared by both the class I and class II (HLA-HD-only) germline-fastq lookups below --
    # every patient must have exactly one germline WES sample (v1 simplification -- see README).
    CFG = config["lcr-modules"]["mhc_hammer"]
    candidates = op.filter_samples(
        CFG["samples"],
        patient_id = wildcards.patient_id,
        seq_type = wildcards.seq_type,
        genome_build = wildcards.genome_build,
    )
    normal = candidates[candidates["tissue_status"].str.lower() == "normal"]
    if len(normal) != 1:
        # oncopipe's own pairing (CFG["paired_runs"], which is what actually requested this
        # patient_id/seq_type/genome_build combination in the first place -- see the
        # pair_status == "matched" narrowing near the top of this file) can consider a tumour
        # "matched" using logic that doesn't line up 1:1 with this direct re-query against
        # CFG["samples"]. Rather than a bare "found 0", dump every sample row this module can see
        # for this patient/seq_type/genome_build regardless of tissue_status, so a genome_build or
        # tissue_status value mismatch on the germline sample's own row is visible immediately
        # instead of requiring a second round-trip.
        all_patient_samples = op.filter_samples(CFG["samples"], patient_id = wildcards.patient_id)
        assert False, (
            f"Expected exactly one germline WES sample for patient '{wildcards.patient_id}' "
            f"({wildcards.seq_type}--{wildcards.genome_build}), found {len(normal)}. This module "
            f"requires exactly one germline WES sample per patient (see README).\n"
            f"Sample(s) matching patient_id + seq_type '{wildcards.seq_type}' + genome_build "
            f"'{wildcards.genome_build}' (any tissue_status): "
            f"{candidates[['sample_id', 'tissue_status']].to_dict('records') if len(candidates) else '(none)'}\n"
            f"ALL sample(s) for patient_id '{wildcards.patient_id}' regardless of seq_type/"
            f"genome_build: "
            f"{all_patient_samples[['sample_id', 'seq_type', 'genome_build', 'tissue_status']].to_dict('records') if len(all_patient_samples) else '(none)'}"
        )
    return normal["sample_id"].tolist()[0]

def _mhc_hammer_get_germline_fqs(wildcards):
    sample_id = _mhc_hammer_get_germline_sample_id(wildcards)
    return {
        "fq1": expand(str(rules._mhc_hammer_generate_fqs.output.fq1), sample_id = sample_id, allow_missing = True),
        "fq2": expand(str(rules._mhc_hammer_generate_fqs.output.fq2), sample_id = sample_id, allow_missing = True),
    }

def _mhc_hammer_get_hla2_germline_fqs(wildcards):
    sample_id = _mhc_hammer_get_germline_sample_id(wildcards)
    return {
        "fq1": expand(str(rules._mhc_hammer_hla2_generate_fqs.output.fq1), sample_id = sample_id, allow_missing = True),
        "fq2": expand(str(rules._mhc_hammer_hla2_generate_fqs.output.fq2), sample_id = sample_id, allow_missing = True),
    }

def _mhc_hammer_reference_dir_for_sample(wildcards):
    patient_id = _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id)
    return expand(str(rules._mhc_hammer_generate_references.output.patient_dir), patient_id = patient_id, allow_missing = True)

def _mhc_hammer_reference_dir_for_tumour(wildcards):
    patient_id = _mhc_hammer_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type)
    return expand(str(rules._mhc_hammer_generate_references.output.patient_dir), patient_id = patient_id, allow_missing = True)

# Battenberg's cellularity_ploidy.txt is genuinely optional (pairs missing it just skip
# copy-number output, matching upstream's own missing_purity_ploidy fallback) -- returning []
# rather than a hard path lets the DAG track it as a real dependency *when present* (so a
# Battenberg rerun correctly invalidates this rule's output) without requiring every pair to
# have one before mhc_hammer can run at all.
def _mhc_hammer_get_battenberg_cp_input(wildcards):
    CFG = config["lcr-modules"]["mhc_hammer"]
    pattern = CFG["options"]["battenberg_cellularity_ploidy"]
    if not pattern:
        return []
    path = pattern.format(
        seq_type = wildcards.seq_type,
        genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id,
        normal_id = wildcards.normal_id
    )
    return [path] if os.path.exists(path) else []

# Parses purity/ploidy out of whatever _mhc_hammer_get_battenberg_cp_input resolved (or "", ""
# if it resolved to nothing). Column names ("cellularity", "ploidy"), not position, matching how
# modules/cnaqc/1.0/src/run_cnaqc.R itself reads the exact same file
# (read.table(header = TRUE, sep = "\t"); cp$cellularity[1]/cp$ploidy[1]).
def _mhc_hammer_parse_cellularity_ploidy(wildcards, input):
    if not input.cellularity_ploidy:
        return "", ""
    with open(input.cellularity_ploidy[0]) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        values = fh.readline().rstrip("\n").split("\t")
    row = dict(zip(header, values))
    return row.get("cellularity", ""), row.get("ploidy", "")

# Cohort-wide equivalent for _mhc_hammer_generate_inventory (no per-pair wildcards there) --
# gathers whichever cellularity_ploidy.txt files actually exist right now across every run in
# CFG["paired_runs"], so the DAG tracks all of them as real dependencies.
def _mhc_hammer_get_all_battenberg_cp_inputs(wildcards):
    CFG = config["lcr-modules"]["mhc_hammer"]
    pattern = CFG["options"]["battenberg_cellularity_ploidy"]
    if not pattern:
        return []
    runs = CFG["paired_runs"]
    paths = []
    for _, row in runs.iterrows():
        path = pattern.format(
            seq_type = row["tumour_seq_type"],
            genome_build = row["tumour_genome_build"],
            tumour_id = row["tumour_sample_id"],
            normal_id = row["normal_sample_id"]
        )
        if os.path.exists(path):
            paths.append(path)
    return paths

def _mhc_hammer_ref_cache_pattern(wildcards):
    CFG = config["lcr-modules"]["mhc_hammer"]
    cache_dir = os.path.abspath(CFG["dirs"]["mhc_reference"] + f"ref_cache/{wildcards.genome_build}")
    return cache_dir + "/%2s/%2s/%s"

# Only pulls in _mhc_hammer_filter_kmers (and transitively _mhc_hammer_build_genome_kmer_index,
# a slow whole-genome jellyfish count+dump over the entire reference) as a real Snakemake
# dependency when fish_reads is actually enabled. Previously _mhc_hammer_subset_bam depended on
# it unconditionally regardless of fish_reads, so disabling kmer-fishing at the options level
# only skipped subset_bam_opt.sh's own internal step -- it never actually skipped the expensive
# upstream build, defeating the point of turning it off for fast iteration. Returns a list (not
# a bare string) so input.kmer_file is consistently indexable either way.
def _mhc_hammer_get_kmer_file_input(wildcards):
    CFG = config["lcr-modules"]["mhc_hammer"]
    if CFG["options"]["fish_reads"]:
        return [str(rules._mhc_hammer_filter_kmers.output.filtered_kmers)]
    return []


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mhc_hammer_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        cram = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.cram",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        # Some "sample_bam" inputs are actually CRAM content symlinked/named as .bam (htslib
        # content-sniffs the real format regardless of extension) -- samtools then looks for a
        # .crai-named sidecar, not .bam.bai, and fails with cram_index_load even though a valid
        # index exists under the "wrong" name. Providing both naming conventions for both the
        # data file and its index covers either case. Mirrors modules/pathseq/1.0/pathseq.smk's
        # _pathseq_input_bam.
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bam, output.cram)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Downloads and unpacks the pinned IMGT-derived MHC Hammer reference bundle from Zenodo
# (docs/mhc_reference_files.md upstream). Two separate zip files per release. Cohort-wide,
# no wildcards -- built once and reused by every patient.
rule _mhc_hammer_download_reference:
    output:
        gtf = CFG["dirs"]["mhc_reference"] + "mhc_references/gtf/mhc.gtf",
        genome_fasta = CFG["dirs"]["mhc_reference"] + "mhc_references/genome/mhc_genome_strand.fasta",
        kmer_file = CFG["dirs"]["mhc_reference"] + "kmer_files/imgt_30mers.fa"
    log:
        stdout = CFG["logs"]["mhc_reference"] + "download_reference.log"
    params:
        mhc_references_url = IMGT_URLS["mhc_references"],
        kmer_files_url = IMGT_URLS["kmer_files"],
        outdir = CFG["dirs"]["mhc_reference"]
    conda:
        CFG["conda_envs"]["mhc_hammer_download"]
    container:
        CFG["container_envs"]["mhc_hammer_download"]
    threads:
        CFG["threads"]["download_reference"]
    resources:
        **CFG["resources"]["download_reference"]
    shell:
        op.as_one_line("""
        (
        mkdir -p {params.outdir} &&
        wget -qO {params.outdir}/mhc_references.zip "{params.mhc_references_url}" &&
        wget -qO {params.outdir}/kmer_files.zip "{params.kmer_files_url}" &&
        unzip -oq {params.outdir}/mhc_references.zip -d {params.outdir} &&
        unzip -oq {params.outdir}/kmer_files.zip -d {params.outdir} &&
        rm {params.outdir}/mhc_references.zip {params.outdir}/kmer_files.zip
        ) > {log.stdout} 2>&1
        """)


# Builds a genome-wide 30-mer occurrence index from a user-supplied reference genome, used by
# _mhc_hammer_filter_kmers below to identify which of the IMGT-derived kmers aren't actually
# HLA-specific (see design note there). Cohort-wide, no wildcards, genome-build-independent by
# design (see inputs.reference_genome_fasta above) -- built once and reused for every patient.
#
# Deliberately NOT canonical (no -C): upstream's own kmer file was built without -C (confirmed
# by reading docs/mhc_reference_files.md), and subset_bam_opt.sh's actual matching is a plain
# `grep -F` literal substring search against SAM SEQ fields (which SAM always reports relative
# to the forward/reference strand regardless of alignment strand) -- so both the kmer file and
# this occurrence count need to use the same strand-specific (non-canonical) convention for the
# comparison to mean anything. Using -C here would count reverse-complement occurrences grep
# could never actually match against, inflating every count.
rule _mhc_hammer_build_genome_kmer_index:
    input:
        genome_fasta = CFG["inputs"]["reference_genome_fasta"]
    output:
        index = CFG["dirs"]["mhc_reference"] + "genome_kmer_index/genome_30mers.jf"
    log:
        stdout = CFG["logs"]["mhc_reference"] + "build_genome_kmer_index.log"
    conda:
        CFG["conda_envs"]["jellyfish"]
    container:
        CFG["container_envs"]["jellyfish"]
    threads:
        CFG["threads"]["build_genome_kmer_index"]
    resources:
        **CFG["resources"]["build_genome_kmer_index"]
    shell:
        op.as_one_line("""
        (
        mkdir -p $(dirname {output.index}) &&
        jellyfish count -m 30 -s 3G -t {threads} -o {output.index} {input.genome_fasta}
        ) > {log.stdout} 2>&1
        """)


# Drops IMGT-derived kmers that occur more than options.kmer_max_genome_occurrences times
# genome-wide (per _mhc_hammer_build_genome_kmer_index above) before they're used for BAM
# subsetting -- a kmer genuinely specific to HLA-A/B/C should occur only a handful of times in
# a reference genome; kmers occurring far more than that are most likely low-complexity/
# repetitive-element contamination rather than real HLA signal, and were observed (2026-07-30,
# on a real capture WES BAM) producing hundreds of thousands of spurious matches per chromosome,
# roughly proportional to chromosome length -- the signature of genome-wide non-specific
# matching rather than true positives. Cohort-wide, no wildcards, replaces the raw downloaded
# kmer file everywhere downstream.
#
# NOTE: correct syntax is `jellyfish query <index> -s <file>` (index positional first, then
# -s) -- an earlier version of this rule had the args swapped ("jellyfish query -s <file>
# <index>"), which ran without error but silently produced zero output lines (confirmed on a
# real run, 2026-07-31). The <kmer> <count>-per-line output format (assumed by the awk below)
# still hasn't been independently verified beyond "it's non-empty and awk's numeric filter
# behaves sanely" -- double check imgt_30mers.genome_occurrences.tsv's actual contents.
rule _mhc_hammer_filter_kmers:
    input:
        index = str(rules._mhc_hammer_build_genome_kmer_index.output.index),
        raw_kmers = str(rules._mhc_hammer_download_reference.output.kmer_file)
    output:
        filtered_kmers = CFG["dirs"]["mhc_reference"] + "kmer_files/imgt_30mers.filtered.fa",
        counts = CFG["dirs"]["mhc_reference"] + "kmer_files/imgt_30mers.genome_occurrences.tsv"
    log:
        stdout = CFG["logs"]["mhc_reference"] + "filter_kmers.log"
    params:
        max_occurrences = CFG["options"]["kmer_max_genome_occurrences"]
    conda:
        CFG["conda_envs"]["jellyfish"]
    container:
        CFG["container_envs"]["jellyfish"]
    threads:
        CFG["threads"]["filter_kmers"]
    resources:
        **CFG["resources"]["filter_kmers"]
    shell:
        # `jellyfish query` (used originally) fails on this whole-genome-scale index:
        # "terminate called after throwing an instance of 'std::runtime_error' what(): Unsupported
        # format" -- confirmed on a real run (2026-07-31) to be a scale limitation in query's
        # random-access-lookup loader, not a corrupt file or config problem: `jellyfish stats`
        # reads the exact same 30GB/~2.52-billion-distinct-kmer index fine, and `query` round-trips
        # correctly on a small test index built the same way. jellyfish's own docs claim count's
        # auto-resize/merge produces output identical in format to a single-pass count, so
        # under-sizing `-s` isn't a fixable explanation either -- this looks like a genuine
        # large-file bug/limitation in query specifically.
        #
        # `jellyfish dump` is a much simpler sequential read (also confirmed working, same as
        # `stats`, at this scale) rather than query's random-access structure, so this rule now
        # streams the *entire* genome dump through awk and intersects it against the (small --
        # thousands of entries) set of IMGT kmers, instead of doing random-access lookups per
        # query kmer. Functionally equivalent to `query -s`, just implemented as a streaming join.
        # This does mean awk has to scan every one of the ~2.5 billion lines jellyfish dump
        # produces -- a real but one-time (cohort-wide, not per-sample) cost.
        op.as_one_line("""
        (
        wd=$(dirname {output.counts}) &&
        grep -v '^>' {input.raw_kmers} | sort -u > $wd/want_kmers.txt &&
        jellyfish dump -c {input.index} | awk -v want_file=$wd/want_kmers.txt
          'BEGIN {{ while ((getline line < want_file) > 0) want[line] = 1 }}
           ($1 in want) {{ print }}' > {output.counts} &&
        rm -f $wd/want_kmers.txt &&
        awk -v max={params.max_occurrences} '$2 <= max {{print $1}}' {output.counts} > {output.filtered_kmers}
        ) > {log.stdout} 2>&1
        """)


# Computes library size (used to normalise depth downstream) on the full input BAM,
# mirroring upstream's FLAGSTAT process. Filename uses the literal MHC_SEQ ("wxs") token, not
# {seq_type} -- see the MHC_SEQ note near the top of this file -- because
# make_cohort_overview_table.R locates this file via a literal `_wxs.library_size.txt$` match.
rule _mhc_hammer_flagstat:
    input:
        bam = str(rules._mhc_hammer_input_bam.output.bam)
    output:
        library_size = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_" + MHC_SEQ + ".library_size.txt",
        flagstat = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_" + MHC_SEQ + ".flagstat"
    log:
        stdout = CFG["logs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/flagstat.log"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    threads:
        CFG["threads"]["flagstat"]
    resources:
        **CFG["resources"]["flagstat"]
    shell:
        # -@ {threads}: this rule reserves {threads} cores via Snakemake but neither samtools
        # subcommand used it by default -- both flagstat and view -c support -@.
        op.as_one_line("""
        (
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat} &&
        samtools view -@ {threads} -c -f 1 -F 2308 {input.bam} > {output.library_size}
        ) > {log.stdout} 2>&1
        """)


# Builds a local htslib MD5-keyed reference cache (via samtools' own seq_cache_populate.pl) from
# the matching reference_files() genome FASTA for {genome_build}. Needed because some
# "sample_bam" inputs are actually CRAM content (see _mhc_hammer_input_bam) -- htslib needs the
# reference sequence to reconstruct read bases from CRAM's reference-based compression, and
# without a local cache it falls back to fetching each reference sequence one at a time from
# EBI's ENA CRAM registry over the network (real symptom hit: `[W::find_file_url] ... Input/output
# error`, silently degrading read extraction rather than hard-failing -- likely the real cause of
# reads being missing/incomplete downstream, separate from the earlier jellyfish kmer-query bug).
# Cohort-wide per genome_build, built once and reused by every sample of that build.
#
# NOTE: seq_cache_populate.pl's exact output layout/behaviour hasn't been independently verified
# by actually running it -- if REF_CACHE-based lookups still fail after this, check that MD5-named
# files actually exist under {params.cache_dir} after this rule runs.
rule _mhc_hammer_build_ref_cache:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        done = touch(CFG["dirs"]["mhc_reference"] + "ref_cache/{genome_build}/.done")
    log:
        stdout = CFG["logs"]["mhc_reference"] + "{genome_build}/build_ref_cache.log"
    params:
        cache_dir = CFG["dirs"]["mhc_reference"] + "ref_cache/{genome_build}"
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["build_ref_cache"]
    resources:
        **CFG["resources"]["build_ref_cache"]
    shell:
        op.as_one_line("""
        (
        mkdir -p {params.cache_dir} &&
        seq_cache_populate.pl -root {params.cache_dir} -subdirs 2 {input.fasta}
        ) > {log.stdout} 2>&1
        """)


# Subsets the input BAM to reads that could plausibly align to the HLA class I genes: the MHC
# region on chr6, unmapped reads, and (fish_reads=True) reads matching the (genome-occurrence-
# filtered -- see _mhc_hammer_filter_kmers) IMGT 30-mer kmer set. Mirrors upstream's SUBSET_BAMS
# process, with the kmer file filtered rather than used raw.
rule _mhc_hammer_subset_bam:
    input:
        bam = str(rules._mhc_hammer_input_bam.output.bam),
        # Declared even though subset_bam_opt.sh only takes -b: guarantees the .bam.crai sidecar
        # (needed when the "bam" is actually CRAM content, see _mhc_hammer_input_bam) exists
        # before this rule's region-based samtools view runs.
        crai = str(rules._mhc_hammer_input_bam.output.crai),
        ref_cache_done = str(rules._mhc_hammer_build_ref_cache.output.done),
        kmer_file = _mhc_hammer_get_kmer_file_input
    output:
        bam = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.subset.sorted.bam",
        bai = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.subset.sorted.bam.bai",
        read_counts = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.read_counts.csv"
    log:
        stdout = CFG["logs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/subset_bam.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        # Coordinate portion only, no chr-prefix -- the shell block below detects that at runtime
        # from this sample's own BAM header (see the long comment on options.mhc_coords_by_system
        # in default.yaml for why).
        mhc_region_coords = lambda w: _mhc_hammer_get_mhc_region_coords(w.genome_build),
        # subset_bam_opt.sh does literal `[ "$var" == true ]` string comparisons (lowercase) --
        # str(True)/str(False) from a Python/YAML bool renders as "True"/"False" (capitalised)
        # when substituted into the shell command, which never matches. Must be lowercased here.
        fish_reads = str(CFG["options"]["fish_reads"]).lower(),
        unmapped_reads = str(CFG["options"]["unmapped_reads"]).lower(),
        # Per-genome_build (not per-seq_type -- see the long comment on options.contig_reads in
        # default.yaml). Falls back to True (upstream's own default) for any genome_build not
        # explicitly listed, same fallback convention as elsewhere in this module. References the
        # CONTIG_READS module-level constant, not CFG directly (CFG-closure gotcha).
        contig_reads = lambda w: str(CONTIG_READS.get(w.genome_build, True)).lower(),
        sort_mem = lambda wildcards, resources: max(1, int(resources.mem_mb / 1000 * 0.8)),
        bam_abs = lambda wildcards, input: os.path.abspath(input.bam),
        # Empty when fish_reads is off (_mhc_hammer_get_kmer_file_input returns []) -- the shell
        # block below creates an empty placeholder file for -k in that case, rather than this
        # trying to abspath a nonexistent input.
        kmer_file_abs = lambda wildcards, input: os.path.abspath(input.kmer_file[0]) if input.kmer_file else "",
        ref_cache_pattern = _mhc_hammer_ref_cache_pattern
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["subset_bam"]
    resources:
        **CFG["resources"]["subset_bam"]
    shell:
        # subset_bam_opt.sh builds its own internal relative paths from `-p` (e.g.
        # tmp_dir="${out_prefix}_tmpDir", bam_header="${tmp_dir}/${out_prefix}.header") -- it
        # expects `-p` to be a bare filename prefix, not a full path, and assumes it's being run
        # from within an already-isolated working directory (same CWD-relative assumption as
        # several other upstream scripts -- see the design note near the top of this file).
        # Passing a full path as `-p` (as an earlier version of this rule did) makes it
        # concatenate the path with itself. `-d`'s filename must be exactly
        # "contigs_placeholder.txt" (not just any empty file) -- the script special-cases that
        # literal name to auto-detect non-standard contigs from the BAM header; any other
        # filename is instead treated as a real (here, empty) user-supplied contig list.
        # subset_bam_opt.sh (bash -eu) skips its own final `rm -rf ${tmp_dir}` cleanup if any
        # step fails partway through, and that tmp_dir isn't a Snakemake-tracked output, so nothing
        # else cleans it up either -- a later retry can then fail on e.g. `samtools merge` refusing
        # to overwrite a stale intermediate file left by the earlier failed attempt. Removing
        # {sample_id}_tmpDir up front makes every invocation of this rule start from a clean slate
        # regardless of prior failed attempts.
        op.as_one_line("""
        wd=$(dirname {output.bam}) && mkdir -p $wd &&
        if samtools view -H {params.bam_abs} | awk -F'\\t' '$1=="@SQ"{{for(i=1;i<=NF;i++) if($i=="SN:chr6") f=1}} END{{exit !f}}'; then
            chr_prefix=chr;
        else
            chr_prefix="";
        fi &&
        echo "${{chr_prefix}}{params.mhc_region_coords}" > $wd/mhc_coords.txt &&
        touch $wd/contigs_placeholder.txt &&
        rm -rf $wd/{wildcards.sample_id}_tmpDir &&
        export REF_CACHE={params.ref_cache_pattern} &&
        kmer_file_path="{params.kmer_file_abs}" &&
        if [ -z "$kmer_file_path" ]; then
            kmer_file_path=$wd/empty_kmers.fa;
            touch $kmer_file_path;
        fi &&
        (
        cd $wd &&
        {params.scripts_dir}/bin/subset_bam_opt.sh
        -b {params.bam_abs}
        -k $kmer_file_path
        -f {params.fish_reads} -c {params.contig_reads} -d contigs_placeholder.txt
        -u {params.unmapped_reads} -h mhc_coords.txt
        -p {wildcards.sample_id} -t {threads}
        -m {params.sort_mem}G -o false
        ) > {log.stdout} 2>&1
        """)


##### HLA CLASS II GERMLINE TYPING (separate, parallel path) #####


# Subsets the input BAM to reads plausibly overlapping the HLA class II region and unmapped
# reads -- deliberately no kmer-fishing here (fish_reads is hardcoded off): the IMGT kmer set
# _mhc_hammer_filter_kmers builds is derived from MHC Hammer's own class-I-only reference bundle,
# so there is no class-II-specific kmer file to fish with. Writes into its own hla2/ subdirectory
# (not the same directory as _mhc_hammer_subset_bam's output, even though both are per-sample)
# specifically to avoid a real collision risk: subset_bam_opt.sh's -d contigs_placeholder.txt
# argument must be that exact literal filename (see the note on _mhc_hammer_subset_bam above), so
# if this rule wrote into the same directory, both rules would race to create/read the
# identically-named contigs_placeholder.txt/mhc_coords.txt/{sample_id}_tmpDir for the same sample
# whenever Snakemake happens to schedule them concurrently (they're independent DAG branches with
# no dependency on each other).
rule _mhc_hammer_hla2_subset_bam:
    input:
        bam = str(rules._mhc_hammer_input_bam.output.bam),
        crai = str(rules._mhc_hammer_input_bam.output.crai),
        ref_cache_done = str(rules._mhc_hammer_build_ref_cache.output.done)
    output:
        bam = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/hla2/{sample_id}.subset.sorted.bam",
        bai = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/hla2/{sample_id}.subset.sorted.bam.bai",
        read_counts = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/hla2/{sample_id}.read_counts.csv"
    log:
        stdout = CFG["logs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/hla2/subset_bam.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        # Coordinate portion only, no chr-prefix -- same runtime BAM-header detection as the
        # class I rule above.
        hla2_region_coords = lambda w: _mhc_hammer_get_hla2_region_coords(w.genome_build),
        unmapped_reads = str(CFG["options"]["unmapped_reads"]).lower(),
        contig_reads = lambda w: str(CONTIG_READS.get(w.genome_build, True)).lower(),
        sort_mem = lambda wildcards, resources: max(1, int(resources.mem_mb / 1000 * 0.8)),
        bam_abs = lambda wildcards, input: os.path.abspath(input.bam),
        ref_cache_pattern = _mhc_hammer_ref_cache_pattern
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["hla2_subset_bam"]
    resources:
        **CFG["resources"]["hla2_subset_bam"]
    shell:
        op.as_one_line("""
        wd=$(dirname {output.bam}) && mkdir -p $wd &&
        if samtools view -H {params.bam_abs} | awk -F'\\t' '$1=="@SQ"{{for(i=1;i<=NF;i++) if($i=="SN:chr6") f=1}} END{{exit !f}}'; then
            chr_prefix=chr;
        else
            chr_prefix="";
        fi &&
        echo "${{chr_prefix}}{params.hla2_region_coords}" > $wd/mhc_coords.txt &&
        touch $wd/contigs_placeholder.txt &&
        rm -rf $wd/{wildcards.sample_id}_tmpDir &&
        export REF_CACHE={params.ref_cache_pattern} &&
        kmer_file_path=$wd/empty_kmers.fa &&
        touch $kmer_file_path &&
        (
        cd $wd &&
        {params.scripts_dir}/bin/subset_bam_opt.sh
        -b {params.bam_abs}
        -k $kmer_file_path
        -f false -c {params.contig_reads} -d contigs_placeholder.txt
        -u {params.unmapped_reads} -h mhc_coords.txt
        -p {wildcards.sample_id} -t {threads}
        -m {params.sort_mem}G -o false
        ) > {log.stdout} 2>&1
        """)


# Converts the class II subset BAM to paired FASTQs -- mirrors _mhc_hammer_generate_fqs exactly.
rule _mhc_hammer_hla2_generate_fqs:
    input:
        bam = str(rules._mhc_hammer_hla2_subset_bam.output.bam)
    output:
        fq1 = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/hla2/{sample_id}.1.fq.gz",
        fq2 = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/hla2/{sample_id}.2.fq.gz"
    log:
        stdout = CFG["logs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/hla2/generate_fqs.log"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    threads:
        CFG["threads"]["hla2_generate_fqs"]
    resources:
        **CFG["resources"]["hla2_generate_fqs"]
    shell:
        op.as_one_line("""
        (
        samtools collate -@ {threads} -u -O {input.bam} |
        samtools fastq -@ {threads} -1 {output.fq1} -2 {output.fq2} -s /dev/null -0 /dev/null -n
        ) > {log.stdout} 2>&1
        """)


# Runs the user-supplied HLA-HD install on the patient's germline sample to type the classical
# class II genes plus HLA-DM/-DO (options.hla2_genes). Patient-level, mirroring
# _mhc_hammer_hlahd's own structure exactly, but: (1) no personalised-reference/novoalign/CN-AIB/
# mutation-calling path follows this -- germline genotyping only, see the module design note at
# the top of this section; (2) uses parse_hlahd_output.R (this module's own script, under src/ --
# not a copy of MHC Hammer's bin/hlahd_parse_output.R) instead of upstream's parser, since that
# script requires a GTF cross-reference this module has no class-II equivalent of. See the long
# comment on that script for what it does differently.
rule _mhc_hammer_hla2_hlahd:
    input:
        unpack(_mhc_hammer_get_hla2_germline_fqs)
    output:
        hla_alleles = CFG["dirs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/hla2/{patient_id}_hla2_alleles.csv",
        result_dir = directory(CFG["dirs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/hla2/result")
    log:
        stdout = CFG["logs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/hla2/hlahd.log"
    params:
        parse_script = HLA2_PARSE_SCRIPT,
        hlahd_dir = HLAHD_DIR,
        workdir = lambda wildcards, output: os.path.dirname(output.hla_alleles),
        keep_intermediates = str(CFG["options"]["keep_hlahd_intermediates"]).lower(),
        min_read_length = CFG["options"]["hlahd_min_read_length"],
        # awk OR-chain selecting exactly the configured gene list from HLA-HD's own stock
        # HLA_gene.split.txt, same technique as the class I rule's A/B/C filter above.
        gene_awk_filter = " || ".join(f'$1 == "{gene}"' for gene in HLA2_GENES),
        genes = " ".join(HLA2_GENES)
    conda:
        CFG["conda_envs"]["mhc_hammer_hlahd"]
    container:
        None
    threads:
        CFG["threads"]["hla2_hlahd"]
    resources:
        **CFG["resources"]["hla2_hlahd"]
    shell:
        # Same workdir-wipe rationale as _mhc_hammer_hlahd above (a partial-failure run can leave
        # a non-empty flattened directory that breaks a rerun's mv-flatten step).
        #
        # Unlike the class I rule, there's no upfront "do the expected *.est.txt files exist"
        # gate before calling the parser -- parse_hlahd_output.R already tolerates missing/
        # "No candidate." est.txt files per gene (marking that gene "not typed"), so it's always
        # safe to call unconditionally. The final content-validation check (below) still catches
        # the same class of silent, cascading failure the class I rule guards against: if HLA-HD
        # itself never ran or typed nothing at all, every gene comes back "not typed" and this
        # fails loudly instead of producing a silently-useless all-"not typed" CSV. Partial
        # results (e.g. DRB3/DRB4/DRB5 legitimately "not typed" for a patient whose haplotypes
        # don't carry that gene) are expected and not treated as failures.
        #
        # Also unlike the class I rule: no `cd {params.workdir}` before invoking the parser.
        # hlahd_parse_output.R (upstream's own script, not ours to change) writes its output to a
        # relative path internally and needs a `cd` + a follow-up `mv` step; parse_hlahd_output.R
        # is ours, so it just takes the real --hlahd_folder/--output paths directly (Snakemake's
        # own {output.X} references, already correctly resolved regardless of shell CWD) -- no
        # relative-path indirection needed at all. Confirmed on a real run (2026-08-11) that
        # copying the class I rule's `cd` pattern here was unnecessary and actively broke a
        # relative script path (see the long comment on HLA2_PARSE_SCRIPT above).
        op.as_one_line("""
        rm -rf {params.workdir} &&
        mkdir -p {params.workdir} &&
        awk '{params.gene_awk_filter}' {params.hlahd_dir}/HLA_gene.split.txt
            > {params.workdir}/hla_class_ii_genes.txt &&
        (
        export PATH=${{PATH}}:{params.hlahd_dir}/bin &&
        bash {params.hlahd_dir}/bin/hlahd.sh -m {params.min_read_length} -c 1.0 -t {threads}
        -f {params.hlahd_dir}/freq_data
        {input.fq1} {input.fq2}
        {params.workdir}/hla_class_ii_genes.txt
        {params.hlahd_dir}/dictionary
        {wildcards.patient_id} {params.workdir} &&
        if [ "{params.keep_intermediates}" = "false" ]; then
            rm -rf {params.workdir}/{wildcards.patient_id}/exon {params.workdir}/{wildcards.patient_id}/intron
                   {params.workdir}/{wildcards.patient_id}/mapfile {params.workdir}/{wildcards.patient_id}/maplist;
        fi &&
        rm -f {params.workdir}/{wildcards.patient_id}/pickup.sh {params.workdir}/{wildcards.patient_id}/estimation.sh &&
        mv {params.workdir}/{wildcards.patient_id}/* {params.workdir}/ &&
        rmdir {params.workdir}/{wildcards.patient_id}
        ) > {log.stdout} 2>&1 &&
        Rscript {params.parse_script}
        --hlahd_folder {output.result_dir}
        --sample_id {wildcards.patient_id}
        --genes {params.genes}
        --output {output.hla_alleles} >> {log.stdout} 2>&1 &&
        if ! awk -F',' '$2 != "not typed" && $3 != "not typed"' {output.hla_alleles} | grep -q .; then
            echo "ERROR: HLA-HD ran for patient {wildcards.patient_id} but failed to confidently type any HLA class II allele pair -- every gene is 'not typed' in {output.hla_alleles}. See {log.stdout}." | tee -a {log.stdout} >&2 &&
            exit 1;
        fi
        """)


# Converts the subsetted BAM to paired FASTQs (samtools collate/fastq), mirroring upstream's
# GENERATE_HLA_FQS process. These are the reads aligned against the personalised reference.
rule _mhc_hammer_generate_fqs:
    input:
        bam = str(rules._mhc_hammer_subset_bam.output.bam)
    output:
        fq1 = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.1.fq.gz",
        fq2 = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.2.fq.gz"
    log:
        stdout = CFG["logs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/generate_fqs.log"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    threads:
        CFG["threads"]["generate_fqs"]
    resources:
        **CFG["resources"]["generate_fqs"]
    shell:
        # -@ {threads}: this rule reserves {threads} cores via Snakemake but neither samtools
        # subcommand used it by default -- both collate and fastq support -@.
        op.as_one_line("""
        (
        samtools collate -@ {threads} -u -O {input.bam} |
        samtools fastq -@ {threads} -1 {output.fq1} -2 {output.fq2} -s /dev/null -0 /dev/null -n
        ) > {log.stdout} 2>&1
        """)


# Runs the user-supplied HLA-HD install on the patient's one germline WES sample to type
# HLA-A/B/C. Patient-level (once per patient, not once per sample). Mirrors upstream's
# HLAHD_LOCAL process. Fails loudly (rather than gracefully excluding the patient, as upstream
# does) if A/B/C are not all typed -- see README v1 simplifications.
#
# Restricts HLA-HD's own gene list to A/B/C (v1 is class-I-only scope; hlahd_parse_output.R is
# already only called with --genes A B C) rather than the ~25-gene stock HLA_gene.split.txt
# (full class I + II + pseudogenes), observed on a real run typing every gene one at a time at
# real (if not individually huge) per-gene cost. The restricted list is derived at runtime from
# the user's own installed HLA_gene.split.txt (grep/awk on their file) rather than bundled here,
# for the same reason novoalign_dir/hlahd_dir/mhc_hammer_scripts_dir aren't bundled -- this file
# is part of HLA-HD's own distribution, not ours to redistribute. Each gene is typed
# independently against its own IMGT/HLA reference alleles, so restricting the list doesn't
# affect A/B/C typing accuracy -- this is what the tool's own gene-list argument is for.
rule _mhc_hammer_hlahd:
    input:
        unpack(_mhc_hammer_get_germline_fqs),
        gtf = str(rules._mhc_hammer_download_reference.output.gtf)
    output:
        hla_alleles = CFG["dirs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}_hla_alleles.csv",
        result_dir = directory(CFG["dirs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/result"),
        # HLA-HD's own consolidated result file, in its native standard IMGT/HLA allele naming
        # (unlike the mhc_hammer-internal personalised-reference contig names used elsewhere in
        # this module's own outputs, e.g. DNA analysis's "hla_a_01_01_01_01"-style allele
        # identifiers) -- HLA-HD already writes this as part of the shell block above, just not
        # previously tracked as its own output. See _mhc_hammer_output_hla_final_result below.
        hla_final_result = CFG["dirs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/result/{patient_id}_final.result.txt"
    log:
        stdout = CFG["logs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/hlahd.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        hlahd_dir = HLAHD_DIR,
        workdir = lambda wildcards, output: os.path.dirname(output.hla_alleles),
        gtf_abs = lambda wildcards, input: os.path.abspath(input.gtf),
        keep_intermediates = str(CFG["options"]["keep_hlahd_intermediates"]).lower(),
        min_read_length = CFG["options"]["hlahd_min_read_length"]
    conda:
        CFG["conda_envs"]["mhc_hammer_hlahd"]
    container:
        None
    threads:
        CFG["threads"]["hlahd"]
    resources:
        **CFG["resources"]["hlahd"]
    shell:
        # Wipes the *entire* workdir up front, not just the nested {patient_id}/{patient_id}
        # HLA-HD-internal dir: an earlier partial-failure run can succeed at the mv-flatten step
        # (leaving e.g. log/ directly in {workdir}) before failing later, and mv refuses to merge
        # a freshly-flattened directory into a same-named non-empty one on a rerun ("Directory
        # not empty"). This workdir is exclusively owned by this rule for this patient, so a full
        # wipe on every invocation is safe.
        op.as_one_line("""
        rm -rf {params.workdir} &&
        mkdir -p {params.workdir} &&
        awk '$1 == "A" || $1 == "B" || $1 == "C"' {params.hlahd_dir}/HLA_gene.split.txt
            > {params.workdir}/hla_class_i_genes.txt &&
        (
        export PATH=${{PATH}}:{params.hlahd_dir}/bin &&
        bash {params.hlahd_dir}/bin/hlahd.sh -m {params.min_read_length} -c 1.0 -t {threads}
        -f {params.hlahd_dir}/freq_data
        {input.fq1} {input.fq2}
        {params.workdir}/hla_class_i_genes.txt
        {params.hlahd_dir}/dictionary
        {wildcards.patient_id} {params.workdir} &&
        if [ "{params.keep_intermediates}" = "false" ]; then
            rm -rf {params.workdir}/{wildcards.patient_id}/exon {params.workdir}/{wildcards.patient_id}/intron
                   {params.workdir}/{wildcards.patient_id}/mapfile {params.workdir}/{wildcards.patient_id}/maplist;
        fi &&
        rm -f {params.workdir}/{wildcards.patient_id}/pickup.sh {params.workdir}/{wildcards.patient_id}/estimation.sh &&
        mv {params.workdir}/{wildcards.patient_id}/* {params.workdir}/ &&
        rmdir {params.workdir}/{wildcards.patient_id}
        ) > {log.stdout} 2>&1 &&
        if [ -f {output.result_dir}/{wildcards.patient_id}_A.est.txt ] &&
           [ -f {output.result_dir}/{wildcards.patient_id}_B.est.txt ] &&
           [ -f {output.result_dir}/{wildcards.patient_id}_C.est.txt ]; then
            (cd {params.workdir} &&
             Rscript {params.scripts_dir}/bin/hlahd_parse_output.R
             --hlahd_folder result --gtf_path {params.gtf_abs}
             --sample_id {wildcards.patient_id} --genes A B C &&
             mv result/{wildcards.patient_id}_hla_alleles.csv {wildcards.patient_id}_hla_alleles.csv) >> {log.stdout} 2>&1;
            if ! awk -F',' '$2 != "not typed" && $3 != "not typed"' {params.workdir}/{wildcards.patient_id}_hla_alleles.csv | grep -q .; then
                echo "ERROR: HLA-HD produced A/B/C estimate files for patient {wildcards.patient_id} but failed to confidently type any allele pair -- every gene is 'not typed' in {wildcards.patient_id}_hla_alleles.csv. Downstream scripts (make_bed_file.R) cannot handle this and would crash with a much more confusing error. See {log.stdout}." | tee -a {log.stdout} >&2 &&
                exit 1;
            fi;
        else
            echo "ERROR: HLA-HD failed to produce HLA A, B and C estimates for patient {wildcards.patient_id}. See {log.stdout}." | tee -a {log.stdout} >&2 &&
            exit 1;
        fi
        """)


# Generates the patient's personalised HLA reference (genome FASTA/GTF/novoindex/SNP positions)
# from their HLA-HD genotype. Patient-level. DNA-only scope -- skips the transcriptome reference
# upstream also builds for RNA analysis. Mirrors upstream's GENERATE_REFERENCES process.
#
# Output is a single directory containing both the fixed-name files referenced explicitly by
# downstream rules (below) and the per-allele *_genome.snp_pos.bed files whose names depend on
# the runtime HLA genotype (see design note near the top of this file).
#   Fixed-name contents: {patient_id}.gtf, {patient_id}_mhc_genome_reference.fa(.fai),
#   {patient_id}_mhc_genome_reference.dict, {patient_id}.mhc_genome_fasta.nix,
#   {patient_id}_sorted.gtf.gz(.tbi), {patient_id}_genome_allele_table.csv, {patient_id}.bed
rule _mhc_hammer_generate_references:
    input:
        genotype = str(rules._mhc_hammer_hlahd.output.hla_alleles),
        mhc_gtf = str(rules._mhc_hammer_download_reference.output.gtf),
        mhc_fasta = str(rules._mhc_hammer_download_reference.output.genome_fasta)
    output:
        patient_dir = directory(CFG["dirs"]["patient_reference"] + "{seq_type}--{genome_build}/{patient_id}")
    log:
        stdout = CFG["logs"]["patient_reference"] + "{seq_type}--{genome_build}/{patient_id}/generate_references.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        novoalign_dir = NOVOALIGN_DIR,
        genotype_abs = lambda wildcards, input: os.path.abspath(input.genotype),
        mhc_gtf_abs = lambda wildcards, input: os.path.abspath(input.mhc_gtf),
        mhc_fasta_abs = lambda wildcards, input: os.path.abspath(input.mhc_fasta)
    conda:
        CFG["conda_envs"]["mhc_hammer_reference"]
    container:
        None
    threads:
        CFG["threads"]["generate_references"]
    resources:
        **CFG["resources"]["generate_references"]
    shell:
        op.as_one_line("""
        mkdir -p {output.patient_dir} &&
        (
        cd {output.patient_dir} &&
        export PATH=${{PATH}}:{params.novoalign_dir} &&
        Rscript {params.scripts_dir}/bin/make_bed_file.R
          --patient_hla_alleles {params.genotype_abs} --patient_id {wildcards.patient_id} --mhc_gtf {params.mhc_gtf_abs} &&
        Rscript {params.scripts_dir}/bin/allele_mismatch.R
          --genome_or_transcriptome genome --patient_id {wildcards.patient_id}
          --patient_hla_alleles {params.genotype_abs} --mhc_fasta {params.mhc_fasta_abs} &&
        Rscript {params.scripts_dir}/bin/create_gtf_and_ref.R
          --gtf_out_path {wildcards.patient_id}.gtf --fa_out_path {wildcards.patient_id}_mhc_genome_reference.fa
          --genome_out_path {wildcards.patient_id}_genome_size.txt --hla_path {params.genotype_abs}
          --mhc_fasta {params.mhc_fasta_abs} --mhc_gtf {params.mhc_gtf_abs} &&
        sort -k1,1 -k4,4n -k5,5n -t$'\\t' {wildcards.patient_id}.gtf > {wildcards.patient_id}_sorted.gtf &&
        bgzip {wildcards.patient_id}_sorted.gtf &&
        tabix -p gff {wildcards.patient_id}_sorted.gtf.gz &&
        samtools faidx {wildcards.patient_id}_mhc_genome_reference.fa &&
        picard CreateSequenceDictionary -R {wildcards.patient_id}_mhc_genome_reference.fa &&
        {params.novoalign_dir}/novoindex {wildcards.patient_id}.mhc_genome_fasta.nix {wildcards.patient_id}_mhc_genome_reference.fa
        ) > {log.stdout} 2>&1
        """)


# Aligns a sample's HLA-region reads against its patient's personalised reference with the
# user-supplied Novoalign. Sample-level. Deliberately split from the sort/filter/index
# postprocessing below (unlike upstream's single NOVOALIGN process, which bundles both) so each
# half of the work can request the threads it can actually use: unlicensed Novoalign itself is
# hard-capped at 1 thread, but samtools sort genuinely benefits from more -- bundled into one
# rule, the whole job would be stuck reserving only 1 thread (or over-reserving cores that sit
# idle for however long novoalign itself runs) rather than sizing each step independently.
rule _mhc_hammer_novoalign:
    input:
        unpack(lambda wildcards: {
            "fq1": str(rules._mhc_hammer_generate_fqs.output.fq1),
            "fq2": str(rules._mhc_hammer_generate_fqs.output.fq2)
        }),
        patient_dir = _mhc_hammer_reference_dir_for_sample
    output:
        bam = temp(CFG["dirs"]["novoalign"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.raw.bam")
    log:
        stderr = CFG["logs"]["novoalign"] + "{seq_type}--{genome_build}/{sample_id}/novoalign.stderr.log"
    params:
        novoalign_dir = NOVOALIGN_DIR,
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id),
        novoindex = lambda wildcards, input: f"{input.patient_dir[0]}/{_mhc_hammer_get_patient_id_for_sample(wildcards.sample_id)}.mhc_genome_fasta.nix"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        None
    threads:
        CFG["threads"]["novoalign"]
    resources:
        **CFG["resources"]["novoalign"]
    shell:
        # `A && B && C &` backgrounds the *entire* A&&B&&C chain as one subshell (a bare `&`
        # terminates whatever AND-OR list precedes it) -- an earlier version of this rule had
        # `mkfifo ... &&` directly followed by the first backgrounded `gzip ... &`, which swept
        # export/wd=/mkfifo into that same backgrounded subshell, so `wd` never propagated back to
        # the foreground novoalign/samtools commands ("wd: unbound variable", confirmed on a real
        # run). Fixed by terminating the foreground setup chain with `;` before the two gzip jobs
        # (which must stay backgrounded -- novoalign needs to start reading both FIFOs concurrently
        # while gzip is still writing, or the FIFOs deadlock). Also self-cleans stale FIFOs left by
        # any earlier partial-failure run of this rule (mkfifo refuses to overwrite an existing
        # file, unlike `>` redirection).
        op.as_one_line("""
        export PATH=${{PATH}}:{params.novoalign_dir} &&
        wd=$(dirname {output.bam}) &&
        rm -f $wd/fq1_uncompressed $wd/fq2_uncompressed &&
        mkfifo $wd/fq1_uncompressed $wd/fq2_uncompressed ;
        gzip -cdf {input.fq1} > $wd/fq1_uncompressed &
        gzip -cdf {input.fq2} > $wd/fq2_uncompressed &
        {params.novoalign_dir}/novoalign -d {params.novoindex} -f $wd/fq1_uncompressed $wd/fq2_uncompressed
          -F STDFQ -R 0 -r All 9999 -o SAM -o FullNW 1> $wd/{wildcards.sample_id}.sam 2> {log.stderr} &&
        samtools view -b -o {output.bam} $wd/{wildcards.sample_id}.sam &&
        rm $wd/fq1_uncompressed $wd/fq2_uncompressed $wd/{wildcards.sample_id}.sam
        """)


# Sorts, filters to properly-paired reads, reheaders, and indexes the raw Novoalign alignment --
# split out from _mhc_hammer_novoalign above specifically so this can request more threads than
# Novoalign itself is able to use. Mirrors the tail end of upstream's single NOVOALIGN process.
rule _mhc_hammer_novoalign_postprocess:
    input:
        bam = str(rules._mhc_hammer_novoalign.output.bam)
    output:
        bam = CFG["dirs"]["novoalign"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.hla.rehead.bam",
        bai = CFG["dirs"]["novoalign"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.hla.rehead.bam.bai"
    log:
        stdout = CFG["logs"]["novoalign"] + "{seq_type}--{genome_build}/{sample_id}/novoalign_postprocess.log"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    threads:
        CFG["threads"]["novoalign_postprocess"]
    resources:
        **CFG["resources"]["novoalign_postprocess"]
    shell:
        op.as_one_line("""
        (
        wd=$(dirname {output.bam}) &&
        samtools sort -@ {threads} -o $wd/{wildcards.sample_id}.sorted.bam {input.bam} &&
        samtools view -@ {threads} -f 2 -b -o $wd/{wildcards.sample_id}.hla.bam $wd/{wildcards.sample_id}.sorted.bam &&
        samtools addreplacerg -r ID:{wildcards.sample_id} -r SM:{wildcards.sample_id} -o {output.bam} $wd/{wildcards.sample_id}.hla.bam &&
        samtools index -@ {threads} {output.bam} &&
        rm $wd/{wildcards.sample_id}.sorted.bam $wd/{wildcards.sample_id}.hla.bam
        ) > {log.stdout} 2>&1
        """)


# Splits the Novoalign BAM into per-allele filtered BAMs, and records which HLA genes/alleles
# passed filtering (heterozygous or not) via `make_hla_bams.sh`. Sample-level. Mirrors upstream's
# MAKE_HLA_ALLELE_BAMS process.
#
# The 4 outputs below are the fixed-name files this module relies on; the per-allele
# *.sorted.filtered.bam files the script also writes into the same directory are not
# individually tracked by Snakemake (see design note near the top of this file) -- downstream
# rules discover them via the personalised reference FASTA headers, exactly as upstream does.
rule _mhc_hammer_make_allele_bams:
    input:
        bam = str(rules._mhc_hammer_novoalign_postprocess.output.bam),
        patient_dir = _mhc_hammer_reference_dir_for_sample
    output:
        passed_hla_genes = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_passed_hla_genes.txt",
        passed_heterozygous_hla_genes = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_passed_heterozygous_hla_genes.txt",
        passed_hla_alleles = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_passed_hla_alleles.txt",
        passed_heterozygous_hla_alleles = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_passed_heterozygous_hla_alleles.txt",
        hla_bam_read_count = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_" + MHC_SEQ + "_novoalign.hla_bam_read_count.csv"
    log:
        stdout = CFG["logs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/make_allele_bams.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        max_mismatch = CFG["options"]["max_mismatch"],
        mhc_seq = MHC_SEQ,
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id),
        bam_abs = lambda wildcards, input: os.path.abspath(input.bam),
        fasta_abs = lambda wildcards, input: os.path.abspath(
            f"{input.patient_dir[0]}/{_mhc_hammer_get_patient_id_for_sample(wildcards.sample_id)}_mhc_genome_reference.fa"
        )
    conda:
        CFG["conda_envs"]["mhc_hammer_alignment"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["make_allele_bams"]
    resources:
        **CFG["resources"]["make_allele_bams"]
    shell:
        # make_hla_bams.sh's own internal `grep ^'>' *reference.fa` (upstream source, not
        # modifiable -- see licensing note) is a bare CWD glob with no corresponding CLI argument
        # for it. In the original Nextflow pipeline this "just works" because Nextflow
        # auto-stages every declared process input into the task's isolated working directory,
        # implicitly placing the reference FASTA there; Snakemake has no equivalent auto-staging,
        # so it's symlinked into $wd explicitly here (confirmed by reading the real upstream
        # script + its calling Nextflow process, modules/local/alignment_modules.nf).
        #
        # make_hla_bams.sh itself only produces the per-allele BAMs + hla_bam_read_count.csv --
        # the passed_hla_genes/passed_heterozygous_*/passed_hla_alleles outputs are generated by
        # additional inline shell logic in upstream's own MAKE_HLA_ALLELE_BAMS Nextflow process
        # script block, not by any separate script file. That's orchestration code (like this
        # rule itself), not one of the license-restricted bin/*.sh scripts, so it's ported/adapted
        # here rather than called externally. Ported behaviour verbatim, including a known
        # upstream inconsistency: the allele_count==1 branch's BAM filename omits the
        # "_novoalign" aligner component that make_hla_bams.sh's own output actually includes, so
        # a homozygous/hemizygous gene's single-allele BAM never matches there and that gene is
        # never added to passed_hla_genes.txt via that branch. Preserved as-is to match upstream's
        # actual (if buggy) behaviour rather than silently diverging from it.
        op.as_one_line("""
        wd=$(dirname {output.passed_hla_genes}) &&
        mkdir -p $wd &&
        ln -sf {params.fasta_abs} $wd/{params.patient_id}_mhc_genome_reference.fa &&
        (cd $wd &&
         {params.scripts_dir}/bin/make_hla_bams.sh {params.bam_abs} {params.scripts_dir}/bin/
           {params.max_mismatch} {wildcards.sample_id} {params.mhc_seq} novoalign &&
         touch {wildcards.sample_id}_passed_heterozygous_hla_genes.txt
               {wildcards.sample_id}_passed_heterozygous_hla_alleles.txt
               {wildcards.sample_id}_passed_hla_genes.txt
               {wildcards.sample_id}_passed_hla_alleles.txt &&
         hla_genes=$(grep '^>' {params.fasta_abs} | sed 's/^>//' | cut -d '_' -f 1,2 | sort -u) &&
         for hla_gene in $hla_genes; do
             alleles=($(grep '^>' {params.fasta_abs} | sed 's/^>//' | grep $hla_gene)) ;
             allele_count=${{#alleles[@]}} ;
             if [ $allele_count -eq "2" ]; then
                 allele1_bam={wildcards.sample_id}_{params.mhc_seq}_novoalign.${{alleles[0]}}.sorted.filtered.bam ;
                 allele2_bam={wildcards.sample_id}_{params.mhc_seq}_novoalign.${{alleles[1]}}.sorted.filtered.bam ;
                 if [ -f $allele1_bam ] && [ -f $allele2_bam ]; then
                     echo $hla_gene >> {wildcards.sample_id}_passed_heterozygous_hla_genes.txt ;
                     echo $hla_gene >> {wildcards.sample_id}_passed_hla_genes.txt ;
                 fi ;
                 if [ -f $allele1_bam ]; then
                     echo ${{alleles[0]}} >> {wildcards.sample_id}_passed_heterozygous_hla_alleles.txt ;
                 fi ;
                 if [ -f $allele2_bam ]; then
                     echo ${{alleles[1]}} >> {wildcards.sample_id}_passed_heterozygous_hla_alleles.txt ;
                 fi ;
             elif [ $allele_count -eq "1" ]; then
                 allele1_bam={wildcards.sample_id}_{params.mhc_seq}.${{alleles[0]}}.sorted.filtered.bam ;
                 if [ -f $allele1_bam ]; then
                     echo $hla_gene >> {wildcards.sample_id}_passed_hla_genes.txt ;
                 fi ;
             fi ;
         done &&
         alleles_all=$(grep '^>' {params.fasta_abs} | sed 's/^>//' | sort -u) &&
         for hla_allele in $alleles_all; do
             bam_path={wildcards.sample_id}_{params.mhc_seq}_novoalign.${{hla_allele}}.sorted.filtered.bam ;
             if [ -f $bam_path ]; then
                 echo $hla_allele >> {wildcards.sample_id}_passed_hla_alleles.txt ;
             fi ;
         done) > {log.stdout} 2>&1
        """)


# Computes per-allele coverage with mosdepth against the patient's personalised reference BED.
# Sample-level. Mirrors upstream's MOSDEPTH process.
rule _mhc_hammer_mosdepth:
    input:
        allele_bams_marker = str(rules._mhc_hammer_make_allele_bams.output.passed_hla_alleles),
        patient_dir = _mhc_hammer_reference_dir_for_sample
    output:
        bed = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}." + MHC_SEQ + ".novoalign.mosdepth.bed"
    log:
        stdout = CFG["logs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/mosdepth.log"
    params:
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id),
        bed_file = lambda wildcards, input: f"{input.patient_dir[0]}/{_mhc_hammer_get_patient_id_for_sample(wildcards.sample_id)}.bed",
        mhc_seq = MHC_SEQ
    conda:
        CFG["conda_envs"]["mosdepth"]
    container:
        CFG["container_envs"]["mosdepth"]
    threads:
        CFG["threads"]["mosdepth"]
    resources:
        **CFG["resources"]["mosdepth"]
    shell:
        op.as_one_line("""
        (
        wd=$(dirname {input.allele_bams_marker}) &&
        > {output.bed} &&
        for bam_file in $wd/*_{params.mhc_seq}_novoalign.*.sorted.filtered.bam; do
            [ -e "$bam_file" ] || continue;
            allele=$(basename $bam_file | sed -E 's/^.*_novoalign\\.//; s/\\.sorted\\.filtered\\.bam$//');
            mosdepth --no-per-base --chrom $allele --by {params.bed_file} $wd/{wildcards.sample_id}.$allele $bam_file &&
            gunzip -f $wd/{wildcards.sample_id}.$allele.regions.bed.gz &&
            cat $wd/{wildcards.sample_id}.$allele.regions.bed >> {output.bed};
        done
        ) > {log.stdout} 2>&1
        """)


##### PAIR-LEVEL HELPER FUNCTIONS #####


def _mhc_hammer_pair_inputs(wildcards, second_allele_marker):
    # Shared plumbing for the two pair-level rules below (detect_cn_aib, detect_muts): both need
    # the tumour+normal allele-bam directories, the patient reference, and (to compute the
    # gene/allele intersection) one marker file from each sample.
    patient_dir = expand(
        str(rules._mhc_hammer_generate_references.output.patient_dir),
        patient_id = _mhc_hammer_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type),
        allow_missing = True
    )
    return {
        "tumour_marker": expand(second_allele_marker, sample_id = wildcards.tumour_id, allow_missing = True),
        "normal_marker": expand(second_allele_marker, sample_id = wildcards.normal_id, allow_missing = True),
        # get_cn.R/get_expected_depth.R/get_logr_aib.R's --tumour_library_size_path/
        # --gl_library_size_path expect a file containing a single number
        # (as.numeric(readLines(path))) -- that's _mhc_hammer_flagstat's own "library_size"
        # output (a bare read-pair count), not its "flagstat" output (samtools flagstat's
        # multi-line human-readable summary). Passing the wrong one here previously made
        # as.numeric() return a vector of NAs, which made these scripts' own
        # `if (is.na(...) | is.na(...))` check error under R >= 4.2's stricter length>1
        # condition handling -- confirmed on a real run (2026-07-31).
        "tumour_library_size": expand(str(rules._mhc_hammer_flagstat.output.library_size), sample_id = wildcards.tumour_id, allow_missing = True),
        "normal_library_size": expand(str(rules._mhc_hammer_flagstat.output.library_size), sample_id = wildcards.normal_id, allow_missing = True),
        "patient_dir": patient_dir
    }


# Detects HLA allele-specific copy number and allelic imbalance (LOH) by comparing tumour and
# germline coverage at germline-heterozygous SNP positions within each HLA gene the tumour and
# normal both passed filtering for. Pair-level. Mirrors upstream's DETECT_CN_AND_AIB process
# ("all_snps" variant only -- upstream's "exon_snps" variant is disabled/commented out upstream
# too). Copy-number output is skipped (but allelic-imbalance/LOH detection still runs) for
# pairs missing a Battenberg `cellularity_ploidy.txt` file, matching upstream's own
# `missing_purity_ploidy` flag.
rule _mhc_hammer_detect_cn_aib:
    input:
        unpack(lambda wildcards: _mhc_hammer_pair_inputs(wildcards, str(rules._mhc_hammer_make_allele_bams.output.passed_heterozygous_hla_genes))),
        cellularity_ploidy = _mhc_hammer_get_battenberg_cp_input
    output:
        # Filename must match concatenate_dna_analysis_tables.R's own
        # paste0(sample_name,"_",snp_type,"_",aligner,"_dna_analysis.csv") convention exactly --
        # confirmed by reading that script.
        dna_analysis = CFG["dirs"]["dna_analysis"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_all_snps_novoalign_dna_analysis.csv"
    log:
        stdout = CFG["logs"]["dna_analysis"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/detect_cn_aib.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        # Per-seq_type -- see the long comment on options.min_depth in default.yaml for why (WGS
        # depth at this locus is structurally much lower than WES's on-target-enriched depth).
        # References the MIN_DEPTH module-level constant, not CFG directly (CFG-closure gotcha).
        min_depth = lambda wildcards: MIN_DEPTH[wildcards.seq_type],
        mhc_seq = MHC_SEQ,
        purity_ploidy = _mhc_hammer_parse_cellularity_ploidy,
        tumour_marker_abs = lambda wildcards, input: os.path.abspath(input.tumour_marker[0]),
        normal_marker_abs = lambda wildcards, input: os.path.abspath(input.normal_marker[0]),
        tumour_library_size_abs = lambda wildcards, input: os.path.abspath(input.tumour_library_size[0]),
        normal_library_size_abs = lambda wildcards, input: os.path.abspath(input.normal_library_size[0]),
        patient_dir_abs = lambda wildcards, input: os.path.abspath(input.patient_dir[0]),
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type)
    conda:
        CFG["conda_envs"]["mhc_hammer_cn_aib"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["detect_cn_aib"]
    resources:
        **CFG["resources"]["detect_cn_aib"]
    shell:
        op.as_one_line("""
        wd=$(dirname {output.dna_analysis}) && mkdir -p $wd &&
        (
        cd $wd &&
        fasta={params.patient_dir_abs}/{params.patient_id}_mhc_genome_reference.fa &&
        tumour_dir=$(dirname {params.tumour_marker_abs}) &&
        normal_dir=$(dirname {params.normal_marker_abs}) &&
        purity="{params.purity_ploidy[0]}" && ploidy="{params.purity_ploidy[1]}" &&
        genes=$(comm -12 <(sort {params.tumour_marker_abs}) <(sort {params.normal_marker_abs})) &&
        for gene in $genes; do
            alleles=($(grep '^>' $fasta | sed 's/^>//' | grep $gene | sort -u)) &&
            allele1=${{alleles[0]}} && allele2=${{alleles[1]}} &&
            allele1_snp_path=${{allele1}}_genome.snp_pos.bed && allele2_snp_path=${{allele2}}_genome.snp_pos.bed &&
            cp {params.patient_dir_abs}/$allele1_snp_path . 2>/dev/null;
            cp {params.patient_dir_abs}/$allele2_snp_path . 2>/dev/null;
            allele1_tumour_bam=$tumour_dir/{wildcards.tumour_id}_{params.mhc_seq}_novoalign.$allele1.sorted.filtered.bam &&
            allele2_tumour_bam=$tumour_dir/{wildcards.tumour_id}_{params.mhc_seq}_novoalign.$allele2.sorted.filtered.bam &&
            allele1_gl_bam=$normal_dir/{wildcards.normal_id}_{params.mhc_seq}_novoalign.$allele1.sorted.filtered.bam &&
            allele2_gl_bam=$normal_dir/{wildcards.normal_id}_{params.mhc_seq}_novoalign.$allele2.sorted.filtered.bam &&
            samtools mpileup $allele1_tumour_bam -f $fasta > {wildcards.tumour_id}.$allele1.mpileup;
            samtools mpileup $allele1_gl_bam -f $fasta > {wildcards.normal_id}.$allele1.mpileup;
            samtools mpileup $allele2_tumour_bam -f $fasta > {wildcards.tumour_id}.$allele2.mpileup;
            samtools mpileup $allele2_gl_bam -f $fasta > {wildcards.normal_id}.$allele2.mpileup;
            Rscript {params.scripts_dir}/bin/mpileup_to_csv.R --mpileup_path {wildcards.tumour_id}.$allele1.mpileup --hla_ref_path $fasta --allele $allele1 --out_csv_name {wildcards.tumour_id}.$allele1.coverage.csv;
            Rscript {params.scripts_dir}/bin/mpileup_to_csv.R --mpileup_path {wildcards.normal_id}.$allele1.mpileup --hla_ref_path $fasta --allele $allele1 --out_csv_name {wildcards.normal_id}.$allele1.coverage.csv;
            Rscript {params.scripts_dir}/bin/mpileup_to_csv.R --mpileup_path {wildcards.tumour_id}.$allele2.mpileup --hla_ref_path $fasta --allele $allele2 --out_csv_name {wildcards.tumour_id}.$allele2.coverage.csv;
            Rscript {params.scripts_dir}/bin/mpileup_to_csv.R --mpileup_path {wildcards.normal_id}.$allele2.mpileup --hla_ref_path $fasta --allele $allele2 --out_csv_name {wildcards.normal_id}.$allele2.coverage.csv;
            Rscript {params.scripts_dir}/bin/get_filtered_pos_bed.R --coverage_path {wildcards.normal_id}.$allele1.coverage.csv --min_depth {params.min_depth} --out_bed_name {wildcards.normal_id}.$allele1.filtered_positions.bed --allele $allele1;
            Rscript {params.scripts_dir}/bin/get_filtered_pos_bed.R --coverage_path {wildcards.normal_id}.$allele2.coverage.csv --min_depth {params.min_depth} --out_bed_name {wildcards.normal_id}.$allele2.filtered_positions.bed --allele $allele2;
            Rscript {params.scripts_dir}/bin/get_filtered_coverage.R --coverage_path {wildcards.tumour_id}.$allele1.coverage.csv --filtered_positions_path {wildcards.normal_id}.$allele1.filtered_positions.bed --out_csv {wildcards.tumour_id}.$allele1.coverage_at_filtered_positions.csv;
            Rscript {params.scripts_dir}/bin/get_filtered_coverage.R --coverage_path {wildcards.tumour_id}.$allele2.coverage.csv --filtered_positions_path {wildcards.normal_id}.$allele2.filtered_positions.bed --out_csv {wildcards.tumour_id}.$allele2.coverage_at_filtered_positions.csv;
            Rscript {params.scripts_dir}/bin/get_filtered_coverage.R --coverage_path {wildcards.normal_id}.$allele1.coverage.csv --filtered_positions_path {wildcards.normal_id}.$allele1.filtered_positions.bed --out_csv {wildcards.normal_id}.$allele1.coverage_at_filtered_positions.csv;
            Rscript {params.scripts_dir}/bin/get_filtered_coverage.R --coverage_path {wildcards.normal_id}.$allele2.coverage.csv --filtered_positions_path {wildcards.normal_id}.$allele2.filtered_positions.bed --out_csv {wildcards.normal_id}.$allele2.coverage_at_filtered_positions.csv;
            Rscript {params.scripts_dir}/bin/get_gene_filtered_snp_positions.R --allele1_snp_bed $allele1_snp_path --allele2_snp_bed $allele2_snp_path --allele1 $allele1 --allele2 $allele2 --allele1_filtered_pos_bed {wildcards.normal_id}.$allele1.filtered_positions.bed --allele2_filtered_pos_bed {wildcards.normal_id}.$allele2.filtered_positions.bed --sample_name {wildcards.normal_id};
            allele1_gene_filtered_snps={wildcards.normal_id}.$allele1.filtered_snp_positions.bed && allele2_gene_filtered_snps={wildcards.normal_id}.$allele2.filtered_snp_positions.bed &&
            if [ -n "$purity" ] && [ -n "$ploidy" ]; then
                Rscript {params.scripts_dir}/bin/get_cn.R --allele1 $allele1 --allele2 $allele2 --allele1_gl_coverage_file {wildcards.normal_id}.$allele1.coverage_at_filtered_positions.csv --allele2_gl_coverage_file {wildcards.normal_id}.$allele2.coverage_at_filtered_positions.csv --allele1_tumour_coverage_file {wildcards.tumour_id}.$allele1.coverage_at_filtered_positions.csv --allele2_tumour_coverage_file {wildcards.tumour_id}.$allele2.coverage_at_filtered_positions.csv --allele1_snp_bed $allele1_gene_filtered_snps --allele2_snp_bed $allele2_gene_filtered_snps --purity $purity --ploidy $ploidy --tumour_library_size_path {params.tumour_library_size_abs} --gl_library_size_path {params.normal_library_size_abs} --gtf_path {params.patient_dir_abs}/{params.patient_id}.gtf --cn_output_path {wildcards.tumour_id}.$gene.all_snps.cn.csv --cn_plots_prefix {wildcards.tumour_id}.$gene.all_snps.cn --scripts_dir {params.scripts_dir}/bin/;
                Rscript {params.scripts_dir}/bin/get_expected_depth.R --allele1 $allele1 --allele2 $allele2 --allele1_gl_coverage_file {wildcards.normal_id}.$allele1.coverage_at_filtered_positions.csv --allele2_gl_coverage_file {wildcards.normal_id}.$allele2.coverage_at_filtered_positions.csv --allele1_snp_bed $allele1_gene_filtered_snps --allele2_snp_bed $allele2_gene_filtered_snps --purity $purity --tumour_library_size_path {params.tumour_library_size_abs} --gl_library_size_path {params.normal_library_size_abs} --expected_depth_output_path {wildcards.tumour_id}.$gene.all_snps.expected_depth.csv --expected_depth_plots_prefix {wildcards.tumour_id}.$gene.all_snps.expected_depth --scripts_dir {params.scripts_dir}/bin/;
                missing_purity_ploidy=FALSE;
            else
                missing_purity_ploidy=TRUE;
            fi;
            for role_allele_bam in "{wildcards.tumour_id} $allele1 $allele1_tumour_bam" "{wildcards.tumour_id} $allele2 $allele2_tumour_bam" "{wildcards.normal_id} $allele1 $allele1_gl_bam" "{wildcards.normal_id} $allele2 $allele2_gl_bam"; do
                set -- $role_allele_bam; role=$1; allele=$2; bam=$3;
                snp_bed=${{allele}}_genome.snp_pos.bed;
                bedtools intersect -loj -bed -b $bam -a $snp_bed | cut -f1-7 > $role.$allele.snp_reads_overlap_7cols.bed;
                Rscript {params.scripts_dir}/bin/count_reads_once.R --snp_reads_overlap_bed $role.$allele.snp_reads_overlap_7cols.bed --snp_path $snp_bed --out_csv $role.$allele.coverage_at_filtered_snps_reads_count_once.csv;
            done;
            Rscript {params.scripts_dir}/bin/get_logr_aib.R --allele1 $allele1 --allele2 $allele2 --allele1_gl_reads_count_once_coverage {wildcards.normal_id}.$allele1.coverage_at_filtered_snps_reads_count_once.csv --allele2_gl_reads_count_once_coverage {wildcards.normal_id}.$allele2.coverage_at_filtered_snps_reads_count_once.csv --allele1_tumour_reads_count_once_coverage {wildcards.tumour_id}.$allele1.coverage_at_filtered_snps_reads_count_once.csv --allele2_tumour_reads_count_once_coverage {wildcards.tumour_id}.$allele2.coverage_at_filtered_snps_reads_count_once.csv --allele1_snp_bed $allele1_gene_filtered_snps --allele2_snp_bed $allele2_gene_filtered_snps --tumour_library_size_path {params.tumour_library_size_abs} --gl_library_size_path {params.normal_library_size_abs} --logr_aib_output_path {wildcards.tumour_id}.$gene.all_snps.logr_aib.csv --logr_aib_plots_prefix {wildcards.tumour_id}.$gene.all_snps.logr_aib --scripts_dir {params.scripts_dir}/bin/;
        done &&
        Rscript {params.scripts_dir}/bin/concatenate_dna_analysis_tables.R --genes $genes --snp_type all_snps --sample_name {wildcards.tumour_id} --aligner novoalign --missing_purity_ploidy ${{missing_purity_ploidy:-TRUE}}
        ) > {log.stdout} 2>&1
        """)
        # `--genes` is deliberately unquoted: concatenate_dna_analysis_tables.R's `--genes`
        # uses argparse nargs="+", so it must receive each gene as a separate word-split CLI
        # token, not one quoted string. Argument names/order in this rule were cross-checked
        # against a local read-only clone of mhc-hammer's bin/*.R scripts.


# Calls somatic mutations in each HLA allele the tumour and normal both passed filtering for,
# using Mutect2 in tumour-normal mode followed by VEP annotation against the personalised
# per-patient HLA GTF/FASTA (not a standard genome-build VEP cache). Pair-level. Mirrors
# upstream's DETECT_MUTS process. VEP is user-supplied (vep_path/vep_cache), not conda-managed
# -- see licensing/VEP note near the top of this file.
#
# Output is a directory (not individually-named files) because the set of *.vep.txt files
# depends on the runtime HLA genotype and on which alleles actually carried a mutation; a
# `{tumour_id}.empty.vep.txt` sentinel is always written if none did, so the directory is never
# empty and downstream globbing (rule _mhc_hammer_parse_mutations) is always well-defined.
rule _mhc_hammer_detect_muts:
    input:
        unpack(lambda wildcards: _mhc_hammer_pair_inputs(wildcards, str(rules._mhc_hammer_make_allele_bams.output.passed_hla_alleles)))
    output:
        vep_dir = directory(CFG["dirs"]["mutations"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/vep")
    log:
        stdout = CFG["logs"]["mutations"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/detect_muts.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        vep_path = VEP_PATH,
        vep_cache = VEP_CACHE,
        tumour_marker_abs = lambda wildcards, input: os.path.abspath(input.tumour_marker[0]),
        normal_marker_abs = lambda wildcards, input: os.path.abspath(input.normal_marker[0]),
        patient_dir_abs = lambda wildcards, input: os.path.abspath(input.patient_dir[0]),
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type),
        mem_gb = lambda wildcards, resources: max(1, int(resources.mem_mb / 1000 * 0.8)),
        mhc_seq = MHC_SEQ
    conda:
        CFG["conda_envs"]["mhc_hammer_detect_muts"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir and vep_path -- see note near the top of this file
    threads:
        CFG["threads"]["detect_muts"]
    resources:
        **CFG["resources"]["detect_muts"]
    shell:
        op.as_one_line("""
        mkdir -p {output.vep_dir} &&
        (
        cd {output.vep_dir} &&
        export PATH={params.vep_path}:${{PATH}} &&
        fasta={params.patient_dir_abs}/{params.patient_id}_mhc_genome_reference.fa &&
        tumour_dir=$(dirname {params.tumour_marker_abs}) &&
        normal_dir=$(dirname {params.normal_marker_abs}) &&
        bgzip -dc {params.patient_dir_abs}/{params.patient_id}_sorted.gtf.gz > {wildcards.tumour_id}.gtf &&
        alleles=$(comm -12 <(sort {params.tumour_marker_abs}) <(sort {params.normal_marker_abs})) &&
        for allele in $alleles; do
            normal_bam=$normal_dir/{wildcards.normal_id}_{params.mhc_seq}_novoalign.$allele.sorted.filtered.bam;
            tumour_bam=$tumour_dir/{wildcards.tumour_id}_{params.mhc_seq}_novoalign.$allele.sorted.filtered.bam;
            gatk --java-options '-Xmx{params.mem_gb}g -Xms1g' Mutect2 -R $fasta -I $normal_bam -I $tumour_bam -normal {wildcards.normal_id} --f1r2-tar-gz {wildcards.tumour_id}.$allele.f1r2.tar.gz --output {wildcards.tumour_id}.$allele.vcf;
            if tail -n 1 {wildcards.tumour_id}.$allele.vcf | grep -q CHROM; then continue; fi;
            gatk --java-options '-Xmx{params.mem_gb}g -Xms1g' LearnReadOrientationModel -I {wildcards.tumour_id}.$allele.f1r2.tar.gz -O {wildcards.tumour_id}.$allele.read-orientation-model.tar.gz;
            gatk --java-options '-Xmx{params.mem_gb}g -Xms1g' FilterMutectCalls -V {wildcards.tumour_id}.$allele.vcf -R $fasta --ob-priors {wildcards.tumour_id}.$allele.read-orientation-model.tar.gz -O {wildcards.tumour_id}.$allele.filt.vcf;
            bcftools norm -m-any {wildcards.tumour_id}.$allele.filt.vcf --output {wildcards.tumour_id}.$allele.norm.filt.vcf;
            Rscript {params.scripts_dir}/bin/make_vep_gtf.R --gtf_path {wildcards.tumour_id}.gtf --gtf_vep_path {wildcards.tumour_id}.$allele.vep.gtf --allele $allele;
            bgzip {wildcards.tumour_id}.$allele.vep.gtf && tabix -p gff {wildcards.tumour_id}.$allele.vep.gtf.gz;
            vep -i {wildcards.tumour_id}.$allele.norm.filt.vcf --gtf {wildcards.tumour_id}.$allele.vep.gtf.gz --fasta $fasta --dir_cache {params.vep_cache} --vcf -o {wildcards.tumour_id}.$allele.norm.filt.vep.vcf --fields "Allele,Consequence,IMPACT,Feature_type,Feature,EXON,INTRON,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,FLAGS";
            gatk --java-options '-Xmx{params.mem_gb}g -Xms1g' VariantsToTable -V {wildcards.tumour_id}.$allele.norm.filt.vep.vcf -O {wildcards.tumour_id}.$allele.norm.filt.vep.txt -F CHROM -F POS -F REF -F ALT -F FILTER -GF AD -GF DP -F CSQ --show-filtered;
        done;
        if ! ls {wildcards.tumour_id}.*.vep.txt >/dev/null 2>&1; then touch {wildcards.tumour_id}.empty.vep.txt; fi
        ) > {log.stdout} 2>&1
        """)
        # make_vep_gtf.R's arguments (--gtf_path/--gtf_vep_path/--allele) were cross-checked
        # against a local read-only clone of mhc-hammer's bin/*.R.


# Synthesises the small cohort-level "inventory" table that several of upstream's cohort-level
# R scripts require (make_mutation_table.R, make_cohort_overview_table.R) -- upstream's own
# pipeline builds this from the user's *.csv samplesheet, which lcr-modules has no equivalent of;
# this rule derives the same columns from CFG["samples"]/CFG["paired_runs"]/Battenberg's
# cellularity_ploidy.txt instead. Confirmed against the actual column names those scripts read (`patient`, `sample_name`,
# `sample_type`, `sequencing_type`, `purity`, `ploidy`, `normal_sample_name`) by reading a local
# clone of mhc-hammer's bin/*.R. Cohort-wide, no wildcards, stdlib-only (no pandas dependency).
# `temp()` alone (kept below as defense in depth) assumes every invocation reaches a clean,
# successful completion, so Snakemake's own post-run cleanup actually deletes these files before
# the next run starts. Real-world runs of this pipeline routinely fail or get interrupted partway
# through (--keep-going), which means that cleanup never happens, and a stale inventory.csv from
# an earlier, differently-scoped invocation (a wider cohort, or an entirely different `subset=`
# patient) is left on disk. With `--rerun-triggers mtime` (this module's target GAMBL environment
# always sets this) and no real tracked file input reflecting "which patients/samples are in scope
# right now", Snakemake has no way to tell that leftover file is wrong, and silently reuses it --
# confirmed on a real run (2026-08-02): a `subset=DLBCL-RICOVER_1013` invocation reused a
# 50-minutes-stale inventory.csv containing six *other* patients and no DLBCL-RICOVER_1013 row at
# all, cascading into make_mutation_table.R's "Why are some tumours missing matched germlines?"
# (no row at all for that tumour --> no match in its own merge, not even an empty one). Force a
# genuine, mtime-triggered rebuild on *every* invocation, independent of rerun-triggers settings
# or prior-run cleanup success: touch a tiny sentinel file stamped with the current wall-clock
# time at parse time (this executes once per `snakemake` process, dry-run or not -- the same as
# the directory-creation side effects op.setup_module() already performs unconditionally) and
# declare it as this rule's real input, so its mtime is always newer than any previous run's
# inventory.csv/hlahd_germline_samples.txt, however old or however that prior run ended.
_mhc_hammer_invocation_marker = CFG["dirs"]["cohort_tables"] + ".invocation_marker"
os.makedirs(os.path.dirname(_mhc_hammer_invocation_marker), exist_ok = True)
with open(_mhc_hammer_invocation_marker, "w") as _mhc_hammer_marker_fh:
    _mhc_hammer_marker_fh.write(str(time.time()))

# Best-effort cohort table: only include pairs whose full upstream chain has actually completed
# on disk, rather than unconditionally requiring every matched pair in CFG["paired_runs"] to
# succeed. Real GAMBL runs routinely have a handful of samples fail deep in preprocessing (e.g.
# _mhc_hammer_subset_bam) for reasons unrelated to the cohort table -- under --keep-going,
# unrelated per-patient/per-pair branches of the DAG already complete independently, but
# _mhc_hammer_cohort_table's inputs used to be a single unconditional expand() over ALL matched
# pairs, so ONE failing pair anywhere in a large cohort made it impossible to ever produce a
# cohort table at all, even after every OTHER pair succeeded (real symptom, 2026-08-02: excluding
# two chronically-failing samples was the only way to get any cohort table out at all). Used by
# both _mhc_hammer_generate_inventory (so inventory.csv/hlahd_germline_samples.txt never
# reference a patient whose data isn't actually staged -- otherwise make_cohort_overview_table.R
# would expect data for a patient the staged directory doesn't have) and _mhc_hammer_cohort_table
# itself, so the two stay consistent with each other.
#
# IMPORTANT: Snakemake decides a rule's input list at DAG-build time, before any of THIS
# invocation's own jobs run -- so this only sees a pair as "completed" if it finished on some
# *earlier* invocation. On a first, from-scratch invocation this will see zero completed pairs
# (nothing has run yet) and produce an empty-but-valid cohort table; a *second* invocation of the
# same command, after --keep-going has already attempted everything and left real, stable,
# on-disk evidence of what succeeded vs failed, is required to get a real cohort table. This
# mirrors the same os.path.exists()-based optional-input pattern already used for Battenberg
# cellularity_ploidy elsewhere in this file, applied here to this module's own previously-
# generated outputs instead of an externally-produced file.
#
# A pair is included only if EVERY one of its required files exists (dna_analysis, patient_dir,
# and both samples' mosdepth/library_size/hla_bam_read_count) -- all-or-nothing per pair, so a
# partially-completed pair (e.g. mosdepth succeeded but detect_cn_aib later failed) never
# contributes inconsistent partial data to the cohort table, the same failure mode
# _mhc_hammer_generate_inventory's own matched-pair restriction (5ad560d5) was fixed for.
def _mhc_hammer_get_completed_pairs():
    CFG = config["lcr-modules"]["mhc_hammer"]
    completed = []
    for _, row in CFG["paired_runs"].iterrows():
        seq_type = row["tumour_seq_type"]
        genome_build = row["tumour_genome_build"]
        tumour_id = row["tumour_sample_id"]
        normal_id = row["normal_sample_id"]
        pair_status = row["pair_status"]
        patient_id = row["tumour_patient_id"]
        dna_analysis_path = str(rules._mhc_hammer_detect_cn_aib.output.dna_analysis).format(
            seq_type = seq_type, genome_build = genome_build,
            tumour_id = tumour_id, normal_id = normal_id, pair_status = pair_status
        )
        patient_dir_path = str(rules._mhc_hammer_generate_references.output.patient_dir).format(
            seq_type = seq_type, genome_build = genome_build, patient_id = patient_id
        )
        per_sample = {}
        for sample_id in (tumour_id, normal_id):
            per_sample[sample_id] = {
                "mosdepth": str(rules._mhc_hammer_mosdepth.output.bed).format(
                    seq_type = seq_type, genome_build = genome_build, sample_id = sample_id
                ),
                "library_size": str(rules._mhc_hammer_flagstat.output.library_size).format(
                    seq_type = seq_type, genome_build = genome_build, sample_id = sample_id
                ),
                "hla_bam_read_count": str(rules._mhc_hammer_make_allele_bams.output.hla_bam_read_count).format(
                    seq_type = seq_type, genome_build = genome_build, sample_id = sample_id
                ),
            }
        required_paths = [dna_analysis_path, patient_dir_path] + [
            path for sample_paths in per_sample.values() for path in sample_paths.values()
        ]
        if not all(os.path.exists(path) for path in required_paths):
            continue
        completed.append({
            "tumour_id": tumour_id,
            "normal_id": normal_id,
            "dna_analysis": dna_analysis_path,
            "patient_dir": patient_dir_path,
            "mosdepth": [per_sample[tumour_id]["mosdepth"], per_sample[normal_id]["mosdepth"]],
            "library_size": [per_sample[tumour_id]["library_size"], per_sample[normal_id]["library_size"]],
            "hla_bam_read_count": [per_sample[tumour_id]["hla_bam_read_count"], per_sample[normal_id]["hla_bam_read_count"]],
        })
    return completed

_mhc_hammer_completed_pairs = _mhc_hammer_get_completed_pairs()
_mhc_hammer_completed_dna_analysis = [p["dna_analysis"] for p in _mhc_hammer_completed_pairs]
_mhc_hammer_completed_patient_dirs = sorted(set(p["patient_dir"] for p in _mhc_hammer_completed_pairs))
_mhc_hammer_completed_mosdepth = sorted(set(
    path for p in _mhc_hammer_completed_pairs for path in p["mosdepth"]
))
_mhc_hammer_completed_library_size = sorted(set(
    path for p in _mhc_hammer_completed_pairs for path in p["library_size"]
))
_mhc_hammer_completed_hla_bam_read_count = sorted(set(
    path for p in _mhc_hammer_completed_pairs for path in p["hla_bam_read_count"]
))
print(
    f"INFO [mhc_hammer]: {len(_mhc_hammer_completed_pairs)} of {len(CFG['paired_runs'])} matched "
    f"tumour/normal pair(s) have a complete upstream chain on disk right now and will be included "
    f"in this invocation's cohort table (the rest, if any, need a rerun of this same command after "
    f"they finish or are excluded)."
)

rule _mhc_hammer_generate_inventory:
    input:
        cellularity_ploidy = _mhc_hammer_get_all_battenberg_cp_inputs,
        invocation_marker = _mhc_hammer_invocation_marker
    output:
        # temp() deletes these once nothing in the current run still needs them, so a clean
        # successful run always regenerates them next time even without the marker above --
        # cheap to regenerate (a quick Python script, no real computation), so paying that cost
        # on every run is preferable to silently stale cohort-wide state.
        inventory = temp(CFG["dirs"]["cohort_tables"] + "inventory.csv"),
        hlahd_germline_samples = temp(CFG["dirs"]["cohort_tables"] + "hlahd_germline_samples.txt")
    run:
        import csv
        CFG = config["lcr-modules"]["mhc_hammer"]
        samples = CFG["samples"]
        runs = CFG["paired_runs"]
        cp_pattern = CFG["options"]["battenberg_cellularity_ploidy"]
        purity_ploidy = {}
        if cp_pattern:
            for _, run_row in runs.iterrows():
                cp_path = cp_pattern.format(
                    seq_type = run_row["tumour_seq_type"],
                    genome_build = run_row["tumour_genome_build"],
                    tumour_id = run_row["tumour_sample_id"],
                    normal_id = run_row["normal_sample_id"]
                )
                if cp_path not in input.cellularity_ploidy:
                    continue
                with open(cp_path) as fh:
                    header = fh.readline().rstrip("\n").split("\t")
                    values = fh.readline().rstrip("\n").split("\t")
                cp_row = dict(zip(header, values))
                purity_ploidy[run_row["tumour_sample_id"]] = (cp_row.get("cellularity", ""), cp_row.get("ploidy", ""))
        normal_for_tumour = dict(zip(runs["tumour_sample_id"], runs["normal_sample_id"]))
        # Restrict to samples belonging to a real MATCHED pair (runs == CFG["paired_runs"],
        # already narrowed to pair_status == "matched" near the top of this file) -- deliberately
        # NOT the stricter _mhc_hammer_get_completed_pairs() (defined above, used by
        # _mhc_hammer_cohort_table) that also requires dna_analysis/mosdepth/library_size/
        # hla_bam_read_count/patient_dir to already exist on disk. inventory.csv's only real
        # consumer of this pairing data is make_mutation_table.R's tumour_sample_name ->
        # germline_sample_name lookup (invoked by _mhc_hammer_parse_mutations), which needs
        # nothing but the pairing relationship itself -- a pure CFG["paired_runs"] fact, always
        # immediately known, completely independent of whether OTHER, unrelated branches of that
        # patient's DAG (detect_cn_aib, mosdepth) have finished. Using the "completed" restriction
        # here (2026-08-02, commit `281d2430`) was wrong: any pair where mutation-calling finishes
        # before some unrelated step (e.g. mosdepth) got zero rows in inventory.csv, so
        # `inventory[sample_name %in% tumour_bam_read_counts$tumour_sample_name]` returned a
        # 0-row, *0-column* result -- data.table can't infer a column type from zero matches on an
        # otherwise-populated table the same way it does for a genuinely empty file, but the
        # *subsequent* merge (`by.x = "tumour_sample_name", by.y = "sample_name"`) still failed
        # with "Incompatible join types: x.sample_name (logical) and i.tumour_sample_name
        # (character)" -- confirmed on a real run (2026-08-04/05): a WGS pair (SP192798/SP192796)
        # whose detect_muts/VEP output had already completed, but whose mosdepth output for that
        # pair had not, had zero rows in inventory.csv even though its pairing was never in
        # question. Restricting to just "matched" (this fix) fixed the original genuinely-unpaired
        # case (commit `5ad560d5`) without needing every unrelated branch to also finish first;
        # _mhc_hammer_cohort_table's OWN inputs (below) correctly keep the stricter "completed"
        # restriction, since that rule genuinely does need dna_analysis/mosdepth/etc. to exist.
        matched_sample_ids = set(runs["tumour_sample_id"]) | set(runs["normal_sample_id"])
        os.makedirs(os.path.dirname(output.inventory), exist_ok=True)
        with open(output.inventory, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["patient", "sample_name", "sample_type", "sequencing_type", "purity", "ploidy", "normal_sample_name"])
            for _, row in samples.iterrows():
                sample_id = row["sample_id"]
                if sample_id not in matched_sample_ids:
                    continue
                tissue_status = str(row["tissue_status"]).lower()
                sample_type = "tumour" if tissue_status in ("tumour", "tumor") else "normal"
                purity, ploidy = purity_ploidy.get(sample_id, ("", ""))
                writer.writerow([row["patient_id"], sample_id, sample_type, MHC_SEQ, purity, ploidy, normal_for_tumour.get(sample_id, "")])
        germline_ids = sorted(set(runs["normal_sample_id"]))
        with open(output.hlahd_germline_samples, "w") as fh:
            for gid in germline_ids:
                fh.write(gid + "\n")


# Aggregates one patient's per-tumour VEP mutation tables (across all of that patient's runs)
# into a single per-patient mutation CSV. Patient-level. Mirrors upstream's PARSE_MUTATIONS
# process. Skipped (no output) for patients where no run produced any real (non-empty) VEP
# table -- matching upstream's own `mutations_detected_check` branch.
#
# make_mutation_table.R needs to run from a directory containing bare-basename copies of the
# vep.txt tables and the tumour/normal allele-filtered BAMs it reads (it looks up VCF sample
# columns and BAM files by bare sample_id, and reads BAMs by bare relative filename) --
# confirmed by reading its source. This mirrors the flat per-task staging directory every
# Nextflow process gets for free; here we build it explicitly via symlinks.
def _mhc_hammer_get_patient_mutation_inputs(wildcards):
    CFG = config["lcr-modules"]["mhc_hammer"]
    patient_runs = op.filter_samples(
        CFG["paired_runs"],
        tumour_patient_id = wildcards.patient_id,
        tumour_seq_type = wildcards.seq_type,
        tumour_genome_build = wildcards.genome_build
    )
    return {
        "vep_dirs": expand(
            str(rules._mhc_hammer_detect_muts.output.vep_dir),
            zip,
            tumour_id = patient_runs["tumour_sample_id"],
            normal_id = patient_runs["normal_sample_id"],
            pair_status = patient_runs["pair_status"],
            allow_missing = True
        ),
        "tumour_marker": expand(
            str(rules._mhc_hammer_make_allele_bams.output.passed_hla_alleles),
            zip, sample_id = patient_runs["tumour_sample_id"], allow_missing = True
        ),
        "normal_marker": expand(
            str(rules._mhc_hammer_make_allele_bams.output.passed_hla_alleles),
            zip, sample_id = sorted(set(patient_runs["normal_sample_id"])), allow_missing = True
        )
    }

def _mhc_hammer_get_patient_run_ids(wildcards):
    CFG = config["lcr-modules"]["mhc_hammer"]
    patient_runs = op.filter_samples(
        CFG["paired_runs"],
        tumour_patient_id = wildcards.patient_id,
        tumour_seq_type = wildcards.seq_type,
        tumour_genome_build = wildcards.genome_build
    )
    tumour_ids = " ".join(patient_runs["tumour_sample_id"].tolist())
    normal_ids = " ".join(sorted(set(patient_runs["normal_sample_id"].tolist())))
    return tumour_ids, normal_ids

rule _mhc_hammer_parse_mutations:
    input:
        unpack(_mhc_hammer_get_patient_mutation_inputs),
        inventory = str(rules._mhc_hammer_generate_inventory.output.inventory)
    output:
        mutations = CFG["dirs"]["mutations"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}_mutations.csv"
    log:
        stdout = CFG["logs"]["mutations"] + "{seq_type}--{genome_build}/{patient_id}/parse_mutations.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        mhc_seq = MHC_SEQ,
        vep_dirs_abs = lambda wildcards, input: " ".join(os.path.abspath(d) for d in input.vep_dirs),
        bam_dirs_abs = lambda wildcards, input: " ".join(sorted(set(
            os.path.dirname(os.path.abspath(f)) for f in list(input.tumour_marker) + list(input.normal_marker)
        ))),
        tumour_normal_ids = lambda wildcards: _mhc_hammer_get_patient_run_ids(wildcards),
        inventory_abs = lambda wildcards, input: os.path.abspath(input.inventory),
        mutations_abs = lambda wildcards, output: os.path.abspath(output.mutations)
    conda:
        CFG["conda_envs"]["mhc_hammer_r"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["parse_mutations"]
    resources:
        **CFG["resources"]["parse_mutations"]
    shell:
        op.as_one_line("""
        wd=$(dirname {output.mutations}) && mkdir -p $wd &&
        stage=$wd/staged && rm -rf $stage && mkdir -p $stage &&
        (
        for d in {params.vep_dirs_abs}; do
            for f in $d/*.vep.txt; do
                case "$f" in *empty.vep.txt) continue ;; esac;
                [ -e "$f" ] && ln -sf $(realpath $f) $stage/$(basename $f);
            done;
        done;
        for d in {params.bam_dirs_abs}; do
            for f in $d/*_{params.mhc_seq}_novoalign.*.sorted.filtered.bam*; do
                [ -e "$f" ] && ln -sf $(realpath $f) $stage/$(basename $f);
            done;
        done;
        cd $stage &&
        vep_tables=$(ls *.vep.txt 2>/dev/null || true) &&
        if [ -z "$vep_tables" ]; then
            echo "No mutations detected for patient {wildcards.patient_id} -- skipping mutation table";
            touch {params.mutations_abs};
        else
            tumour_ids="{params.tumour_normal_ids[0]}" && normal_ids="{params.tumour_normal_ids[1]}" &&
            tumour_bams=$(for tid in $tumour_ids; do ls ${{tid}}_{params.mhc_seq}_novoalign.*.sorted.filtered.bam 2>/dev/null; done) &&
            gl_bams=$(for nid in $normal_ids; do ls ${{nid}}_{params.mhc_seq}_novoalign.*.sorted.filtered.bam 2>/dev/null; done | sort -u) &&
            Rscript {params.scripts_dir}/bin/make_mutation_table.R --vep_tables $vep_tables --wxs_tumour_bam_files $tumour_bams --wxs_gl_bam_files $gl_bams --mutation_save_path {params.mutations_abs} --scripts_dir {params.scripts_dir}/bin/ --inventory_path {params.inventory_abs};
        fi
        ) > {log.stdout} 2>&1 &&
        rm -rf $stage
        """)


# Reformats a patient's {patient_id}_mutations.csv into a MAF-like table for tools that expect
# MAF input (e.g. maftools) but don't need real genomic coordinates for HLA -- mhc_hammer calls
# mutations against each patient's own personalised, allele-specific reference (one contig per
# typed HLA allele), so there's no shared genomic coordinate system across alleles/samples to
# report. See the header comment in src/mutations_to_maf.R for exactly which columns are real
# (Hugo_Symbol, Variant_Classification, allele calls, VEP annotations, read counts, sample IDs)
# vs. synthetic placeholders (Chromosome/Start_Position/End_Position/NCBI_Build). Module-owned
# script, not an upstream MHC Hammer bin/ file -- postprocesses this module's own output only, so
# (unlike most rules in this file) it needs no `container: None` for licensing/host-path reasons;
# left as None anyway since no container image is built for the mhc_hammer_r env in this module.
rule _mhc_hammer_mutations_to_maf:
    input:
        mutations = str(rules._mhc_hammer_parse_mutations.output.mutations)
    output:
        maf = CFG["dirs"]["mutations"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}_mutations.maf"
    log:
        stdout = CFG["logs"]["mutations"] + "{seq_type}--{genome_build}/{patient_id}/mutations_to_maf.log"
    params:
        script = MUTATIONS_TO_MAF_SCRIPT
    conda:
        CFG["conda_envs"]["mhc_hammer_r"]
    container:
        None
    threads:
        CFG["threads"]["mutations_to_maf"]
    resources:
        **CFG["resources"]["mutations_to_maf"]
    shell:
        op.as_one_line("""
        Rscript {params.script}
        --mutations_csv {input.mutations}
        --patient_id {wildcards.patient_id}
        --output_maf {output.maf}
        > {log.stdout} 2>&1
        """)


# Builds the cohort-wide gene-level summary table (docs/mhc_hammer_outputs.md upstream,
# DNA-only column subset), combining every patient's HLA-HD genotype, allele-specific coverage
# (mosdepth), DNA copy-number/allelic-imbalance, and mutation calls into one CSV. Cohort-wide,
# no wildcards. Mirrors upstream's CREATE_MHC_HAMMER_TABLE process.
#
# Like _mhc_hammer_parse_mutations, make_cohort_overview_table.R's `--csv_tables_path` is a
# manifest of *bare basenames* it reads relative to its own working directory (matches upstream's
# own `.map { it.getName() }.collectFile(...)` -- confirmed by reading the script and
# workflows/mhc_hammer.nf), not a list of absolute paths -- so inputs are staged flat with
# symlinks first, exactly as in _mhc_hammer_parse_mutations.
#
# Two fixes for real bugs in make_cohort_overview_table.R itself (upstream's own script,
# unmodifiable -- see licensing note near the top of this file), worked around from this rule's
# side instead:
# (1) `ls -1 . > input_csvs.txt` lists its own (already-created-by-the-redirect) output file,
#     feeding a bogus self-referencing row into the manifest -- piping through `grep -v` before
#     the final redirect avoids this (confirmed empirically: the file doesn't exist yet when `ls`
#     runs inside the pipeline, only once the whole pipeline's output is redirected at the end).
# (2) `transcriptome_hlahd_tables`'s loop (`for(line_idx in 1:nrow(transcriptome_hlahd_tables))`)
#     is the one table-type loop in the whole script *not* guarded by `if(nrow(...) > 0)` first --
#     every other category is. DNA-only v1 never produces a real *_transcriptome_allele_table.csv,
#     so nrow is always 0 there, and R's `1:0` evaluates to `c(1, 0)` (not an empty range),
#     crashing with `fread(NA)` ("missing value where TRUE/FALSE needed") -- confirmed on a real
#     run. Staging a placeholder CSV per patient gives the loop nrow==1 so it never hits the `1:0`
#     case. A *header-only* placeholder isn't enough, though (also confirmed on a real run, via a
#     follow-on "Incompatible join types: x.patient (logical) and i.patient (character)" error):
#     fread() on zero data rows has nothing to type-infer from and defaults every column to
#     logical. One dummy data row gives fread() real values to infer correct types from. The
#     `gene`/`allele1`/`allele2` fields are a fixed "NONE" sentinel (safe: no real HLA gene is ever
#     named "NONE", so this can never collide with genuine genome-derived rows for that patient
#     downstream), but the `patient` field must be that directory's own real patient_id, NOT a
#     shared fixed sentinel -- confirmed on a real multi-patient run: `hlahd_alleles <- merge(
#     genome_cohort_alleles, transcriptome_cohort_alleles, by = c("patient", "gene", "allele1",
#     "allele2", "homozygous"), all = TRUE)` is a *full outer join* (not the left join used
#     elsewhere in this same script for the wxs/rna library-size sections), so every patient's
#     placeholder row survives into the merged table rather than being dropped. A fixed sentinel
#     patient value (e.g. "__no_transcriptome_data__") would then make every patient's placeholder
#     row collide into the exact same (patient, gene) group once more than one patient exists,
#     tripping the script's own post-merge duplicate check ("Different HLA alleles from the genome
#     and transcriptome allele table."). Using the real patient_id keeps each placeholder row
#     unique across the cohort while still never matching real data (gene/allele1/allele2 stay
#     "NONE").
rule _mhc_hammer_cohort_table:
    input:
        dna_analysis = _mhc_hammer_completed_dna_analysis,
        mosdepth = _mhc_hammer_completed_mosdepth,
        library_size = _mhc_hammer_completed_library_size,
        hla_bam_read_count = _mhc_hammer_completed_hla_bam_read_count,
        patient_dirs = _mhc_hammer_completed_patient_dirs,
        inventory = str(rules._mhc_hammer_generate_inventory.output.inventory),
        hlahd_germline_samples = str(rules._mhc_hammer_generate_inventory.output.hlahd_germline_samples)
    output:
        cohort_table = CFG["dirs"]["cohort_tables"] + "cohort_mhc_hammer_gene_table.csv"
    log:
        stdout = CFG["logs"]["cohort_tables"] + "cohort_table.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        min_number_of_snps = CFG["options"]["min_number_of_snps"],
        max_copy_number_range = CFG["options"]["max_copy_number_range"],
        min_expected_depth = CFG["options"]["min_expected_depth"],
        min_frac_mapping_uniquely = CFG["options"]["min_frac_mapping_uniquely"],
        max_frac_mapping_multi_gene = CFG["options"]["max_frac_mapping_multi_gene"],
        # min_depth is per-seq_type (see default.yaml), but this rule is cohort-wide (not
        # wildcarded on seq_type) and can aggregate a mix of capture/genome pairs in one run.
        # make_cohort_overview_table.R's own use of --dna_snp_min_depth is purely a provenance
        # column (`min_dna_snp_depth_param`, recorded once per output row, never used to filter
        # anything -- the real per-seq_type filtering already happened upstream, correctly, in
        # _mhc_hammer_detect_cn_aib's get_filtered_pos_bed.R call, which does have a real
        # wildcards.seq_type). Recording the capture value here since that's this tool's original,
        # validated default; for a genome-seq_type row in a mixed cohort this column will read 30
        # even though 5 (or whatever options.min_depth.genome is set to) was the value actually
        # used to filter that row's SNPs -- cosmetic inaccuracy only, not a correctness issue.
        min_depth = CFG["options"]["min_depth"]["capture"],
        flat_files_abs = lambda wildcards, input: " ".join(os.path.abspath(f) for f in (
            list(input.dna_analysis) + list(input.mosdepth) + list(input.library_size) + list(input.hla_bam_read_count)
        )),
        patient_dirs_abs = lambda wildcards, input: " ".join(sorted(set(os.path.abspath(d) for d in input.patient_dirs))),
        inventory_abs = lambda wildcards, input: os.path.abspath(input.inventory),
        hlahd_germline_samples_abs = lambda wildcards, input: os.path.abspath(input.hlahd_germline_samples),
        cohort_table_abs = lambda wildcards, output: os.path.abspath(output.cohort_table)
    conda:
        CFG["conda_envs"]["mhc_hammer_r"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["cohort_table"]
    resources:
        **CFG["resources"]["cohort_table"]
    shell:
        op.as_one_line("""
        wd=$(dirname {output.cohort_table}) && mkdir -p $wd &&
        stage=$wd/staged && rm -rf $stage && mkdir -p $stage &&
        (
        for f in {params.flat_files_abs}; do ln -sf $f $stage/$(basename $f); done;
        for d in {params.patient_dirs_abs}; do
            for f in $d/*_genome_allele_table.csv; do
                [ -e "$f" ] && ln -sf $f $stage/$(basename $f);
            done;
            patient=$(basename $d) ;
            echo 'gene,allele1,allele2,patient,num_snps,homozygous' > $stage/${{patient}}_transcriptome_allele_table.csv ;
            echo "NONE,NONE,NONE,${{patient}},0,TRUE" >> $stage/${{patient}}_transcriptome_allele_table.csv ;
        done;
        cd $stage &&
        ( ls -1 . | grep -v '^input_csvs\.txt$' > input_csvs.txt || true ) &&
        if [ ! -s input_csvs.txt ]; then
            echo "No completed tumour/normal pair(s) available yet (input_csvs.txt is empty) -- writing an empty cohort table. This is expected on a first, from-scratch invocation; rerun the same command once --keep-going has finished attempting everything to get real results." ;
            touch {params.cohort_table_abs} ;
        else
            Rscript {params.scripts_dir}/bin/make_cohort_overview_table.R
              --inventory_path {params.inventory_abs}
              --csv_tables_path input_csvs.txt
              --hlahd_germline_samples_path {params.hlahd_germline_samples_abs}
              --max_cn_range {params.max_copy_number_range}
              --min_n_snps {params.min_number_of_snps}
              --min_expected_depth {params.min_expected_depth}
              --min_frac_mapping_uniquely {params.min_frac_mapping_uniquely}
              --max_frac_mapping_multi_gene {params.max_frac_mapping_multi_gene}
              --dna_snp_min_depth {params.min_depth} &&
            mv cohort_mhc_hammer_gene_table.csv {params.cohort_table_abs} ;
        fi
        ) > {log.stdout} 2>&1 &&
        rm -rf $stage
        """)


# Symlinks the final per-run output files into the module results directory (under '99-outputs/')
rule _mhc_hammer_output_dna_analysis:
    input:
        dna_analysis = str(rules._mhc_hammer_detect_cn_aib.output.dna_analysis)
    output:
        dna_analysis = CFG["dirs"]["outputs"] + "dna_analysis/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dna_analysis.csv"
    run:
        op.relative_symlink(input.dna_analysis, output.dna_analysis, in_module = True)

rule _mhc_hammer_output_mutations:
    input:
        mutations = str(rules._mhc_hammer_parse_mutations.output.mutations)
    output:
        mutations = CFG["dirs"]["outputs"] + "mutations/{seq_type}--{genome_build}/{patient_id}.mutations.csv"
    run:
        op.relative_symlink(input.mutations, output.mutations, in_module = True)

rule _mhc_hammer_output_mutations_maf:
    input:
        maf = str(rules._mhc_hammer_mutations_to_maf.output.maf)
    output:
        maf = CFG["dirs"]["outputs"] + "mutations_maf/{seq_type}--{genome_build}/{patient_id}.mutations.maf"
    run:
        op.relative_symlink(input.maf, output.maf, in_module = True)

rule _mhc_hammer_output_cohort_table:
    input:
        cohort_table = str(rules._mhc_hammer_cohort_table.output.cohort_table)
    output:
        cohort_table = CFG["dirs"]["outputs"] + "cohort_tables/cohort_mhc_hammer_gene_table.csv"
    run:
        op.relative_symlink(input.cohort_table, output.cohort_table, in_module = True)

rule _mhc_hammer_output_hla_final_result:
    input:
        hla_final_result = str(rules._mhc_hammer_hlahd.output.hla_final_result)
    output:
        hla_final_result = CFG["dirs"]["outputs"] + "hla_alleles/{seq_type}--{genome_build}/{patient_id}.hla_alleles.txt"
    run:
        op.relative_symlink(input.hla_final_result, output.hla_final_result, in_module = True)

rule _mhc_hammer_output_hla2_alleles:
    input:
        hla_alleles = str(rules._mhc_hammer_hla2_hlahd.output.hla_alleles)
    output:
        hla_alleles = CFG["dirs"]["outputs"] + "hla2_alleles/{seq_type}--{genome_build}/{patient_id}.hla2_alleles.csv"
    run:
        op.relative_symlink(input.hla_alleles, output.hla_alleles, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks. Uses
# CFG["paired_runs"] (not CFG["runs"]) so that tumour samples without a matched germline WES
# sample -- which this module cannot process at all, since HLA typing and the personalised
# reference both depend on the patient's own germline sample -- are never requested as targets.
rule _mhc_hammer_all:
    input:
        expand(
            [
                str(rules._mhc_hammer_output_dna_analysis.output.dna_analysis)
            ],
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            tumour_id = CFG["paired_runs"]["tumour_sample_id"],
            normal_id = CFG["paired_runs"]["normal_sample_id"],
            pair_status = CFG["paired_runs"]["pair_status"]
        ),
        expand(
            [
                str(rules._mhc_hammer_output_mutations.output.mutations)
            ],
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            patient_id = CFG["paired_runs"]["tumour_patient_id"]
        ),
        expand(
            [
                str(rules._mhc_hammer_output_mutations_maf.output.maf)
            ],
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            patient_id = CFG["paired_runs"]["tumour_patient_id"]
        ),
        str(rules._mhc_hammer_output_cohort_table.output.cohort_table),
        # HLA-HD's own class I result file, in standard nomenclature -- same patient-level
        # CFG["paired_runs"] scoping as the mutations/cohort targets above.
        expand(
            [
                str(rules._mhc_hammer_output_hla_final_result.output.hla_final_result)
            ],
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            patient_id = CFG["paired_runs"]["tumour_patient_id"]
        ),
        # HLA class II germline typing -- same CFG["paired_runs"] patient scoping as the
        # mutations target above (only patients with a real matched germline get requested).
        expand(
            [
                str(rules._mhc_hammer_output_hla2_alleles.output.hla_alleles)
            ],
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            patient_id = CFG["paired_runs"]["tumour_patient_id"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
