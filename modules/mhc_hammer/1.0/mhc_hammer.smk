#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Tom Jones & Clare Puttick (McGranahanLab/mhc-hammer)
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import os
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
VEP_CACHE = CFG["options"]["vep_cache"]
MHC_COORDS = CFG["options"]["mhc_coords"]

# Upstream's own scripts hardcode the literal tokens "wxs"/"rnaseq" (its internal DNA-vs-RNA
# vocabulary) into filenames and into the regexes several cohort-level R scripts use to find and
# parse those filenames (e.g. `make_cohort_overview_table.R` matches `_wxs.library_size.txt$` and
# `_wxs_novoalign.hla_bam_read_count.csv$` literally) -- this is independent of and NOT the same
# as lcr-modules' own `seq_type` wildcard (`capture`/`genome`). Confirmed by reading the actual
# scripts in a local clone. Every filename segment that flows into one of those scripts' pattern
# matching must use this literal constant, never `{wildcards.seq_type}` directly.
MHC_SEQ = "wxs"
IMGT_URLS = CFG["options"]["imgt_release_urls"][CFG["options"]["imgt_release"]]

# Define rules to be run locally when using a compute cluster
localrules:
    _mhc_hammer_input_bam,
    _mhc_hammer_download_reference,
    _mhc_hammer_generate_inventory,
    _mhc_hammer_output_dna_analysis,
    _mhc_hammer_output_mutations,
    _mhc_hammer_output_cohort_table,
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

def _mhc_hammer_get_germline_fqs(wildcards):
    # Every patient must have exactly one germline WES sample (v1 simplification -- see README).
    CFG = config["lcr-modules"]["mhc_hammer"]
    normal = op.filter_samples(
        CFG["samples"],
        patient_id = wildcards.patient_id,
        seq_type = wildcards.seq_type,
        genome_build = wildcards.genome_build,
    )
    normal = normal[normal["tissue_status"].str.lower() == "normal"]
    assert len(normal) == 1, (
        f"Expected exactly one germline WES sample for patient '{wildcards.patient_id}' "
        f"({wildcards.seq_type}--{wildcards.genome_build}), found {len(normal)}. This module "
        f"requires exactly one germline WES sample per patient (see README)."
    )
    sample_id = normal["sample_id"].tolist()[0]
    return {
        "fq1": expand(str(rules._mhc_hammer_generate_fqs.output.fq1), sample_id = sample_id, allow_missing = True),
        "fq2": expand(str(rules._mhc_hammer_generate_fqs.output.fq2), sample_id = sample_id, allow_missing = True),
    }

def _mhc_hammer_reference_dir_for_sample(wildcards):
    patient_id = _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id)
    return expand(str(rules._mhc_hammer_generate_references.output.patient_dir), patient_id = patient_id, allow_missing = True)

def _mhc_hammer_reference_dir_for_tumour(wildcards):
    patient_id = _mhc_hammer_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type)
    return expand(str(rules._mhc_hammer_generate_references.output.patient_dir), patient_id = patient_id, allow_missing = True)

def _mhc_hammer_get_purity_ploidy(wildcards):
    CFG = config["lcr-modules"]["mhc_hammer"]
    purity_ploidy_file = CFG["options"]["purity_ploidy_file"]
    if not purity_ploidy_file:
        return "", ""
    with open(purity_ploidy_file) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] == wildcards.tumour_id:
                return fields[2], fields[1] # purity, ploidy
    return "", ""


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mhc_hammer_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


# Downloads and unpacks the pinned IMGT-derived MHC Hammer reference bundle from Zenodo
# (docs/mhc_reference_files.md upstream). Two separate zip files per release. Cohort-wide,
# no wildcards -- built once and reused by every patient.
rule _mhc_hammer_download_reference:
    output:
        gtf = CFG["dirs"]["mhc_reference"] + "mhc_references/gtf/mhc.gtf",
        genome_fasta = CFG["dirs"]["mhc_reference"] + "mhc_references/genome/mhc_genome_strand.fasta",
        kmer_file = CFG["dirs"]["mhc_reference"] + "kmer_files/imgt_30mers.fa"
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
        mkdir -p {params.outdir} &&
        wget -qO {params.outdir}/mhc_references.zip "{params.mhc_references_url}" &&
        wget -qO {params.outdir}/kmer_files.zip "{params.kmer_files_url}" &&
        unzip -oq {params.outdir}/mhc_references.zip -d {params.outdir} &&
        unzip -oq {params.outdir}/kmer_files.zip -d {params.outdir} &&
        rm {params.outdir}/mhc_references.zip {params.outdir}/kmer_files.zip
        """)


# Builds a genome-wide 30-mer occurrence index from a user-supplied reference genome, used by
# _mhc_hammer_filter_kmers below to identify which of the IMGT-derived kmers aren't actually
# HLA-specific (see design note there). Cohort-wide, no wildcards, genome-build-independent by
# design (see options.reference_genome_fasta above) -- built once and reused for every patient.
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
        genome_fasta = CFG["options"]["reference_genome_fasta"]
    output:
        index = CFG["dirs"]["mhc_reference"] + "genome_kmer_index/genome_30mers.jf"
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
        mkdir -p $(dirname {output.index}) &&
        jellyfish count -m 30 -s 3G -t {threads} -o {output.index} {input.genome_fasta}
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
# NOTE: `jellyfish query -s <file> <index>` is assumed to batch-query every line of
# imgt_30mers.fa (one bare 30-mer per line, no FASTA headers) and print "<kmer> <count>" per
# line -- this wasn't independently verified against a real jellyfish install; check the actual
# output format on first real run.
rule _mhc_hammer_filter_kmers:
    input:
        index = str(rules._mhc_hammer_build_genome_kmer_index.output.index),
        raw_kmers = str(rules._mhc_hammer_download_reference.output.kmer_file)
    output:
        filtered_kmers = CFG["dirs"]["mhc_reference"] + "kmer_files/imgt_30mers.filtered.fa",
        counts = CFG["dirs"]["mhc_reference"] + "kmer_files/imgt_30mers.genome_occurrences.tsv"
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
        op.as_one_line("""
        jellyfish query -s {input.raw_kmers} {input.index} > {output.counts} &&
        awk -v max={params.max_occurrences} '$2 <= max {{print $1}}' {output.counts} > {output.filtered_kmers}
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
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat} &&
        samtools view -@ {threads} -c -f 1 -F 2308 {input.bam} > {output.library_size}
        """)


# Subsets the input BAM to reads that could plausibly align to the HLA class I genes: the MHC
# region on chr6, unmapped reads, and (fish_reads=True) reads matching the (genome-occurrence-
# filtered -- see _mhc_hammer_filter_kmers) IMGT 30-mer kmer set. Mirrors upstream's SUBSET_BAMS
# process, with the kmer file filtered rather than used raw.
rule _mhc_hammer_subset_bam:
    input:
        bam = str(rules._mhc_hammer_input_bam.output.bam),
        kmer_file = str(rules._mhc_hammer_filter_kmers.output.filtered_kmers)
    output:
        bam = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.subset.sorted.bam",
        bai = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.subset.sorted.bam.bai",
        read_counts = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.read_counts.csv"
    params:
        scripts_dir = SCRIPTS_DIR,
        mhc_chrom = lambda w: MHC_COORDS[w.genome_build],
        # subset_bam_opt.sh does literal `[ "$var" == true ]` string comparisons (lowercase) --
        # str(True)/str(False) from a Python/YAML bool renders as "True"/"False" (capitalised)
        # when substituted into the shell command, which never matches. Must be lowercased here.
        fish_reads = str(CFG["options"]["fish_reads"]).lower(),
        unmapped_reads = str(CFG["options"]["unmapped_reads"]).lower(),
        contig_reads = str(CFG["options"]["contig_reads"]).lower(),
        sort_mem = lambda wildcards, resources: max(1, int(resources.mem_mb / 1000 * 0.8)),
        bam_abs = lambda wildcards, input: os.path.abspath(input.bam),
        kmer_file_abs = lambda wildcards, input: os.path.abspath(input.kmer_file)
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
        echo {params.mhc_chrom} > $wd/mhc_coords.txt &&
        touch $wd/contigs_placeholder.txt &&
        rm -rf $wd/{wildcards.sample_id}_tmpDir &&
        (
        cd $wd &&
        {params.scripts_dir}/bin/subset_bam_opt.sh
        -b {params.bam_abs}
        -k {params.kmer_file_abs}
        -f {params.fish_reads} -c {params.contig_reads} -d contigs_placeholder.txt
        -u {params.unmapped_reads} -h mhc_coords.txt
        -p {wildcards.sample_id} -t {threads}
        -m {params.sort_mem}G -o false
        )
        """)


# Converts the subsetted BAM to paired FASTQs (samtools collate/fastq), mirroring upstream's
# GENERATE_HLA_FQS process. These are the reads aligned against the personalised reference.
rule _mhc_hammer_generate_fqs:
    input:
        bam = str(rules._mhc_hammer_subset_bam.output.bam)
    output:
        fq1 = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.1.fq.gz",
        fq2 = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.2.fq.gz"
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
        samtools collate -@ {threads} -u -O {input.bam} |
        samtools fastq -@ {threads} -1 {output.fq1} -2 {output.fq2} -s /dev/null -0 /dev/null -n
        """)


# Runs the user-supplied HLA-HD install on the patient's one germline WES sample to type
# HLA-A/B/C. Patient-level (once per patient, not once per sample). Mirrors upstream's
# HLAHD_LOCAL process. Fails loudly (rather than gracefully excluding the patient, as upstream
# does) if A/B/C are not all typed -- see README v1 simplifications.
rule _mhc_hammer_hlahd:
    input:
        unpack(_mhc_hammer_get_germline_fqs),
        gtf = str(rules._mhc_hammer_download_reference.output.gtf)
    output:
        hla_alleles = CFG["dirs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}_hla_alleles.csv",
        result_dir = directory(CFG["dirs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/result")
    log:
        stdout = CFG["logs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/hlahd.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        hlahd_dir = HLAHD_DIR,
        workdir = lambda wildcards, output: os.path.dirname(output.hla_alleles),
        gtf_abs = lambda wildcards, input: os.path.abspath(input.gtf)
    conda:
        CFG["conda_envs"]["mhc_hammer_hlahd"]
    container:
        None
    threads:
        CFG["threads"]["hlahd"]
    resources:
        **CFG["resources"]["hlahd"]
    shell:
        op.as_one_line("""
        rm -rf {params.workdir}/{wildcards.patient_id} &&
        mkdir -p {params.workdir} &&
        (
        export PATH=${{PATH}}:{params.hlahd_dir}/bin &&
        bash {params.hlahd_dir}/bin/hlahd.sh -m 100 -c 1.0 -t {threads}
        -f {params.hlahd_dir}/freq_data
        {input.fq1} {input.fq2}
        {params.hlahd_dir}/HLA_gene.split.txt
        {params.hlahd_dir}/dictionary
        {wildcards.patient_id} {params.workdir} &&
        rm -rf {params.workdir}/{wildcards.patient_id}/exon {params.workdir}/{wildcards.patient_id}/intron
                {params.workdir}/{wildcards.patient_id}/mapfile {params.workdir}/{wildcards.patient_id}/maplist &&
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
        else
            echo "ERROR: HLA-HD failed to produce HLA A, B and C estimates for patient {wildcards.patient_id}. See {log.stdout}." >&2 &&
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
        op.as_one_line("""
        export PATH=${{PATH}}:{params.novoalign_dir} &&
        wd=$(dirname {output.bam}) &&
        mkfifo $wd/fq1_uncompressed $wd/fq2_uncompressed &&
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
        wd=$(dirname {output.bam}) &&
        samtools sort -@ {threads} -o $wd/{wildcards.sample_id}.sorted.bam {input.bam} &&
        samtools view -@ {threads} -f 2 -b -o $wd/{wildcards.sample_id}.hla.bam $wd/{wildcards.sample_id}.sorted.bam &&
        samtools addreplacerg -r ID:{wildcards.sample_id} -r SM:{wildcards.sample_id} -o {output.bam} $wd/{wildcards.sample_id}.hla.bam &&
        samtools index -@ {threads} {output.bam} &&
        rm $wd/{wildcards.sample_id}.sorted.bam $wd/{wildcards.sample_id}.hla.bam
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
    params:
        scripts_dir = SCRIPTS_DIR,
        max_mismatch = CFG["options"]["max_mismatch"],
        mhc_seq = MHC_SEQ,
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id),
        bam_abs = lambda wildcards, input: os.path.abspath(input.bam)
    conda:
        CFG["conda_envs"]["mhc_hammer_alignment"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["make_allele_bams"]
    resources:
        **CFG["resources"]["make_allele_bams"]
    shell:
        op.as_one_line("""
        wd=$(dirname {output.passed_hla_genes}) &&
        mkdir -p $wd &&
        (cd $wd &&
         {params.scripts_dir}/bin/make_hla_bams.sh {params.bam_abs} {params.scripts_dir}/
           {params.max_mismatch} {wildcards.sample_id} {params.mhc_seq} novoalign)
        """)


# Computes per-allele coverage with mosdepth against the patient's personalised reference BED.
# Sample-level. Mirrors upstream's MOSDEPTH process.
rule _mhc_hammer_mosdepth:
    input:
        allele_bams_marker = str(rules._mhc_hammer_make_allele_bams.output.passed_hla_alleles),
        patient_dir = _mhc_hammer_reference_dir_for_sample
    output:
        bed = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}." + MHC_SEQ + ".novoalign.mosdepth.bed"
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
        wd=$(dirname {input.allele_bams_marker}) &&
        > {output.bed} &&
        for bam_file in $wd/*_{params.mhc_seq}_novoalign.*.sorted.filtered.bam; do
            [ -e "$bam_file" ] || continue;
            allele=$(basename $bam_file | sed -E 's/^.*_novoalign\\.//; s/\\.sorted\\.filtered\\.bam$//');
            mosdepth --no-per-base --chrom $allele --by {params.bed_file} $wd/{wildcards.sample_id}.$allele $bam_file &&
            gunzip -f $wd/{wildcards.sample_id}.$allele.regions.bed.gz &&
            cat $wd/{wildcards.sample_id}.$allele.regions.bed >> {output.bed};
        done
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
        "tumour_flagstat": expand(str(rules._mhc_hammer_flagstat.output.flagstat), sample_id = wildcards.tumour_id, allow_missing = True),
        "normal_flagstat": expand(str(rules._mhc_hammer_flagstat.output.flagstat), sample_id = wildcards.normal_id, allow_missing = True),
        "patient_dir": patient_dir
    }


# Detects HLA allele-specific copy number and allelic imbalance (LOH) by comparing tumour and
# germline coverage at germline-heterozygous SNP positions within each HLA gene the tumour and
# normal both passed filtering for. Pair-level. Mirrors upstream's DETECT_CN_AND_AIB process
# ("all_snps" variant only -- upstream's "exon_snps" variant is disabled/commented out upstream
# too). Copy-number output is skipped (but allelic-imbalance/LOH detection still runs) for
# samples missing a `purity_ploidy_file` entry, matching upstream's `missing_purity_ploidy` flag.
rule _mhc_hammer_detect_cn_aib:
    input:
        unpack(lambda wildcards: _mhc_hammer_pair_inputs(wildcards, str(rules._mhc_hammer_make_allele_bams.output.passed_heterozygous_hla_genes)))
    output:
        # Filename must match concatenate_dna_analysis_tables.R's own
        # paste0(sample_name,"_",snp_type,"_",aligner,"_dna_analysis.csv") convention exactly --
        # confirmed by reading that script.
        dna_analysis = CFG["dirs"]["dna_analysis"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_all_snps_novoalign_dna_analysis.csv"
    log:
        stdout = CFG["logs"]["dna_analysis"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/detect_cn_aib.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        min_depth = CFG["options"]["min_depth"],
        mhc_seq = MHC_SEQ,
        purity_ploidy = lambda wildcards: _mhc_hammer_get_purity_ploidy(wildcards),
        tumour_marker_abs = lambda wildcards, input: os.path.abspath(input.tumour_marker[0]),
        normal_marker_abs = lambda wildcards, input: os.path.abspath(input.normal_marker[0]),
        tumour_flagstat_abs = lambda wildcards, input: os.path.abspath(input.tumour_flagstat[0]),
        normal_flagstat_abs = lambda wildcards, input: os.path.abspath(input.normal_flagstat[0]),
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
                Rscript {params.scripts_dir}/bin/get_cn.R --allele1 $allele1 --allele2 $allele2 --allele1_gl_coverage_file {wildcards.normal_id}.$allele1.coverage_at_filtered_positions.csv --allele2_gl_coverage_file {wildcards.normal_id}.$allele2.coverage_at_filtered_positions.csv --allele1_tumour_coverage_file {wildcards.tumour_id}.$allele1.coverage_at_filtered_positions.csv --allele2_tumour_coverage_file {wildcards.tumour_id}.$allele2.coverage_at_filtered_positions.csv --allele1_snp_bed $allele1_gene_filtered_snps --allele2_snp_bed $allele2_gene_filtered_snps --purity $purity --ploidy $ploidy --tumour_library_size_path {params.tumour_flagstat_abs} --gl_library_size_path {params.normal_flagstat_abs} --gtf_path {params.patient_dir_abs}/{params.patient_id}.gtf --cn_output_path {wildcards.tumour_id}.$gene.all_snps.cn.csv --cn_plots_prefix {wildcards.tumour_id}.$gene.all_snps.cn --scripts_dir {params.scripts_dir}/bin/;
                Rscript {params.scripts_dir}/bin/get_expected_depth.R --allele1 $allele1 --allele2 $allele2 --allele1_gl_coverage_file {wildcards.normal_id}.$allele1.coverage_at_filtered_positions.csv --allele2_gl_coverage_file {wildcards.normal_id}.$allele2.coverage_at_filtered_positions.csv --allele1_snp_bed $allele1_gene_filtered_snps --allele2_snp_bed $allele2_gene_filtered_snps --purity $purity --tumour_library_size_path {params.tumour_flagstat_abs} --gl_library_size_path {params.normal_flagstat_abs} --expected_depth_output_path {wildcards.tumour_id}.$gene.all_snps.expected_depth.csv --expected_depth_plots_prefix {wildcards.tumour_id}.$gene.all_snps.expected_depth --scripts_dir {params.scripts_dir}/bin/;
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
            Rscript {params.scripts_dir}/bin/get_logr_aib.R --allele1 $allele1 --allele2 $allele2 --allele1_gl_reads_count_once_coverage {wildcards.normal_id}.$allele1.coverage_at_filtered_snps_reads_count_once.csv --allele2_gl_reads_count_once_coverage {wildcards.normal_id}.$allele2.coverage_at_filtered_snps_reads_count_once.csv --allele1_tumour_reads_count_once_coverage {wildcards.tumour_id}.$allele1.coverage_at_filtered_snps_reads_count_once.csv --allele2_tumour_reads_count_once_coverage {wildcards.tumour_id}.$allele2.coverage_at_filtered_snps_reads_count_once.csv --allele1_snp_bed $allele1_gene_filtered_snps --allele2_snp_bed $allele2_gene_filtered_snps --tumour_library_size_path {params.tumour_flagstat_abs} --gl_library_size_path {params.normal_flagstat_abs} --logr_aib_output_path {wildcards.tumour_id}.$gene.all_snps.logr_aib.csv --logr_aib_plots_prefix {wildcards.tumour_id}.$gene.all_snps.logr_aib --scripts_dir {params.scripts_dir}/bin/;
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
# this rule derives the same columns from CFG["samples"]/CFG["paired_runs"]/purity_ploidy_file
# instead. Confirmed against the actual column names those scripts read (`patient`, `sample_name`,
# `sample_type`, `sequencing_type`, `purity`, `ploidy`, `normal_sample_name`) by reading a local
# clone of mhc-hammer's bin/*.R. Cohort-wide, no wildcards, stdlib-only (no pandas dependency).
rule _mhc_hammer_generate_inventory:
    output:
        inventory = CFG["dirs"]["cohort_tables"] + "inventory.csv",
        hlahd_germline_samples = CFG["dirs"]["cohort_tables"] + "hlahd_germline_samples.txt"
    run:
        import csv
        CFG = config["lcr-modules"]["mhc_hammer"]
        samples = CFG["samples"]
        runs = CFG["paired_runs"]
        purity_ploidy_file = CFG["options"]["purity_ploidy_file"]
        purity_ploidy = {}
        if purity_ploidy_file:
            with open(purity_ploidy_file) as fh:
                for line in fh:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 3:
                        purity_ploidy[fields[0]] = (fields[2], fields[1])  # sample_id -> (purity, ploidy)
        normal_for_tumour = dict(zip(runs["tumour_sample_id"], runs["normal_sample_id"]))
        os.makedirs(os.path.dirname(output.inventory), exist_ok=True)
        with open(output.inventory, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["patient", "sample_name", "sample_type", "sequencing_type", "purity", "ploidy", "normal_sample_name"])
            for _, row in samples.iterrows():
                sample_id = row["sample_id"]
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
rule _mhc_hammer_cohort_table:
    input:
        dna_analysis = expand(
            str(rules._mhc_hammer_detect_cn_aib.output.dna_analysis),
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            tumour_id = CFG["paired_runs"]["tumour_sample_id"],
            normal_id = CFG["paired_runs"]["normal_sample_id"],
            pair_status = CFG["paired_runs"]["pair_status"]
        ),
        mosdepth = expand(
            str(rules._mhc_hammer_mosdepth.output.bed),
            zip,
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"],
            sample_id = CFG["samples"]["sample_id"]
        ),
        library_size = expand(
            str(rules._mhc_hammer_flagstat.output.library_size),
            zip,
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"],
            sample_id = CFG["samples"]["sample_id"]
        ),
        hla_bam_read_count = expand(
            str(rules._mhc_hammer_make_allele_bams.output.hla_bam_read_count),
            zip,
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"],
            sample_id = CFG["samples"]["sample_id"]
        ),
        patient_dirs = expand(
            str(rules._mhc_hammer_generate_references.output.patient_dir),
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            patient_id = CFG["paired_runs"]["tumour_patient_id"]
        ),
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
        min_depth = CFG["options"]["min_depth"],
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
        done;
        cd $stage &&
        ls -1 . > input_csvs.txt &&
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
        mv cohort_mhc_hammer_gene_table.csv {params.cohort_table_abs}
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

rule _mhc_hammer_output_cohort_table:
    input:
        cohort_table = str(rules._mhc_hammer_cohort_table.output.cohort_table)
    output:
        cohort_table = CFG["dirs"]["outputs"] + "cohort_tables/cohort_mhc_hammer_gene_table.csv"
    run:
        op.relative_symlink(input.cohort_table, output.cohort_table, in_module = True)


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
        str(rules._mhc_hammer_output_cohort_table.output.cohort_table)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
