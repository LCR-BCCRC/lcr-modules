#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Hartwig Medical Foundation (hartwigmedical/hmftools, LILAC)
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####

import os
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version = "1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")


# Setup module and store module-specific configuration in `CFG`
CFG = op.setup_module(
    name = "lilac",
    version = "1.0",
    subdirectories = ["inputs", "reference", "pave", "lilac", "outputs"],
)

# HLA typing genuinely needs the patient's own germline sample -- a substitute normal from a
# different patient (oncopipe's "unmatched" pair_status) is meaningless for this purpose, same
# reasoning as modules/mhc_hammer/1.0. CFG["paired_runs"] already exists from op.setup_module()
# (a standard oncopipe field, not something this module derives from scratch); this narrows it in
# place, same as mhc_hammer's own CFG["paired_runs"] = CFG["paired_runs"][... == "matched"].
CFG["paired_runs"] = CFG["paired_runs"][CFG["paired_runs"]["pair_status"] == "matched"]

# genome_build -> LILAC's own -ref_genome_version argument ("37" or "38", confirmed against the
# real installed hmftools-lilac=2.0.1 jar's own registered config item -- NOT "V37"/"V38" as
# LILAC's own README example commands show, though both forms were empirically confirmed to parse
# without error). Mirrors modules/hmftools/1.1's own options.version_map exactly (same problem,
# same genome_build flavour vocabulary), kept as a separate copy rather than a cross-module
# reference -- modules/sage/1.1 does the same rather than sharing hmftools/1.1's copy.
VERSION_MAP = CFG["options"]["version_map"]

possible_genome_builds = ", ".join(list(VERSION_MAP.keys()))
for genome_build in CFG["paired_runs"]["tumour_genome_build"]:
    assert genome_build in VERSION_MAP, (
        f"Samples table includes genome build '{genome_build}', not yet covered by "
        f"options.version_map (currently: {possible_genome_builds}). Add it there before rerunning."
    )

# GRCh38 BAMs must be aligned to a build with no HLA ALT contigs (LILAC's own README warning) --
# mirrors modules/hmftools/1.1's own use_masked_ref switch and masked_string construction exactly.
masked_string = ""
if CFG["options"]["use_masked_ref"]:
    masked_string = "_masked"

# Static passthrough of LILAC's own tunable parameters -- built once at parse time (not per-rule)
# since none of these depend on wildcards. Every key here is checked against the real installed
# hmftools-lilac=2.0.1 jar's own registered config items (`lilac -h`'s error-dump output), not just
# LILAC's own README table -- the two disagree in at least one place: the README documents a
# `-min_base_qual` option that doesn't exist in the real v2.0.1 CLI at all (confirmed on a real
# run: "unregistered config item: -min_base_qual"), so it's deliberately not included here.
# "genes" is included deliberately (restricted to MHC_CLASS_1, not LILAC's own "ALL" default) --
# see the long comment on options.genes in default.yaml for why: a real run crashed trying to
# solve for HLA-DQA1 (a class II gene "ALL" includes in this build) against class-I-only
# reference data.
LILAC_TUNABLE_FLAGS = " ".join(
    f"-{key} {CFG['options'][key]}"
    for key in [
        "genes", "min_evidence_support", "min_evidence_factor",
        "min_high_qual_evidence_factor", "min_fragments_per_allele",
        "min_fragments_to_remove_single", "top_score_threshold",
        "hla_y_threshold", "freq_score_penalty", "write_types"
    ]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _lilac_input_bam,
    _lilac_download_reference,
    _lilac_download_ensembl_cache,
    _lilac_output_tsv,
    _lilac_output_qc,
    _lilac_output_somatic_vcf,
    _lilac_all,


##### HELPER FUNCTIONS #####


# Battenberg/PURPLE-style optional-input pattern (see modules/mhc_hammer/1.0's
# options.battenberg_cellularity_ploidy, now inputs.battenberg_cellularity_ploidy) -- a pair
# missing this file just runs LILAC without it (no CN or somatic-mutation-vs-CN output), rather
# than being required up front. Re-fetches CFG from config[...] inline (not the module-level CFG
# name) since op.cleanup_module(CFG) deletes that name before Snakemake evaluates input functions.
# Uses os.path.isfile (not os.path.exists) deliberately -- a real run showed an unset/misconfigured
# inputs.gene_copy_number resolving to a bare, non-templated directory path (the results root
# itself), which os.path.exists() happily accepted (directories exist too), silently handing
# "-gene_copy_number <a directory>" to LILAC. isfile() can't be fooled the same way.
def _lilac_get_gene_copy_number_input(wildcards):
    CFG = config["lcr-modules"]["lilac"]
    pattern = CFG["inputs"]["gene_copy_number"]
    if not pattern:
        return []
    path = pattern.format(
        seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id, normal_id = wildcards.normal_id, pair_status = wildcards.pair_status
    )
    return [path] if os.path.isfile(path) else []

# The raw, as-supplied somatic VCF (e.g. SAGE's combined output) -- existence-gated the same way as
# gene_copy_number above. This is NOT what _lilac_run consumes directly (see
# _lilac_get_somatic_vcf_input below) -- it's raw SAGE/caller output with no PAVE annotation, and
# LILAC's own SomaticVariantAnnotation requires the VCF INFO tag PAVE writes ("IMPACT") to identify
# each variant's gene/coding-effect. Real, silent, fatal-to-correctness bug confirmed on hundreds of
# real samples: without that tag, VariantContextDecorator.gene() returns "", HlaGene.fromString("")
# throws internally and is caught as null, and every variant ends up assigned to no real gene --
# LILAC completes normally throughout, but SomaticMissense/SomaticNonsenseOrFrameshift/SomaticSplice
# are 0.0 for every allele, every sample, regardless of real mutation burden. _lilac_pave_annotate
# below is the fix -- this function is now just its gate/input.
def _lilac_get_raw_somatic_vcf_input(wildcards):
    CFG = config["lcr-modules"]["lilac"]
    pattern = CFG["inputs"]["somatic_vcf"]
    if not pattern:
        return []
    path = pattern.format(
        seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id, normal_id = wildcards.normal_id, pair_status = wildcards.pair_status
    )
    return [path] if os.path.isfile(path) else []

# What _lilac_run actually consumes as somatic_vcf -- PAVE's IMPACT-annotated rewrite of the raw
# VCF above, not the raw VCF itself. Gated on the same raw-file check (so a pair with no somatic_vcf
# still runs LILAC without one, exactly as before); _lilac_pave_annotate's own output path is
# returned directly (not existence-checked -- Snakemake is responsible for producing it).
def _lilac_get_somatic_vcf_input(wildcards):
    if not _lilac_get_raw_somatic_vcf_input(wildcards):
        return []
    path = str(rules._lilac_pave_annotate.output.vcf).format(
        seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id, normal_id = wildcards.normal_id, pair_status = wildcards.pair_status
    )
    return [path]

# PAVE's Ensembl gene/transcript cache (-ensembl_data_dir) -- NOT something LILAC itself needs.
# modules/hmftools/1.1 already downloads the exact same reference data (same bcgsc.ca mirror, same
# zip) for AMBER/COBALT/PURPLE/LINX -- set inputs.ensembl_data_dir to that module's own cache
# directory (e.g. results/hmftools-1.1/00-inputs/references/{genome_build}/ensembl_cache/) to avoid
# this module downloading a second, byte-identical copy. Leave unset (default) and this module
# downloads its own via _lilac_download_ensembl_cache -- kept optional, not required, so lilac still
# works standalone without hmftools/1.1 also being included in the runner.
# ancient() on the reuse path: this directory isn't produced by any rule in this DAG, so nothing
# should ever be "newer" than it in a way that would trigger a rerun.
def _lilac_get_ensembl_data_dir_input(wildcards):
    CFG = config["lcr-modules"]["lilac"]
    pattern = CFG["inputs"].get("ensembl_data_dir", "")
    if pattern:
        path = pattern.format(seq_type = wildcards.seq_type, genome_build = wildcards.genome_build)
        if os.path.isdir(path):
            return ancient(path)
    return str(rules._lilac_download_ensembl_cache.output.cache).format(genome_build = wildcards.genome_build)

# Both tumour and normal BAMs for a pair, from the sample-level _lilac_input_bam rule --
# allow_missing=True keeps {seq_type}/{genome_build} as real wildcards (filled in automatically
# from the calling rule's own wildcards), only sample_id is set explicitly per role.
def _lilac_get_pair_bams(wildcards):
    return {
        "tumour_bam": expand(str(rules._lilac_input_bam.output.bam), sample_id = wildcards.tumour_id, allow_missing = True),
        "tumour_bai": expand(str(rules._lilac_input_bam.output.bai), sample_id = wildcards.tumour_id, allow_missing = True),
        "normal_bam": expand(str(rules._lilac_input_bam.output.bam), sample_id = wildcards.normal_id, allow_missing = True),
        "normal_bai": expand(str(rules._lilac_input_bam.output.bai), sample_id = wildcards.normal_id, allow_missing = True),
    }


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/'). Provides both
# .bam and .cram naming for the data file and its index -- some "sample_bam" inputs in this repo
# are actually CRAM content symlinked/named as .bam (htslib content-sniffs the real format
# regardless of extension), and CRAM readers then look for a .crai sidecar, not .bam.bai. Mirrors
# modules/mhc_hammer/1.0's own _mhc_hammer_input_bam (itself mirroring modules/pathseq/1.0).
# NOTE: whether LILAC's own htsjdk-based BAM/CRAM reading needs this same treatment hasn't been
# independently confirmed on real CRAM input -- kept as cheap insurance either way.
rule _lilac_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        cram = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.cram",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bam, output.cram)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Downloads LILAC's own 3 required resource CSVs directly from a pre-extracted mirror on
# www.bcgsc.ca/downloads/morinlab/hmftools-references/lilac/ -- the same convention
# modules/hmftools/1.1 and modules/sage/1.1 already use for this exact upstream tool suite
# (confirmed directory listing, 2026-08-26: all 3 files here are plain, build-independent
# filenames, ~30MB total combined -- unlike the ~5.6GB-per-build official HMF/oncoanalyser bundle
# this rule originally pointed at. The mirror also has hla.37.bed/hla.38.bed -- the one file in
# this resource set that *is* build-specific -- but LILAC's own -resource_dir requirement, per its
# real registered CLI config items, is only these 3 CSVs; the bed file isn't fetched here). No
# genome_build/ref_genome_version wildcard needed at all: one cohort-wide download, reused by
# every pair regardless of build.
rule _lilac_download_reference:
    output:
        freq = CFG["dirs"]["reference"] + "lilac_allele_frequencies.csv",
        nuc = CFG["dirs"]["reference"] + "hla_ref_nucleotide_sequences.csv",
        aa = CFG["dirs"]["reference"] + "hla_ref_aminoacid_sequences.csv"
    log:
        stdout = CFG["logs"]["reference"] + "download_reference.log"
    params:
        url = "https://www.bcgsc.ca/downloads/morinlab/hmftools-references/lilac"
    conda:
        CFG["conda_envs"]["wget"]
    container:
        None
    threads:
        CFG["threads"]["download_reference"]
    resources:
        **CFG["resources"]["download_reference"]
    shell:
        op.as_one_line("""
        (
        wget -qO {output.freq} {params.url}/lilac_allele_frequencies.csv &&
        wget -qO {output.nuc} {params.url}/hla_ref_nucleotide_sequences.csv &&
        wget -qO {output.aa} {params.url}/hla_ref_aminoacid_sequences.csv
        ) > {log.stdout} 2>&1
        """)


# Downloads PAVE's Ensembl gene/transcript cache, genome-build-keyed -- same
# www.bcgsc.ca/downloads/morinlab/hmftools-references/ensembl_data_cache/{build}.zip mirror and
# VERSION_MAP values modules/hmftools/1.1's own _hmftools_get_ensembl_cache already uses for this
# exact upstream reference data (software-version-independent, safe to mirror as this module's own
# copy rather than depending cross-module on hmftools/1.1 -- same reasoning as
# _lilac_download_reference above). Only PAVE needs this; LILAC itself doesn't.
# Only actually in the DAG when _lilac_get_ensembl_data_dir_input's reuse check (inputs.ensembl_data_dir)
# doesn't find an existing directory to point at instead -- see that function's comment. Set that
# config value (e.g. to modules/hmftools/1.1's own cache path) to skip this download entirely and
# reuse an existing copy of the same data.
rule _lilac_download_ensembl_cache:
    output:
        cache = directory(CFG["dirs"]["reference"] + "{genome_build}/ensembl_cache/"),
        complete = touch(CFG["dirs"]["reference"] + "{genome_build}/ensembl_cache/cache.complete")
    log:
        stdout = CFG["logs"]["reference"] + "{genome_build}/download_ensembl_cache.log"
    params:
        url = "https://www.bcgsc.ca/downloads/morinlab/hmftools-references/ensembl_data_cache",
        alt_build = lambda w: VERSION_MAP[w.genome_build]
    conda:
        CFG["conda_envs"]["wget"]
    container:
        None
    threads:
        CFG["threads"]["download_reference"]
    resources:
        **CFG["resources"]["download_reference"]
    shell:
        op.as_one_line("""
        (
        wget -qO {output.cache}/{params.alt_build}.zip {params.url}/{params.alt_build}.zip &&
        unzip -qod {output.cache} {output.cache}/{params.alt_build}.zip
        ) > {log.stdout} 2>&1
        """)


# Annotates a raw somatic VCF (e.g. SAGE's combined output) with the VCF INFO tag ("IMPACT") LILAC
# requires to identify each variant's gene/coding-effect -- see the long comment on
# _lilac_get_raw_somatic_vcf_input for the real, silent-failure bug this fixes (every sample
# reporting 0.0 somatic mutations regardless of real mutation burden, confirmed on a real cohort of
# hundreds). Real, verified-minimal PAVE invocation: confirmed by downloading the actual
# pave-v1.9.1 release jar directly from GitHub and running it locally with a plain `java` (bypasses
# conda entirely -- same trick used for LILAC's and PURPLE's own CLI verification earlier), then
# reading its registered CLI config items -- only -sample/-input_vcf/-ref_genome/
# -ref_genome_version/-ensembl_data_dir are REQUIRED; driver_gene_panel/gnomad/clinvar/blacklist/PON
# are all optional and safely omitted here, since the only thing this module needs from PAVE is the
# IMPACT tag, not a full driver-annotation pass.
# Only ever in the DAG for a pair that actually has a real somatic_vcf on disk -- nothing downstream
# requests this rule's output otherwise (same existence-gated-optional-input mechanism already used
# for gene_copy_number/somatic_vcf themselves).
rule _lilac_pave_annotate:
    input:
        vcf = _lilac_get_raw_somatic_vcf_input,
        ensembl_cache = _lilac_get_ensembl_data_dir_input,
        fasta = reference_files("genomes/{genome_build}" + masked_string + "/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["pave"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.pave.vcf.gz"
    log:
        stdout = CFG["logs"]["pave"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pave.log"
    params:
        ref_genome_version = lambda wildcards: VERSION_MAP[wildcards.genome_build],
        # The real, downloaded ensembl_data_cache zip unpacks into a build-versioned subdirectory
        # (confirmed by downloading the real 37.zip and inspecting it directly: its files are all
        # under "37/", e.g. "37/ensembl_gene_data.csv", not at the top level) -- true whether this
        # module downloaded the cache itself or it's an externally-supplied one following the same
        # modules/hmftools/1.1-established layout, so this appends the same subdirectory either way.
        ensembl_dir = lambda wildcards, input: os.path.join(str(input.ensembl_cache), VERSION_MAP[wildcards.genome_build]),
        output_dir = lambda wildcards, output: os.path.dirname(output.vcf),
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["pave"]
    container:
        None
    threads:
        CFG["threads"]["pave"]
    resources:
        **CFG["resources"]["pave"]
    shell:
        op.as_one_line("""
        mkdir -p {params.output_dir} &&
        pave -Xmx{params.jvmheap}m
        -sample {wildcards.tumour_id}
        -input_vcf {input.vcf[0]}
        -ref_genome {input.fasta}
        -ref_genome_version {params.ref_genome_version}
        -ensembl_data_dir {params.ensembl_dir}
        -output_dir {params.output_dir}
        -output_vcf {output.vcf}
        -threads {threads}
        > {log.stdout} 2>&1
        """)


# Core rule: HLA class I typing (always), plus -- when the optional gene_copy_number/somatic_vcf
# inputs are available for this pair -- tumour allele-specific copy number and somatic mutation
# assignment per allele, all in one LILAC invocation (no personalised-reference construction,
# realignment, or allele-specific BAM splitting needed -- unlike modules/mhc_hammer/1.0, LILAC
# does its own HLA-region fragment filtering directly against the full BAM). One rule per matched
# tumour/normal pair.
#
# LILAC writes its own output files directly into -output_dir using -sample as the filename
# prefix; output.somatic_vcf is touched empty when somatic_vcf wasn't provided, since LILAC itself
# won't create that file in that case and Snakemake requires every declared output to exist.
#
# NOTE (from LILAC's own README, not this module's design): if more than ~9% of the HLA-A/B/C
# coding regions have less than 10x coverage, LILAC fails outright rather than typing the sample
# -- a genuinely expected outcome for real low-coverage samples, same class of real, legitimate
# per-sample failure documented for modules/mhc_hammer/1.0, not something to work around here.
rule _lilac_run:
    input:
        unpack(_lilac_get_pair_bams),
        gene_copy_number = _lilac_get_gene_copy_number_input,
        somatic_vcf = _lilac_get_somatic_vcf_input,
        resource_freq = str(rules._lilac_download_reference.output.freq),
        resource_nuc = str(rules._lilac_download_reference.output.nuc),
        resource_aa = str(rules._lilac_download_reference.output.aa),
        fasta = reference_files("genomes/{genome_build}" + masked_string + "/genome_fasta/genome.fa")
    output:
        tsv = CFG["dirs"]["lilac"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.lilac.tsv",
        qc = CFG["dirs"]["lilac"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.lilac.qc.tsv",
        candidates = CFG["dirs"]["lilac"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.lilac.candidates.coverage.tsv",
        somatic_vcf = CFG["dirs"]["lilac"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.lilac.somatic.vcf.gz"
    log:
        stdout = CFG["logs"]["lilac"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lilac.log"
    params:
        ref_genome_version = lambda wildcards: VERSION_MAP[wildcards.genome_build],
        resource_dir = lambda wildcards, input: os.path.dirname(input.resource_freq),
        output_dir = lambda wildcards, output: os.path.dirname(output.tsv),
        gene_copy_number_flag = lambda wildcards, input: f"-gene_copy_number {input.gene_copy_number[0]}" if input.gene_copy_number else "",
        somatic_vcf_flag = lambda wildcards, input: f"-somatic_vcf {input.somatic_vcf[0]}" if input.somatic_vcf else "",
        tunable_flags = LILAC_TUNABLE_FLAGS,
        # The lilac wrapper script defaults to -Xmx1g regardless of how much memory Snakemake/the
        # cluster scheduler actually reserved for this job -- those are two entirely separate
        # settings, and this rule never told the JVM about the real allocation. Real, fatal
        # failure without this: OutOfMemoryError while loading the ~24MB nucleotide reference file
        # alone, well before touching any BAM data. Mirrors modules/hmftools/1.1's own jvmheap
        # pattern exactly (used there for AMBER/COBALT/PURPLE/LINX, the same JVM-based tool
        # family) -- 80% of the reserved memory, leaving headroom for JVM overhead beyond the heap.
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["lilac"]
    container:
        None
    threads:
        CFG["threads"]["lilac"]
    resources:
        **CFG["resources"]["lilac"]
    shell:
        op.as_one_line("""
        mkdir -p {params.output_dir} &&
        lilac
        -Xmx{params.jvmheap}m
        -sample {wildcards.tumour_id}
        -ref_genome {input.fasta}
        -ref_genome_version {params.ref_genome_version}
        -resource_dir {params.resource_dir}
        -reference_bam {input.normal_bam}
        -tumor_bam {input.tumour_bam}
        {params.gene_copy_number_flag}
        {params.somatic_vcf_flag}
        {params.tunable_flags}
        -threads {threads}
        -output_dir {params.output_dir}
        > {log.stdout} 2>&1 &&
        if [ ! -f {output.somatic_vcf} ]; then touch {output.somatic_vcf}; fi
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _lilac_output_tsv:
    input:
        tsv = str(rules._lilac_run.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "lilac/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lilac.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module = True)

rule _lilac_output_qc:
    input:
        qc = str(rules._lilac_run.output.qc)
    output:
        qc = CFG["dirs"]["outputs"] + "lilac/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lilac.qc.tsv"
    run:
        op.relative_symlink(input.qc, output.qc, in_module = True)

# Symlinked unconditionally, but will be a 0-byte file for any pair that had no somatic_vcf input
# (see _lilac_run's own touch-if-missing step above) -- an honest, empty representation rather
# than a missing target, consistent with how this module handles every other optional input.
rule _lilac_output_somatic_vcf:
    input:
        somatic_vcf = str(rules._lilac_run.output.somatic_vcf)
    output:
        somatic_vcf = CFG["dirs"]["outputs"] + "lilac/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lilac.somatic.vcf.gz"
    run:
        op.relative_symlink(input.somatic_vcf, output.somatic_vcf, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks. Uses
# CFG["paired_runs"] (narrowed to pair_status == "matched" above) so tumour samples without a
# matched germline WES sample are never requested as targets -- HLA typing needs the patient's own
# germline sample.
rule _lilac_all:
    input:
        expand(
            [
                str(rules._lilac_output_tsv.output.tsv),
                str(rules._lilac_output_qc.output.qc),
                str(rules._lilac_output_somatic_vcf.output.somatic_vcf)
            ],
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            tumour_id = CFG["paired_runs"]["tumour_sample_id"],
            normal_id = CFG["paired_runs"]["normal_sample_id"],
            pair_status = CFG["paired_runs"]["pair_status"]
        )


##### CLEANUP #####


op.cleanup_module(CFG)
