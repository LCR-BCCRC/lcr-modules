#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Li Song (liulab-dfci/TRUST4)
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
# `CFG` is a shortcut to `config["lcr-modules"]["trust4"]`
CFG = op.setup_module(
    name = "trust4",
    version = "1.0",
    subdirectories = ["inputs", "trust4_reference", "trust4", "outputs"],
)

# op.setup_module() already restricts CFG["samples"] to seq_types declared in pairing_config
# (mrna only, per this module's own config/default.yaml) -- every sample referenced below is
# already guaranteed to be real RNA-seq, no further filtering needed anywhere in this file.

# genome_build -> which of TRUST4's own pre-built, genome-coordinate-aware V/D/J/C gene FASTAs
# to download/use. Mirrors modules/pvactools/1.0's own VEP_ASSEMBLY_MAP pattern exactly (a
# per-module copy of this mapping, not a cross-module reference -- see that module's own
# pvactools.smk for the identical idiom). Fail fast at parse time (not a confusing TRUST4
# runtime error) if any genome_build in the samples table isn't covered.
TRUST4_BCRTCR_MAP = CFG["options"]["trust4_bcrtcr_map"]
_possible_genome_builds = ", ".join(TRUST4_BCRTCR_MAP.keys())
for _genome_build in CFG["samples"]["genome_build"]:
    assert _genome_build in TRUST4_BCRTCR_MAP, (
        f"Samples table includes genome build '{_genome_build}', not yet covered by "
        f"options.trust4_bcrtcr_map (currently: {_possible_genome_builds}). Add it there before rerunning."
    )

# Pinned TRUST4 GitHub commit its own reference downloads are fetched from -- see
# options.trust4_repo_commit's own comment in default.yaml for the real reproducibility
# reasoning.
TRUST4_REPO_COMMIT = CFG["options"]["trust4_repo_commit"]


# Define rules to be run locally when using a compute cluster. Beyond the trivial symlink/
# target rules, this also includes the two reference-download rules -- small, quick network
# I/O, not worth cluster-job scheduling overhead (mirrors modules/mhc_hammer/1.0's own
# _mhc_hammer_download_reference being listed here for the same reason).
localrules:
    _trust4_download_bcrtcr,
    _trust4_download_imgt,
    _trust4_output_report,
    _trust4_output_cdr3,
    _trust4_output_annot,
    _trust4_output_airr,
    _trust4_output_report_reconstructed_fastq,
    _trust4_output_cdr3_reconstructed_fastq,
    _trust4_output_annot_reconstructed_fastq,
    _trust4_output_airr_reconstructed_fastq,
    _trust4_all,


##### RULES #####


# Normalizes the input into a real, usable BAM under '00-inputs/'. Real, confirmed production
# failure (2026-09): a sample's own inputs.sample_bam was actually CRAM content stored/
# symlinked with a ".bam" name (the same real-world storage convention modules/mhc_hammer/1.0's
# own _mhc_hammer_input_bam already documents hitting for DNA BAMs), and TRUST4's own
# bam-extractor crashed with "[bam_header_read] invalid BAM binary header (this is not a BAM
# file)". Confirmed via the Makefile that bam-extractor links against the *bundled*, ancient
# samtools-0.1.19's own libbam.a (`-lbam`), not htslib -- that library predates CRAM entirely,
# so unlike mhc_hammer's own downstream tools (all real htslib-based `samtools`, which content-
# sniffs correctly regardless of extension -- mhc_hammer's own fix was purely about index-
# sidecar naming), no amount of correct sidecar naming fixes this: bam-extractor structurally
# cannot parse CRAM at all, full stop.
#
# Content-sniffed by reading the real first 4 bytes directly (same technique already used in
# this module's own reclaim/diagnostic tooling, not htsfile -- avoids adding a runtime
# dependency just for a 4-byte check): real BAM/BGZF starts with the gzip magic 1f 8b; real
# CRAM starts with the literal ASCII "CRAM", regardless of what the file's own name/extension
# claims. Genuine BAM is symlinked (cheap -- no reason to pay a full re-encode when the input
# was already fine); genuine CRAM is decoded via `samtools view -b -T {input.fasta}` (a real
# htslib tool, correct regardless of what upstream naming convention produced it). `-T` is
# resolved via reference_files() for this exact {genome_build}, matching
# modules/mhc_hammer/1.0's own _mhc_hammer_vep_annotate's identical idiom. Deliberately NOT
# relying on CRAM's own automatic REF_CACHE/URL-based reference resolution: mhc_hammer's own
# CHANGELOG already documents that silently degrading into fetching reference sequences one at
# a time from EBI's ENA CRAM registry over the network on a real cluster run -- passing -T
# explicitly avoids that failure mode entirely rather than risking a repeat of it. Anything
# that's neither real BAM nor real CRAM magic fails loudly here with a clear message, rather
# than silently symlinking a corrupted/truncated/unrecognized file and letting it surface later
# as a confusing bam-extractor crash instead.
#
# Real, requested fix (an earlier version of this rule always re-encoded unconditionally,
# discarding the "already a real BAM" case's own free symlink option, and never temp()-marked
# its own output at all): both output.bam/output.bai are temp() -- _trust4_run is their only
# consumer, so a genuine re-encoded duplicate doesn't sit on disk forever once TRUST4 itself is
# done with it (and for the symlink case, temp() cleanup is free regardless). The ORIGINAL
# input.bam this rule reads from is never written to, moved, or deleted either way -- only
# what lands in 00-inputs/ (symlink or real converted copy) is temp()-managed.
#
# Index output named "{sample_id}.bai", NOT "{sample_id}.bam.bai": a real, confirmed
# AmbiguousRuleException (caught via a real dry-run of demo/run_trust4.smk) -- modules/utils/2.1's
# own _utils_bam_index rule has a completely generic "{out_dir}/{prefix}/{suffix}.bam.bai"
# output pattern with no constraints, so a *.bam.bai output here collides with it the moment any
# Snakefile including both this module and modules/utils needs to build this file (i.e. any real
# trust4 deployment chained to star, since star itself requires utils). Same class of bug
# modules/mhc_hammer/1.0's own CHANGELOG documents for _mhc_hammer_novoalign_postprocess's bai,
# and that module's own conclusion (a ruleorder guard is unreliable, no reliable way to detect
# whether an arbitrary outer Snakefile happens to load modules/utils) applies equally here.
# htslib/samtools accept BOTH "{bam}.bam.bai" and "{bam}.bai" as valid index-sidecar names for a
# file named "{bam}.bam", so this avoids the collision with zero functional downside.
# `samtools index` itself only knows how to write the "{bam}.bai" convention (no explicit -o
# output path in samtools 1.9), so the index is generated under that name first, then renamed.
# Indexing runs unconditionally in both branches (even the cheap symlink one) rather than trying
# to reuse an existing sidecar from the original input -- a fresh `samtools index` on an
# already-sorted BAM is a cheap block-boundary scan, not a full re-encode, so it isn't worth the
# complexity of tracking which sidecar-naming convention (.bam.bai vs .crai) the original might
# already have.
#
# No longer a localrule: the CRAM branch does real, potentially expensive work (a full linear
# re-encode) -- the symlink branch alone would still be cheap enough to run locally, but which
# branch a given sample needs isn't known until this rule actually runs.
rule _trust4_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"),
        bai = temp(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai")
    log:
        stdout = CFG["logs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}/normalize_bam.log"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    threads:
        CFG["threads"]["input_bam"]
    resources:
        **CFG["resources"]["input_bam"]
    # NOTE: this is a shell: block, not run:, even though the logic below is a real if/elif/else
    # -- Snakemake refuses conda:/container: on a run:-based rule ("Conda environments are only
    # allowed with shell, script, notebook, or wrapper directives"), a real constraint caught via
    # a real dry-run (an earlier version of this rule used run: + op.absolute_symlink() for the
    # symlink branch, which is not expressible here as a result). Magic-byte check done with
    # plain `od`/`tr` (POSIX/coreutils, always available) rather than htsfile, to avoid adding a
    # runtime dependency just for a 4-byte check -- matches the same raw-byte technique already
    # used in this module's own reclaim/diagnostic tooling. `op.as_one_line()` only space-joins
    # lines (never inserts `;`), so every if/elif/else/fi keyword below has an explicit trailing
    # `;` on the line before it.
    shell:
        op.as_one_line("""
        (
        magic=$(head -c4 {input.bam} | od -An -tx1 | tr -d ' ') ;
        if [[ "$magic" == 1f8b* ]] ;
        then
        ln -sf $(readlink -f {input.bam}) {output.bam} ;
        elif [[ "$magic" == "4352414d" ]] ;
        then
        samtools view -b -T {input.fasta} -@ {threads} -o {output.bam} {input.bam} ;
        else
        echo "ERROR: {input.bam}: first 4 bytes are neither real BAM/BGZF magic (1f8b) nor real CRAM magic (4352414d = ASCII 'CRAM') -- got $magic. This file may be corrupted/truncated, or in a format this rule does not recognize." >&2 ;
        exit 1 ;
        fi &&
        samtools index -@ {threads} {output.bam} &&
        mv {output.bam}.bai {output.bai}
        ) > {log.stdout} 2>&1
        """)


# Downloads TRUST4's own pre-built, genome-coordinate-aware V/D/J/C gene FASTA for this
# genome_build, pinned to a specific commit (options.trust4_repo_commit) for reproducibility.
# Cohort-wide per genome_build, built once and reused by every sample of that build. Mirrors
# modules/mhc_hammer/1.0's own _mhc_hammer_download_reference in spirit (bare wget, no
# checksum verification -- this repo's established convention for a small, pinned-commit
# reference file; see also modules/vcf2maf/1.3/vcf2maf.smk's own _vcf2maf_gamblr_config for
# the identical pinned-commit-download pattern).
rule _trust4_download_bcrtcr:
    output:
        bcrtcr = CFG["dirs"]["trust4_reference"] + "{genome_build}/bcrtcr.fa"
    log:
        stdout = CFG["logs"]["trust4_reference"] + "{genome_build}/download_bcrtcr.log"
    params:
        url = lambda wildcards: (
            f"https://raw.githubusercontent.com/liulab-dfci/TRUST4/{TRUST4_REPO_COMMIT}/"
            f"{TRUST4_BCRTCR_MAP[wildcards.genome_build]}"
        )
    conda:
        CFG["conda_envs"]["trust4_download"]
    container:
        CFG["container_envs"]["trust4_download"]
    threads:
        CFG["threads"]["download_bcrtcr"]
    resources:
        **CFG["resources"]["download_bcrtcr"]
    shell:
        op.as_one_line("""
        wget -qO {output.bcrtcr} "{params.url}"
        > {log.stdout} 2>&1
        """)


# Downloads TRUST4's own pre-built human IMGT+C.fa reference (species-only, not genome-build-
# specific), pinned to the same commit as above. Cohort-wide, no wildcards at all -- built
# once and reused by every sample regardless of genome_build.
rule _trust4_download_imgt:
    output:
        imgt = CFG["dirs"]["trust4_reference"] + "human_IMGT+C.fa"
    log:
        stdout = CFG["logs"]["trust4_reference"] + "download_imgt.log"
    params:
        url = f"https://raw.githubusercontent.com/liulab-dfci/TRUST4/{TRUST4_REPO_COMMIT}/human_IMGT+C.fa"
    conda:
        CFG["conda_envs"]["trust4_download"]
    container:
        CFG["container_envs"]["trust4_download"]
    threads:
        CFG["threads"]["download_imgt"]
    resources:
        **CFG["resources"]["download_imgt"]
    shell:
        op.as_one_line("""
        wget -qO {output.imgt} "{params.url}"
        > {log.stdout} 2>&1
        """)


# Runs TRUST4 on this sample's own STAR-aligned BAM. Real, confirmed CLI/behaviour
# (trust4=1.1.10, confirmed via a real local conda install + a real run against TRUST4's own
# bundled example.bam, 2026-09): run-trust4 internally calls
# bam-extractor -> trust4 (assembly) -> annotator -> trust-simplerep.pl/trust-airr.pl, and
# `--clean 1` deletes its own intermediate files (*_toassemble*, *_assembled_reads.fa,
# *_final.out, *_raw.out, *_airr_align.tsv) after a successful run -- confirmed directly by
# inspecting the output directory after a real run: exactly *_report.tsv/*_cdr3.out/
# *_annot.fa/*_airr.tsv survived, nothing else. This means none of those intermediate files
# need any Snakemake temp() bookkeeping at all -- they simply aren't declared as outputs here,
# since they won't exist by the time this rule's shell command exits successfully (unlike
# modules/mhc_hammer/1.0, which had to carefully audit and temp()-mark its own intermediates
# by hand since upstream never cleaned them up itself).
rule _trust4_run:
    input:
        bam = str(rules._trust4_input_bam.output.bam),
        bai = str(rules._trust4_input_bam.output.bai),
        bcrtcr = str(rules._trust4_download_bcrtcr.output.bcrtcr),
        imgt = str(rules._trust4_download_imgt.output.imgt)
    output:
        report = CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_report.tsv",
        cdr3 = CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_cdr3.out",
        annot = CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_annot.fa",
        airr = CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_airr.tsv"
    log:
        stdout = CFG["logs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/trust4_run.log"
    params:
        # run-trust4 combines --od + -o internally (prefix = "$outputDirectory/$prefix"), so
        # passing the sample's own subdirectory as --od and a bare sample_id as -o produces
        # exactly {outdir}/{sample_id}_report.tsv etc., matching the declared outputs above.
        # run-trust4 creates --od itself if missing (confirmed from its own source), no mkdir
        # needed here.
        outdir = lambda wildcards, output: os.path.dirname(output.report),
        # --repseq: see options.repseq's own comment in default.yaml for why this defaults off.
        repseq_flag = "--repseq" if CFG["options"]["repseq"] else "",
        # --abnormalUnmapFlag: see options.abnormal_unmap_flag's own comment in default.yaml.
        abnormal_unmap_flag = "--abnormalUnmapFlag" if CFG["options"]["abnormal_unmap_flag"] else ""
    conda:
        CFG["conda_envs"]["trust4"]
    container:
        CFG["container_envs"]["trust4"]
    threads:
        CFG["threads"]["trust4_run"]
    resources:
        **CFG["resources"]["trust4_run"]
    shell:
        op.as_one_line("""
        run-trust4
        -b {input.bam}
        -f {input.bcrtcr}
        --ref {input.imgt}
        -o {wildcards.sample_id}
        --od {params.outdir}
        -t {threads}
        --clean 1
        {params.repseq_flag}
        {params.abnormal_unmap_flag}
        > {log.stdout} 2>&1
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _trust4_output_report:
    input:
        report = str(rules._trust4_run.output.report)
    output:
        report = CFG["dirs"]["outputs"] + "report/{seq_type}--{genome_build}/{sample_id}.trust4_report.tsv"
    run:
        op.relative_symlink(input.report, output.report, in_module = True)


rule _trust4_output_cdr3:
    input:
        cdr3 = str(rules._trust4_run.output.cdr3)
    output:
        cdr3 = CFG["dirs"]["outputs"] + "cdr3/{seq_type}--{genome_build}/{sample_id}.trust4_cdr3.out"
    run:
        op.relative_symlink(input.cdr3, output.cdr3, in_module = True)


rule _trust4_output_annot:
    input:
        annot = str(rules._trust4_run.output.annot)
    output:
        annot = CFG["dirs"]["outputs"] + "annot/{seq_type}--{genome_build}/{sample_id}.trust4_annot.fa"
    run:
        op.relative_symlink(input.annot, output.annot, in_module = True)


rule _trust4_output_airr:
    input:
        airr = str(rules._trust4_run.output.airr)
    output:
        airr = CFG["dirs"]["outputs"] + "airr/{seq_type}--{genome_build}/{sample_id}.trust4_airr.tsv"
    run:
        op.relative_symlink(input.airr, output.airr, in_module = True)


##### OPT-IN: BAM-RECONSTRUCTED-FASTQ COMPARISON PATH #####
#
# Real, requested comparison feature, NOT part of this module's default target list
# (_trust4_all below deliberately does not expand() over these) -- runs TRUST4 in its own
# FASTQ mode (-1/-2, via fastq-extractor) instead of BAM mode (-b, via bam-extractor), on FASTQ
# reconstructed from the same already-normalized BAM _trust4_input_bam produces. This sidesteps
# the entire class of BAM/CRAM-convention issues this module has already hit twice (CRAM
# disguised as .bam; an aligner that doesn't place unmapped-pair mates adjacently) -- fastq-
# extractor is a completely different code path with no such assumptions -- at the real cost of
# not being a truly independent comparison: these FASTQs are reconstructed from an ALREADY-
# ALIGNED bam (via `samtools collate | samtools fastq`), so they inherit whatever alignment-time
# decisions already happened (secondary/supplementary alignment filtering, soft-clip handling),
# not the sample's true original sequencer output. Real cohort has no tracked original FASTQs to
# use instead (confirmed with the user) -- reconstruction from BAM is what's actually available.
#
# Deliberately opt-in rather than doubling every sample's compute by default: request a specific
# sample's own target directly, e.g.
#   snakemake ... results/trust4-1.0/99-outputs/report_reconstructed_fastq/{seqbuild}/{sample_id}.trust4_reconstructed_fastq_report.tsv
# "_reconstructed_fastq" is threaded through every rule/output name here deliberately, so
# provenance is unambiguous wherever one of these files ends up (never confusable with the
# primary, BAM-mode report).
#
# Verified end-to-end against TRUST4's own bundled example.bam before being wired in here: a
# real `samtools collate -O -u example.bam | samtools fastq -1 r1.fq -2 r2.fq ...` followed by a
# real `run-trust4 -f human_IMGT+C.fa --ref human_IMGT+C.fa -1 r1.fq -2 r2.fq` produced the same
# CDR3 calls/counts as the existing BAM-mode smoke test on that same file.
#
# Reuses the existing 02-trust4 numbered directory (nested under a reconstructed_fastq/
# subpath) rather than adding a new top-level subdirectory -- adding one to CFG["dirs"] would
# renumber every subdirectory after it (same reasoning already documented for why the
# completed-samples-style manifest pattern elsewhere in this repo avoids that).


# Reconstructs paired FASTQ from _trust4_input_bam's own (already CRAM-normalized, so this
# never has to deal with CRAM itself) BAM via `samtools collate` (group by read name, required
# for correct R1/R2 pairing) piped into `samtools fastq` (split into R1/R2; -0/-s discard
# unpaired/singleton reads rather than writing them anywhere, -n disables read-name /1 /2
# suffix mangling since TRUST4 doesn't need it). Not shared with _trust4_input_bam's own bam
# output's temp() lifetime concerns -- Snakemake already extends a temp() file's lifetime to
# cover every declared consumer, so adding this as a second consumer doesn't risk it being
# cleaned up before this rule gets to read it.
rule _trust4_bam_to_fastq:
    input:
        bam = str(rules._trust4_input_bam.output.bam)
    output:
        fq1 = temp(CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/reconstructed_fastq/{sample_id}_1.fq"),
        fq2 = temp(CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/reconstructed_fastq/{sample_id}_2.fq")
    log:
        stdout = CFG["logs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/reconstructed_fastq/bam_to_fastq.log"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    threads:
        CFG["threads"]["bam_to_fastq"]
    resources:
        **CFG["resources"]["bam_to_fastq"]
    shell:
        op.as_one_line("""
        (
        samtools collate -O -u -@ {threads} {input.bam} |
        samtools fastq -@ {threads} -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null -n -
        ) > {log.stdout} 2>&1
        """)


# Same shape as _trust4_run, but FASTQ input (-1/-2) instead of BAM (-b). Real, confirmed
# difference from run-trust4's own docs: FASTQ mode has no genome coordinates to extract
# candidate reads by, so -f uses human_IMGT+C.fa here (NOT the genome-coordinate bcrtcr.fa
# _trust4_run's own BAM-mode invocation needs) -- matches the exact example command in
# run-trust4's own README ("./run-trust4 -f human_IMGT+C.fa --ref human_IMGT+C.fa -1 ... -2
# ..."). options.repseq/options.abnormal_unmap_flag still apply here -- both are orthogonal to
# which input mode is used.
rule _trust4_run_fastq:
    input:
        fq1 = str(rules._trust4_bam_to_fastq.output.fq1),
        fq2 = str(rules._trust4_bam_to_fastq.output.fq2),
        imgt = str(rules._trust4_download_imgt.output.imgt)
    output:
        report = CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/reconstructed_fastq/{sample_id}_report.tsv",
        cdr3 = CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/reconstructed_fastq/{sample_id}_cdr3.out",
        annot = CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/reconstructed_fastq/{sample_id}_annot.fa",
        airr = CFG["dirs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/reconstructed_fastq/{sample_id}_airr.tsv"
    log:
        stdout = CFG["logs"]["trust4"] + "{seq_type}--{genome_build}/{sample_id}/reconstructed_fastq/trust4_run.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.report),
        repseq_flag = "--repseq" if CFG["options"]["repseq"] else "",
        abnormal_unmap_flag = "--abnormalUnmapFlag" if CFG["options"]["abnormal_unmap_flag"] else ""
    conda:
        CFG["conda_envs"]["trust4"]
    container:
        CFG["container_envs"]["trust4"]
    threads:
        CFG["threads"]["trust4_run_fastq"]
    resources:
        **CFG["resources"]["trust4_run_fastq"]
    shell:
        op.as_one_line("""
        run-trust4
        -f {input.imgt}
        --ref {input.imgt}
        -1 {input.fq1}
        -2 {input.fq2}
        -o {wildcards.sample_id}
        --od {params.outdir}
        -t {threads}
        --clean 1
        {params.repseq_flag}
        {params.abnormal_unmap_flag}
        > {log.stdout} 2>&1
        """)


rule _trust4_output_report_reconstructed_fastq:
    input:
        report = str(rules._trust4_run_fastq.output.report)
    output:
        report = CFG["dirs"]["outputs"] + "report_reconstructed_fastq/{seq_type}--{genome_build}/{sample_id}.trust4_reconstructed_fastq_report.tsv"
    run:
        op.relative_symlink(input.report, output.report, in_module = True)


rule _trust4_output_cdr3_reconstructed_fastq:
    input:
        cdr3 = str(rules._trust4_run_fastq.output.cdr3)
    output:
        cdr3 = CFG["dirs"]["outputs"] + "cdr3_reconstructed_fastq/{seq_type}--{genome_build}/{sample_id}.trust4_reconstructed_fastq_cdr3.out"
    run:
        op.relative_symlink(input.cdr3, output.cdr3, in_module = True)


rule _trust4_output_annot_reconstructed_fastq:
    input:
        annot = str(rules._trust4_run_fastq.output.annot)
    output:
        annot = CFG["dirs"]["outputs"] + "annot_reconstructed_fastq/{seq_type}--{genome_build}/{sample_id}.trust4_reconstructed_fastq_annot.fa"
    run:
        op.relative_symlink(input.annot, output.annot, in_module = True)


rule _trust4_output_airr_reconstructed_fastq:
    input:
        airr = str(rules._trust4_run_fastq.output.airr)
    output:
        airr = CFG["dirs"]["outputs"] + "airr_reconstructed_fastq/{seq_type}--{genome_build}/{sample_id}.trust4_reconstructed_fastq_airr.tsv"
    run:
        op.relative_symlink(input.airr, output.airr, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _trust4_all:
    input:
        expand(
            [
                str(rules._trust4_output_report.output.report),
                str(rules._trust4_output_cdr3.output.cdr3),
                str(rules._trust4_output_annot.output.annot),
                str(rules._trust4_output_airr.output.airr)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"],
            sample_id = CFG["samples"]["sample_id"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
