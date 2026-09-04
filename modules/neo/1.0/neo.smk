#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Hartwig Medical Foundation (hartwigmedical/hmftools, NEO)
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
    name = "neo",
    version = "1.0",
    subdirectories = ["inputs", "reference", "pave", "epitopes", "scores", "outputs"],
)

# Neoepitope prediction genuinely needs the patient's own germline sample (for the somatic
# mutation/HLA typing chain this module consumes) -- same reasoning as modules/lilac/1.0.
CFG["paired_runs"] = CFG["paired_runs"][CFG["paired_runs"]["pair_status"] == "matched"]

# genome_build -> NEO's own -ref_genome_version argument ("37"/"38"), mirrors
# modules/lilac/1.0's own options.version_map exactly -- deliberately a separate copy, not a
# cross-module reference, matching how every module in this HMF Tools family keeps its own copy.
VERSION_MAP = CFG["options"]["version_map"]

possible_genome_builds = ", ".join(list(VERSION_MAP.keys()))
for genome_build in CFG["paired_runs"]["tumour_genome_build"]:
    assert genome_build in VERSION_MAP, (
        f"Samples table includes genome build '{genome_build}', not yet covered by "
        f"options.version_map (currently: {possible_genome_builds}). Add it there before rerunning."
    )

# Mirrors modules/lilac/1.0's own use_masked_ref switch and masked_string construction exactly.
masked_string = ""
if CFG["options"]["use_masked_ref"]:
    masked_string = "_masked"

# Define rules to be run locally when using a compute cluster
localrules:
    _neo_download_reference,
    _neo_download_tpm_medians,
    _neo_download_ensembl_cache,
    _neo_prep_purple_dir,
    _neo_output_epitopes,
    _neo_output_neoepitope_scores,
    _neo_output_peptide_scores,
    _neo_all,


##### HELPER FUNCTIONS #####


# NOTE: every helper function below re-fetches `CFG` from the persistent `config` global (rather
# than closing over the module-level `CFG` variable) since op.cleanup_module(CFG) deletes the
# module-level shortcut before Snakemake lazily evaluates input/param functions.

# somatic_vcf/lilac_tsv_file/purple_purity_file (below) are REQUIRED (not optional) cross-module
# inputs -- unlike modules/lilac/1.0's own os.path.exists()-gated optional inputs, NEO cannot
# produce anything meaningful without real PURPLE/LILAC data, so there's no "run without it"
# fallback to support for those three. Each pattern points at one real, specific file (not a bare
# directory) so Snakemake's own path-matching recognises it as a real dependency on whatever rule
# actually produces it (in hmftools/1.1 or lilac, if that module is also included in the same
# runner) -- if it isn't, Snakemake correctly, honestly fails with "no rule to produce this file"
# rather than silently producing an empty/wrong result.

# linx_neo_epitope_file is the one exception -- OPTIONAL, existence-gated exactly like
# modules/lilac/1.0's own gene_copy_number/somatic_vcf. Confirmed on a real cluster run that no
# currently-released LINX version's own NeoEpitopeFusion schema includes the "copyNumber" field
# NEO v1.3's own reader unconditionally requires (checked every tagged LINX release through the
# current latest, v2.3 -- none of them write it; only the unreleased hmftools master branch does).
# Until a real LINX release closes this gap, this module skips fusion-derived neoepitopes rather
# than crash or fabricate the missing column -- see CHANGELOG.md for the full real-source trace.
# Point-mutation-derived neoepitopes (the PAVE-annotated somatic_vcf path, the majority for most
# samples) are entirely unaffected -- that's a separate code path in NeoEpitopeFinder.
def _neo_get_linx_epitope_file(wildcards):
    CFG = config["lcr-modules"]["neo"]
    pattern = CFG["inputs"].get("linx_neo_epitope_file", "")
    if not pattern:
        return []
    path = pattern.format(
        seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id, normal_id = wildcards.normal_id, pair_status = wildcards.pair_status
    )
    return [path] if os.path.isfile(path) else []

def _neo_get_raw_somatic_vcf(wildcards):
    CFG = config["lcr-modules"]["neo"]
    return CFG["inputs"]["somatic_vcf"].format(
        seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id, normal_id = wildcards.normal_id, pair_status = wildcards.pair_status
    )

# PAVE reads the somatic VCF via htsjdk, which requires a .tbi index physically alongside it --
# confirmed via a real crash (TribbleException: "An index is required, but none found.") on
# modules/hmftools/1.1's own symlinked purple_somatic_vcf output, which originally dropped PURPLE's
# own companion index on relocation. Tracked here as a real Snakemake input (not passed to any CLI
# flag -- PAVE finds it itself by the standard sibling-.tbi convention) purely so Snakemake enforces
# its existence/ordering before this rule runs, same path-pattern-matching convention as the VCF
# itself.
def _neo_get_raw_somatic_vcf_index(wildcards):
    return _neo_get_raw_somatic_vcf(wildcards) + ".tbi"

def _neo_get_lilac_tsv_file(wildcards):
    CFG = config["lcr-modules"]["neo"]
    return CFG["inputs"]["lilac_tsv_file"].format(
        seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id, normal_id = wildcards.normal_id, pair_status = wildcards.pair_status
    )

def _neo_get_purple_purity_file(wildcards):
    CFG = config["lcr-modules"]["neo"]
    return CFG["inputs"]["purple_purity_file"].format(
        seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id, normal_id = wildcards.normal_id, pair_status = wildcards.pair_status
    )

# PAVE's Ensembl gene/transcript cache (-ensembl_data_dir) -- mirrors modules/lilac/1.0's own
# _lilac_get_ensembl_data_dir_input exactly (same reuse-or-download pattern, same real subdirectory
# quirk: the downloaded zip unpacks into a build-versioned subdirectory, e.g. "37/", confirmed by
# lilac's own real-zip inspection earlier -- appended in params, not here).
def _neo_get_ensembl_data_dir_input(wildcards):
    CFG = config["lcr-modules"]["neo"]
    pattern = CFG["inputs"].get("ensembl_data_dir", "")
    if pattern:
        path = pattern.format(seq_type = wildcards.seq_type, genome_build = wildcards.genome_build)
        if os.path.isdir(path):
            return ancient(path)
    return str(rules._neo_download_ensembl_cache.output.cache).format(genome_build = wildcards.genome_build)


##### RULES #####


# NEO's own bundled scoring reference data (position-weight matrices, likelihood distributions).
# Real filename convention confirmed two ways: (1) reading BindCommon.formFilename() directly --
# "<dir>/neo_train_<fileType>.csv" when no -score_file_id is passed (ScoreConfig.ScoreFileId is an
# optional CLI item; a null id skips the "_<id>" suffix entirely and there is no "." before it --
# both details were wrong in this rule's first version), "<dir>/neo_train_<fileType>_<id>.csv" if
# one is; (2) the real, current hartwigmedical/hmftools pipeline/README_RESOURCES.md NEO section,
# which lists exactly these 8 filenames with no id suffix at all -- i.e. the standard, publicly
# distributed bundle is generated and used id-less. So this module does the same: no
# -score_file_id passed to NeoScorer (see that rule), matching the bundle exactly. Cohort-wide, no
# genome_build wildcard -- this is protein/peptide-level binding data, not genome-coordinate
# reference data, so (like modules/lilac/1.0's own 3 resource CSVs) it's build-independent.
rule _neo_download_reference:
    output:
        pos_weight = CFG["dirs"]["reference"] + "neo_train_pos_weight.csv",
        flank_pos_weight = CFG["dirs"]["reference"] + "neo_train_flank_pos_weight.csv",
        likelihood = CFG["dirs"]["reference"] + "neo_train_likelihood.csv",
        expression_dist = CFG["dirs"]["reference"] + "neo_train_expression_dist.csv",
        recognition = CFG["dirs"]["reference"] + "neo_train_recognition.csv",
        rand_dist = CFG["dirs"]["reference"] + "neo_train_rand_dist.csv",
        likelihood_rand_dist = CFG["dirs"]["reference"] + "neo_train_likelihood_rand_dist.csv",
        exp_likelihood_rand_dist = CFG["dirs"]["reference"] + "neo_train_exp_likelihood_rand_dist.csv"
    log:
        stdout = CFG["logs"]["reference"] + "download_reference.log"
    params:
        url = "https://www.bcgsc.ca/downloads/morinlab/hmftools-references/neo"
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
        wget -qO {output.pos_weight} {params.url}/$(basename {output.pos_weight}) &&
        wget -qO {output.flank_pos_weight} {params.url}/$(basename {output.flank_pos_weight}) &&
        wget -qO {output.likelihood} {params.url}/$(basename {output.likelihood}) &&
        wget -qO {output.expression_dist} {params.url}/$(basename {output.expression_dist}) &&
        wget -qO {output.recognition} {params.url}/$(basename {output.recognition}) &&
        wget -qO {output.rand_dist} {params.url}/$(basename {output.rand_dist}) &&
        wget -qO {output.likelihood_rand_dist} {params.url}/$(basename {output.likelihood_rand_dist}) &&
        wget -qO {output.exp_likelihood_rand_dist} {params.url}/$(basename {output.exp_likelihood_rand_dist})
        ) > {log.stdout} 2>&1
        """)


# Downloads NEO's own cohort TPM-by-cancer-type medians (HMF's real "hmf_3444" cohort --
# isofox.hmf_3444.transcript_medians.<build>.csv, from the same misc/neo/tpm_cohort/ location in
# HMF's own pipeline_resources bundle the 8 reference CSVs above come from). Genome-build-keyed,
# unlike those 8 -- this is real per-transcript expression data, not build-independent
# protein/peptide binding data. Real column headers confirmed directly: HMF's own broad
# tissue-of-origin categories (Stomach, Lung, Liver, Lymphoid tissue, Skin, ... Other), plus "All"
# for the pan-cancer median -- NOT per-disease-subtype, so options.cancer_type only needs one real
# value for this whole lymphoma/leukemia cohort: "Lymphoid tissue" (see its own comment in
# config/default.yaml). Used by _neo_scorer's own -cancer_tpm_medians_file/-cancer_type; without
# either, NeoScorer gracefully falls back to the pan-cancer "All" column rather than failing
# (confirmed in TpmMediansCache.java), so this rule is optional in the same spirit as the other
# reference downloads, not a hard requirement -- _neo_scorer still runs fine (just without
# expression-adjusted TPM columns) if it's skipped.
rule _neo_download_tpm_medians:
    output:
        tpm_medians = CFG["dirs"]["reference"] + "{genome_build}/neo_tpm_medians.csv"
    log:
        stdout = CFG["logs"]["reference"] + "{genome_build}/download_tpm_medians.log"
    params:
        url = "https://www.bcgsc.ca/downloads/morinlab/hmftools-references/neo",
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
        wget -qO {output.tpm_medians} {params.url}/neo_tpm_medians_{params.alt_build}.csv
        > {log.stdout} 2>&1
        """)


# Downloads PAVE's Ensembl gene/transcript cache, genome-build-keyed -- exact same
# www.bcgsc.ca/downloads/morinlab/hmftools-references/ensembl_data_cache/{build}.zip mirror and
# VERSION_MAP values modules/hmftools/1.1 and modules/lilac/1.0 both already use for this exact
# upstream reference data. Only actually in the DAG when inputs.ensembl_data_dir's reuse check
# (_neo_get_ensembl_data_dir_input) doesn't find an existing directory to point at instead.
rule _neo_download_ensembl_cache:
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


# Annotates the required somatic_vcf input with the VCF INFO tag ("IMPACT") NEO itself requires to
# identify each variant's gene/coding-effect -- the exact same real dependency already found and
# fixed once for modules/lilac/1.0 (see that module's own _lilac_pave_annotate), confirmed here
# independently by reading NEO's own source directly: PointMutationData.isRelevantMutation() reads
# impact.WorstCodingEffect straight from com.hartwig.hmftools.common.variant.impact.VariantImpact
# (the same PAVE-written tag), and NeoSampleTask.java checks variantContext.hasAttribute(
# VAR_TRANS_IMPACT_ANNOTATION) before ever considering a variant. Without this, every variant would
# silently be excluded from neoepitope prediction, with no error -- LILAC's own README warning
# about this applies equally here. This module owns its own copy of this step rather than
# cross-referencing lilac's own PAVE output, consistent with how every module in this HMF Tools
# family keeps its own copy of shared logic (version_map, reference downloads, etc.) -- also
# correct regardless of whether modules/lilac/1.0 is even included in a given runner.
rule _neo_pave_annotate:
    input:
        vcf = _neo_get_raw_somatic_vcf,
        vcf_index = _neo_get_raw_somatic_vcf_index,
        ensembl_cache = _neo_get_ensembl_data_dir_input,
        fasta = reference_files("genomes/{genome_build}" + masked_string + "/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["pave"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.pave.vcf.gz"
    log:
        stdout = CFG["logs"]["pave"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pave.log"
    params:
        ref_genome_version = lambda wildcards: VERSION_MAP[wildcards.genome_build],
        # The real, downloaded ensembl_data_cache zip unpacks into a build-versioned subdirectory
        # (confirmed by downloading the real 37.zip and inspecting it directly while building
        # modules/lilac/1.0 -- its files are all under "37/", not at the top level) -- true whether
        # this module downloaded the cache itself or it's an externally-supplied one following the
        # same modules/hmftools/1.1-established layout.
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
        -input_vcf {input.vcf}
        -ref_genome {input.fasta}
        -ref_genome_version {params.ref_genome_version}
        -ensembl_data_dir {params.ensembl_dir}
        -output_dir {params.output_dir}
        -output_vcf {output.vcf}
        -threads {threads}
        > {log.stdout} 2>&1
        """)


# Step 1: identifies candidate neoepitopes from point mutations (the PAVE-annotated somatic VCF)
# and, when available, gene fusions (LINX's own neoepitope file). Real entry point confirmed by
# downloading the actual neo-v1.3 release jar and inspecting its class listing directly: the real
# class is com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder -- the module's own README names a
# stale, nonexistent com.hartwig.hmftools.neo.cohort.NeoEpitopeAnnotator instead. Real required
# config (-linx_dir/-somatic_vcf/-ref_genome_version/-ref_genome) confirmed via the real jar's own
# registered-config-item dump -- -linx_dir must still always be passed (CLI-required), but confirmed
# via NeoSampleTask.getSvFusions() source that it gracefully warns and skips fusion-neoepitope
# lookup (rather than crashing) when the sample's own expected file isn't present in that
# directory -- exploited deliberately below (params.linx_dir) when linx_epitope_file isn't
# configured/found, rather than needing to omit the flag itself. RNA/expression annotation
# deliberately out of scope for v1 (DNA-only mutation predictions only) -- same precedent already
# set for modules/mhc_hammer/1.0 and modules/lilac/1.0.
rule _neo_epitope_finder:
    input:
        linx_epitope_file = _neo_get_linx_epitope_file,
        somatic_vcf = str(rules._neo_pave_annotate.output.vcf),
        fasta = reference_files("genomes/{genome_build}" + masked_string + "/genome_fasta/genome.fa"),
        ensembl_cache = _neo_get_ensembl_data_dir_input
    output:
        neo_data = CFG["dirs"]["epitopes"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.neo.neo_data.tsv"
    log:
        stdout = CFG["logs"]["epitopes"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/epitope_finder.log"
    params:
        ref_genome_version = lambda wildcards: VERSION_MAP[wildcards.genome_build],
        # Falls back to this rule's own output_dir (guaranteed to exist, guaranteed to never
        # contain a real <sample>.linx.neoepitope.tsv) when linx_epitope_file isn't
        # configured/found -- satisfies -linx_dir's CLI-required check while NeoEpitopeFinder
        # gracefully finds no fusion-neoepitope file there and skips that part of its own analysis.
        linx_dir = lambda wildcards, input, output: os.path.dirname(input.linx_epitope_file[0]) if input.linx_epitope_file else os.path.dirname(output.neo_data),
        ensembl_dir = lambda wildcards, input: os.path.join(str(input.ensembl_cache), VERSION_MAP[wildcards.genome_build]),
        output_dir = lambda wildcards, output: os.path.dirname(output.neo_data),
        req_amino_acids = CFG["options"]["req_amino_acids"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["neo"]
    container:
        None
    threads:
        CFG["threads"]["epitope_finder"]
    resources:
        **CFG["resources"]["epitope_finder"]
    shell:
        op.as_one_line("""
        mkdir -p {params.output_dir} &&
        neo -Xmx{params.jvmheap}m
        -sample {wildcards.tumour_id}
        -ref_genome {input.fasta}
        -ref_genome_version {params.ref_genome_version}
        -ensembl_data_dir {params.ensembl_dir}
        -linx_dir {params.linx_dir}
        -somatic_vcf {input.somatic_vcf}
        -req_amino_acids {params.req_amino_acids}
        -threads {threads}
        -output_dir {params.output_dir}
        > {log.stdout} 2>&1
        """)


# NeoScorer's own -purple_dir reads both PURPLE's real, unmodified <sample>.purple.qc (via
# PurpleQCFile, whose own reader uses defaulted getValue() lookups -- confirmed NOT to hit the same
# crash class below) and <sample>.purple.purity.tsv (via PurplePurity, whose reader does an
# unchecked fieldsIndexMap.get() per column). Confirmed via a real crash
# (NullPointerException reading fieldsIndexMap.get("runMode")) that NEO v1.3's own PurplePurity
# reader (checked out at the real neo-v1.3 tag) requires two columns -- "runMode"/"targeted" -- that
# modules/hmftools/1.1's pinned PURPLE 2.54 simply never wrote (confirmed by checking out the real
# purple-v2.54 tag's own writer, PurityContextFile.java, directly: no such columns exist in its own
# header() at all). A large real upstream version gap -- PURPLE has since reached v4.4 -- but unlike
# the LINX neoepitope copyNumber gap (real ADR discussed with the user), both missing values here
# are deterministic run metadata this module already knows, not unknowable biological measurements,
# so they're filled in directly rather than skipped: runMode is always TUMOR_GERMLINE (this module's
# own CFG["paired_runs"] is narrowed to pair_status == "matched" only, so every sample NeoScorer
# ever processes really is a matched tumour/normal run), and targeted is read straight off
# {seq_type} (true for capture/WES, false for genome/WGS) -- both real, already-known facts, not
# guesses. Writes a small "purple_compat" directory (real purity.tsv + 2 appended columns, real QC
# file symlinked in unmodified) for _neo_scorer's own -purple_dir to point at instead of
# hmftools/1.1's real PURPLE directory directly.
rule _neo_prep_purple_dir:
    input:
        purity = _neo_get_purple_purity_file
    output:
        purity = CFG["dirs"]["scores"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/purple_compat/{tumour_id}.purple.purity.tsv",
        qc = CFG["dirs"]["scores"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/purple_compat/{tumour_id}.purple.qc"
    params:
        real_qc = lambda wildcards, input: os.path.join(os.path.dirname(input.purity), f"{wildcards.tumour_id}.purple.qc"),
        run_mode = "TUMOR_GERMLINE",
        targeted = lambda wildcards: "true" if wildcards.seq_type == "capture" else "false"
    run:
        with open(input.purity) as f:
            lines = f.read().rstrip("\n").split("\n")
        lines[0] += "\trunMode\ttargeted"
        for i in range(1, len(lines)):
            lines[i] += f"\t{params.run_mode}\t{params.targeted}"
        with open(output.purity, "w") as f:
            f.write("\n".join(lines) + "\n")
        op.relative_symlink(params.real_qc, output.qc, in_module = True)


# Step 2: scores each candidate neoepitope's peptides for MHC binding/presentation likelihood,
# using this pair's own LILAC allele typing and PURPLE purity/CN context. Real entry point
# confirmed the same way as _neo_epitope_finder: com.hartwig.hmftools.neo.score.NeoScorer (README
# names a stale com.hartwig.hmftools.neo.scorer.NeoScorer). Invoked via the installed `neo` wrapper
# script with the target class as its own first argument (confirmed directly: the wrapper inspects
# its first non-memory argument and switches from `-jar neo.jar` to `-cp neo.jar <class>` when it
# starts with "com.hartwig.hmftools.", so -Xmx placement doesn't matter). NeoScorer's own real
# registered config has zero CLI-flagged-REQUIRED items at all (confirmed via the real jar) -- its
# true requirements are enforced at runtime; -score_file_dir is confirmed functionally required via
# a real crash (NullPointerException in BindScoreMatrix.loadFromCsv) without it. -score_file_id is
# deliberately NOT passed -- ScoreConfig.ScoreFileId is optional, and a null id makes NeoScorer look
# for the plain, id-less filenames _neo_download_reference downloads (matching the standard,
# publicly distributed HMF resource bundle exactly -- see that rule's own comment). RNA/isofox/
# -cancer_tpm_medians_file deliberately omitted for v1 (README confirms graceful degradation to
# un-adjusted presentation likelihood without RNA).
rule _neo_scorer:
    input:
        neo_data = str(rules._neo_epitope_finder.output.neo_data),
        lilac_tsv = _neo_get_lilac_tsv_file,
        purple_purity = str(rules._neo_prep_purple_dir.output.purity),
        purple_qc = str(rules._neo_prep_purple_dir.output.qc),
        score_files = [
            str(rules._neo_download_reference.output.pos_weight),
            str(rules._neo_download_reference.output.flank_pos_weight),
            str(rules._neo_download_reference.output.likelihood),
            str(rules._neo_download_reference.output.expression_dist),
            str(rules._neo_download_reference.output.recognition),
            str(rules._neo_download_reference.output.rand_dist),
            str(rules._neo_download_reference.output.likelihood_rand_dist),
            str(rules._neo_download_reference.output.exp_likelihood_rand_dist)
        ],
        ensembl_cache = _neo_get_ensembl_data_dir_input,
        tpm_medians = str(rules._neo_download_tpm_medians.output.tpm_medians)
    output:
        neoepitopes = CFG["dirs"]["scores"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.neo.neoepitope.tsv",
        peptide_scores = CFG["dirs"]["scores"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.neo.peptide_scores.tsv"
    log:
        stdout = CFG["logs"]["scores"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/scorer.log"
    params:
        ensembl_dir = lambda wildcards, input: os.path.join(str(input.ensembl_cache), VERSION_MAP[wildcards.genome_build]),
        neo_dir = lambda wildcards, input: os.path.dirname(input.neo_data),
        lilac_dir = lambda wildcards, input: os.path.dirname(input.lilac_tsv),
        purple_dir = lambda wildcards, input: os.path.dirname(input.purple_purity),
        score_file_dir = lambda wildcards, input: os.path.dirname(input.score_files[0]),
        peptide_lengths = CFG["options"]["peptide_lengths"],
        cancer_type = CFG["options"]["cancer_type"],
        output_dir = lambda wildcards, output: os.path.dirname(output.neoepitopes),
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["neo"]
    container:
        None
    threads:
        CFG["threads"]["scorer"]
    resources:
        **CFG["resources"]["scorer"]
    shell:
        op.as_one_line("""
        mkdir -p {params.output_dir} &&
        neo com.hartwig.hmftools.neo.score.NeoScorer
        -Xmx{params.jvmheap}m
        -sample {wildcards.tumour_id}
        -ensembl_data_dir {params.ensembl_dir}
        -neo_dir {params.neo_dir}
        -lilac_dir {params.lilac_dir}
        -purple_dir {params.purple_dir}
        -score_file_dir {params.score_file_dir}
        -peptide_lengths {params.peptide_lengths}
        -cancer_tpm_medians_file {input.tpm_medians}
        -cancer_type "{params.cancer_type}"
        -write_peptide_scores
        -threads {threads}
        -output_dir {params.output_dir}
        > {log.stdout} 2>&1
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _neo_output_epitopes:
    input:
        neo_data = str(rules._neo_epitope_finder.output.neo_data)
    output:
        neo_data = CFG["dirs"]["outputs"] + "epitopes/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.neo.neo_data.tsv"
    run:
        op.relative_symlink(input.neo_data, output.neo_data, in_module = True)

rule _neo_output_neoepitope_scores:
    input:
        neoepitopes = str(rules._neo_scorer.output.neoepitopes)
    output:
        neoepitopes = CFG["dirs"]["outputs"] + "neoepitope_scores/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.neo.neoepitope.tsv"
    run:
        op.relative_symlink(input.neoepitopes, output.neoepitopes, in_module = True)

rule _neo_output_peptide_scores:
    input:
        peptide_scores = str(rules._neo_scorer.output.peptide_scores)
    output:
        peptide_scores = CFG["dirs"]["outputs"] + "peptide_scores/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.neo.peptide_scores.tsv"
    run:
        op.relative_symlink(input.peptide_scores, output.peptide_scores, in_module = True)


# Generates the target sentinels for each run
rule _neo_all:
    input:
        expand(
            [
                str(rules._neo_output_epitopes.output.neo_data),
                str(rules._neo_output_neoepitope_scores.output.neoepitopes),
                str(rules._neo_output_peptide_scores.output.peptide_scores)
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
