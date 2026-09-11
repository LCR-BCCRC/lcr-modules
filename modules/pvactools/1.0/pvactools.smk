#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Griffith Lab (griffithlab/pVACtools)
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####

import json
import os
import shlex
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
    name = "pvactools",
    version = "1.0",
    subdirectories = ["inputs", "coding_variants", "mhcflurry_reference", "vep_plugins", "vep", "alleles", "pvacseq", "outputs"],
)

# Neoantigen prediction needs a real tumour/normal somatic comparison -- same reasoning as
# modules/lilac/1.0 and modules/neo/1.0.
CFG["paired_runs"] = CFG["paired_runs"][CFG["paired_runs"]["pair_status"] == "matched"]

# Fail fast at parse time (not a confusing pvacseq runtime error) if any IEDB-gated algorithm is
# requested without a real local IEDB install -- see options.iedb_gated_algorithms' own comment in
# default.yaml for which algorithms this currently covers.
_requested_gated = set(CFG["options"]["prediction_algorithms"]) & set(CFG["options"]["iedb_gated_algorithms"])
assert not _requested_gated or CFG["options"]["iedb_install_directory"], (
    f"options.prediction_algorithms includes IEDB-gated algorithm(s) {_requested_gated} but "
    f"options.iedb_install_directory is not set. Either remove them or provide a real local IEDB "
    f"install directory (see options.iedb_install_directory's own comment in default.yaml)."
)

# genome_build -> VEP's own -assembly argument ("GRCh37"/"GRCh38"). Real, confirmed requirement
# (not just tidiness): a real run against a shared cache directory containing multiple assemblies
# under the same VEP cache version failed with "ERROR: Multiple assemblies found for cache version
# 105 (GRCh38, GRCh37) - specify one using --assembly [assembly]" -- VEP cannot always infer this on
# its own. Mirrors modules/vcf2maf/1.3's own VCF2MAF_VERSION_MAP exactly (same
# "grch37"/"grch38"-style genome_build values, same "grch" -> "GRCh" transform), kept as a separate,
# user-configurable copy rather than a cross-module reference, matching how every other module in
# this repo keeps its own copy of this same genome-build mapping.
VEP_ASSEMBLY_MAP = CFG["options"]["vep_assembly_map"]
_possible_genome_builds = ", ".join(VEP_ASSEMBLY_MAP.keys())
for _genome_build in CFG["paired_runs"]["tumour_genome_build"]:
    assert _genome_build in VEP_ASSEMBLY_MAP, (
        f"Samples table includes genome build '{_genome_build}', not yet covered by "
        f"options.vep_assembly_map (currently: {_possible_genome_builds}). Add it there before rerunning."
    )

# Single, cohort-wide harmonized build for _pvactools_join_neoantigen_maf's own output -- NOT a
# per-pair wildcard (confirmed with the user: one oncoplot wants one shared coordinate system, not
# multiple simultaneous projections). Lowercase "grch37"/"grch38", matching modules/vcf2maf/1.3's
# own internal genome "version" convention (VEP_ASSEMBLY_MAP above uses VEP's own "GRCh37"/"GRCh38"
# --assembly convention instead -- same information, different casing, for a different tool).
NEOANTIGEN_TARGET_BUILD = CFG["options"]["neoantigen_target_build"].lower()
assert NEOANTIGEN_TARGET_BUILD in ("grch37", "grch38"), (
    f"options.neoantigen_target_build must be 'grch37' or 'grch38' (got '{CFG['options']['neoantigen_target_build']}') "
    f"-- matches modules/vcf2maf/1.3's own CrossMap chain-pairing limitation (hg19<->hg38 only, see "
    f"that module's own vcf2maf.smk get_chain())."
)

# Define rules to be run locally when using a compute cluster. Beyond the trivial symlink/target
# rules, this also includes the small, pure Python bookkeeping steps (subsetting positions from an
# already-computed MAF, building the allele-list string) -- fast enough that submitting them as
# separate cluster jobs is pure scheduling overhead, same reasoning as modules/mhc_hammer/1.0's own
# localrules list.
localrules:
    _pvactools_download_mhcflurry_models,
    _pvactools_install_vep_plugins,
    _pvactools_extract_coding_regions,
    _pvactools_prepare_allele_list,
    _pvactools_output_class1_filtered,
    _pvactools_output_class1_all_epitopes,
    _pvactools_output_class2_filtered,
    _pvactools_output_class2_all_epitopes,
    _pvactools_output_combined_filtered,
    _pvactools_output_combined_all_epitopes,
    _pvactools_output_combined_aggregated,
    _pvactools_apply_additional_filter,
    _pvactools_join_neoantigen_maf,
    _pvactools_all,


##### HELPER FUNCTIONS #####


# NOTE: every helper function below re-fetches `CFG` from the persistent `config` global (rather
# than closing over the module-level `CFG` variable) since Snakemake evaluates input/param
# functions lazily, after op.cleanup_module(CFG) has already deleted the module-level name -- same
# pattern already established in modules/mhc_hammer/1.0/mhc_hammer.smk and modules/gridss/2.0.

def _pvactools_format_pair_path(pattern, wildcards):
    return pattern.format(
        seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id, normal_id = wildcards.normal_id, pair_status = wildcards.pair_status
    )

# Required. This module has no variant source of its own -- see inputs.vcf2maf_raw_vcf's own
# comment in default.yaml. No existence gate (unlike modules/lilac/1.0's optional cross-module
# inputs): a pair with no vcf2maf output simply cannot run pvacseq at all, so this is a hard
# dependency, mirroring modules/neo/1.0's own required-input helper pattern (_neo_get_raw_somatic_vcf
# etc.) rather than modules/lilac/1.0's optional one.
def _pvactools_get_vcf2maf_raw_vcf(wildcards):
    CFG = config["lcr-modules"]["pvactools"]
    return _pvactools_format_pair_path(CFG["inputs"]["vcf2maf_raw_vcf"], wildcards)

def _pvactools_get_vcf2maf_maf(wildcards):
    CFG = config["lcr-modules"]["pvactools"]
    return _pvactools_format_pair_path(CFG["inputs"]["vcf2maf_maf"], wildcards)

# Builds the concrete additional_filters/{preset}/... path for a specific preset, for a rule that
# isn't itself wildcarded on filter_preset (_pvactools_apply_additional_filter is). Mirrors that
# rule's own declared output pattern directly rather than trying to coerce rules.X.output.Y (whose
# string form still carries THAT rule's own {filter_preset} wildcard) into a second .format() call.
def _pvactools_get_filtered_for_preset(wildcards, preset):
    CFG = config["lcr-modules"]["pvactools"]
    pattern = (
        CFG["dirs"]["outputs"] +
        "additional_filters/" + preset +
        "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.tsv"
    )
    return _pvactools_format_pair_path(pattern, wildcards)

# Mirrors modules/vcf2maf/1.3's own get_chain() (vcf2maf.smk) exactly: same fixed hg19<->hg38
# chain-file pairing (that module's own documented limitation -- "This only currently supports
# hg38 and hg19" -- not a new one introduced here) and the same reference_files() mechanism this
# module already uses for the genome FASTA in _pvactools_vep_annotate. Reuses VEP_ASSEMBLY_MAP
# (already required, already validated at parse time) for "which build family is this" rather than
# a second, separate genome-build map -- same "parallel copy, not cross-module coupling" convention
# already established for VEP_ASSEMBLY_MAP itself. Returns None if this pair's native build already
# matches NEOANTIGEN_TARGET_BUILD (no lift needed).
def _pvactools_native_build_version(genome_build):
    return VEP_ASSEMBLY_MAP[genome_build].lower()

def _pvactools_get_crossmap_chain(genome_build):
    native_version = _pvactools_native_build_version(genome_build)
    if native_version == NEOANTIGEN_TARGET_BUILD:
        return None
    if native_version == "grch38" and NEOANTIGEN_TARGET_BUILD == "grch37":
        return reference_files("genomes/" + genome_build + "/chains/grch38/hg38ToHg19.over.chain")
    if native_version == "grch37" and NEOANTIGEN_TARGET_BUILD == "grch38":
        return reference_files("genomes/" + genome_build + "/chains/grch37/hg19ToHg38.over.chain")
    raise AssertionError(
        f"No supported CrossMap chain from '{native_version}' (genome_build '{genome_build}') to "
        f"'{NEOANTIGEN_TARGET_BUILD}' -- only grch37<->grch38 lifting is supported, mirroring "
        f"modules/vcf2maf/1.3's own same limitation."
    )

# Mirrors modules/mhc_hammer/1.0's own _mhc_hammer_get_patient_id_for_tumour exactly: recovers
# patient_id for a pair-level rule (wildcarded on tumour_id/normal_id, not patient_id) via
# CFG["runs"].
def _pvactools_get_patient_id_for_tumour(tumour_id, seq_type):
    CFG = config["lcr-modules"]["pvactools"]
    hits = op.filter_samples(CFG["runs"], tumour_sample_id = tumour_id, tumour_seq_type = seq_type)
    return hits["tumour_patient_id"].tolist()[0]

# Required. Class I HLA typing -- modules/mhc_hammer/1.0's own _mhc_hammer_output_hla_final_result
# (HLA-HD's own consolidated per-patient result). Patient-keyed like hla2_alleles below (HLA is
# germline, shared across a patient's samples). No existence gate INSIDE this function -- by the
# time it's called, _pvactools_all's own target-list filtering below has already guaranteed this
# file exists for any pair actually requested, mirroring modules/neo/1.0's own required-input
# helper pattern (a bare, ungated path formatter) rather than an inline check here.
def _pvactools_get_hla_final_result(wildcards):
    CFG = config["lcr-modules"]["pvactools"]
    pattern = CFG["inputs"]["hla_alleles"]
    patient_id = _pvactools_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type)
    return pattern.format(seq_type = wildcards.seq_type, genome_build = wildcards.genome_build, patient_id = patient_id)

# Optional -- Class II HLA typing. Existence-gated like modules/lilac/1.0's own gene_copy_number/
# somatic_vcf: a pair with no Class II typing available (mhc_hammer's own class II path wasn't run,
# or typing simply failed for this patient) still gets a Class-I-only pvacseq run rather than
# failing outright.
def _pvactools_get_hla2_alleles(wildcards):
    CFG = config["lcr-modules"]["pvactools"]
    pattern = CFG["inputs"].get("hla2_alleles", "")
    if not pattern:
        return []
    patient_id = _pvactools_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type)
    path = pattern.format(seq_type = wildcards.seq_type, genome_build = wildcards.genome_build, patient_id = patient_id)
    return [path] if os.path.isfile(path) else []

# Real, expected scenario for a cohort-wide sample list: some pairs will have failed mhc_hammer or
# the upstream SLMS-3/vcf2maf pipeline, so their required inputs (vcf2maf's raw VCF/MAF, mhc_hammer's
# own Class I typing) simply don't exist on disk yet. Without this filter, _pvactools_all would
# request targets for every matched pair regardless, and Snakemake would hard-fail with
# MissingInputException on the very first pair missing anything -- not a graceful per-pair skip.
# This isn't a job for checkpoints (those handle a rule whose own OUTPUT SET isn't known until it
# actually runs, e.g. splitting a file into an unknown number of pieces) -- the candidate pair list
# here is already fully known upfront from the sample table; the only question is whether each
# pair's inputs already exist, which is exactly the same existence-gated-optional-input pattern
# already used throughout this codebase (e.g. hla2_alleles above, modules/lilac/1.0's own
# gene_copy_number/somatic_vcf) -- just applied at the target-list level, since these three are
# hard requirements a pair can't run without at all, rather than ones it can gracefully run without.
def _pvactools_pair_has_required_inputs(row):
    wildcards_dict = dict(
        seq_type = row["tumour_seq_type"], genome_build = row["tumour_genome_build"],
        tumour_id = row["tumour_sample_id"], normal_id = row["normal_sample_id"], pair_status = row["pair_status"]
    )
    raw_vcf = CFG["inputs"]["vcf2maf_raw_vcf"].format(**wildcards_dict)
    maf = CFG["inputs"]["vcf2maf_maf"].format(**wildcards_dict)
    patient_id = _pvactools_get_patient_id_for_tumour(row["tumour_sample_id"], row["tumour_seq_type"])
    hla = CFG["inputs"]["hla_alleles"].format(
        seq_type = wildcards_dict["seq_type"], genome_build = wildcards_dict["genome_build"], patient_id = patient_id
    )
    return os.path.isfile(raw_vcf) and os.path.isfile(raw_vcf + ".tbi") and os.path.isfile(maf) and os.path.isfile(hla)

_pvactools_missing_inputs = ~CFG["paired_runs"].apply(_pvactools_pair_has_required_inputs, axis = 1)
if _pvactools_missing_inputs.any():
    print(
        f"INFO [pvactools]: skipping {_pvactools_missing_inputs.sum()} of {len(CFG['paired_runs'])} "
        f"pair(s) missing required upstream input (vcf2maf raw VCF/MAF or mhc_hammer Class I typing): "
        f"{CFG['paired_runs'][_pvactools_missing_inputs]['tumour_sample_id'].tolist()}"
    )
CFG["paired_runs"] = CFG["paired_runs"][~_pvactools_missing_inputs]


##### RULES #####


# Fetches MHCflurry's own trained models -- a one-time, cohort-wide download (not per-pair).
# MHCnuggetsI/MHCnuggetsII need no equivalent step: confirmed via a real local install that their
# trained models ship bundled directly in the pip/conda package (mhcnuggets/saves/production/*.h5).
#
# Deliberately does NOT set MHCFLURRY_DOWNLOADS_DIR at all -- confirmed via a real run, then traced
# through mhcflurry's own downloads.py source directly: this is an unconditional bug in mhcflurry's
# own configure(), not a usage mistake. configure() only sets the module-global _CURRENT_RELEASE
# (from --release/MHCFLURRY_DOWNLOADS_CURRENT_RELEASE/the package's own default) INSIDE its
# `if not MHCFLURRY_DOWNLOADS_DIR:` branch -- so setting MHCFLURRY_DOWNLOADS_DIR to ANY custom value
# skips that whole block, leaving _CURRENT_RELEASE permanently None regardless of --release. The
# `fetch` subcommand's own get_current_release_downloads() then does
# metadata['releases'][None] -> KeyError: None, unconditionally, every time. There is no
# environment-variable combination that has both a custom download directory and a working
# _CURRENT_RELEASE. Letting mhcflurry manage its own default location (a per-user
# `user_data_dir("mhcflurry", ...)`, e.g. under $HOME) is the only clean fix -- confirmed this
# doesn't break _pvactools_run's own later use of the downloaded models either, since
# get_path()/get_default_class1_models_dir() (what pvacseq actually calls at prediction time) only
# read get_downloads_dir() (the real, resolved directory), never _CURRENT_RELEASE. Flag: this means
# mhcflurry's models land in $HOME rather than this module's own results directory, unlike every
# other reference/resource this repo manages -- worth revisiting if $HOME quota is a real constraint
# on a given cluster.
rule _pvactools_download_mhcflurry_models:
    output:
        complete = touch(CFG["dirs"]["mhcflurry_reference"] + "mhcflurry_downloads.complete")
    log:
        stdout = CFG["logs"]["mhcflurry_reference"] + "download_mhcflurry_models.log"
    params:
        release = CFG["options"]["mhcflurry_release"]
    conda:
        CFG["conda_envs"]["pvactools"]
    container:
        None
    threads:
        CFG["threads"]["download_mhcflurry"]
    resources:
        **CFG["resources"]["download_mhcflurry"]
    shell:
        op.as_one_line("""
        mhcflurry-downloads fetch --release {params.release} > {log.stdout} 2>&1
        """)


# Installs pVACtools' own Frameshift/Wildtype VEP plugins -- a one-time, cohort-wide step (not
# per-pair). Genuinely new pattern in this repo (confirmed via repo-wide grep for
# "install_vep_plugin"/"VEP_plugins": zero existing hits) -- modeled on modules/lilac/1.0's own
# _lilac_download_ensembl_cache (module auto-fetches its own support data rather than requiring
# the user to supply it). Real, confirmed output (pvactools=7.1.3): exactly Frameshift.pm and
# Wildtype.pm, written flat into the given directory, no nesting.
rule _pvactools_install_vep_plugins:
    output:
        plugins_dir = directory(CFG["dirs"]["vep_plugins"]),
        frameshift = CFG["dirs"]["vep_plugins"] + "Frameshift.pm",
        wildtype = CFG["dirs"]["vep_plugins"] + "Wildtype.pm"
    log:
        stdout = CFG["logs"]["vep_plugins"] + "install_vep_plugins.log"
    conda:
        CFG["conda_envs"]["pvactools"]
    container:
        None
    threads:
        CFG["threads"]["install_vep_plugins"]
    resources:
        **CFG["resources"]["install_vep_plugins"]
    shell:
        op.as_one_line("""
        mkdir -p {output.plugins_dir} &&
        pvacseq install_vep_plugin {output.plugins_dir} > {log.stdout} 2>&1
        """)


# Extracts protein-altering variant positions from vcf2maf's own final MAF, reusing its
# already-computed Variant_Classification (from its earlier VEP run) rather than re-annotating and
# re-scoring variants already known to be silent/non-coding -- the large majority for a WGS sample.
# Position-only matching (+/-1bp padding, see src/maf_coding_positions.py), not exact ref/alt, to
# tolerate MAF/VCF indel-representation differences -- not yet independently verified to never drop
# a real candidate at a multiallelic/decomposed site (see this module's CHANGELOG). Pure stdlib
# Python bookkeeping -- localrule, no conda env needed.
rule _pvactools_extract_coding_regions:
    input:
        maf = _pvactools_get_vcf2maf_maf,
        script = CFG["options"]["maf_coding_positions_script"]
    output:
        regions = CFG["dirs"]["coding_variants"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.coding_regions.tsv",
        audit = CFG["dirs"]["outputs"] + "coding_variant_audit/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.coding_variant_audit.tsv"
    # Scheduling hint, not a resource change -- same real, requested pattern already established in
    # modules/mhc_hammer/1.0 (see that module's own _mhc_hammer_flagstat comment): biases Snakemake
    # toward advancing pairs already in flight through this module's own per-pair chain (especially
    # ones needing only the cheap final output/filter rules re-run, e.g. after a code change to
    # _pvactools_apply_additional_filter) before starting brand-new pairs' early, more expensive
    # stages. Values increase monotonically with DAG depth through this module's per-pair chain:
    # this rule/_pvactools_prepare_allele_list=10 (parallel branches, same depth), the grouped
    # "vep_and_run" job (_pvactools_subset_to_coding_variants/_pvactools_vep_annotate/_pvactools_run
    # -- all three share one priority value since a group job's own priority is the max of its
    # members' regardless)=20, the six 99-outputs/ symlink rules=30,
    # _pvactools_apply_additional_filter=40. Default priority (0) still applies to every other rule
    # (including the cohort-wide, one-time _pvactools_download_mhcflurry_models/
    # _pvactools_install_vep_plugins setup rules), so this whole per-pair chain already outranks
    # everything else too.
    priority: 10
    log:
        stdout = CFG["logs"]["coding_variants"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/extract_coding_regions.log"
    params:
        classifications = " ".join(CFG["options"]["coding_variant_classifications"])
    shell:
        op.as_one_line("""
        python3 {input.script}
        --maf {input.maf}
        --classifications {params.classifications}
        --output-regions {output.regions}
        --output-audit {output.audit}
        > {log.stdout} 2>&1
        """)


# Subsets vcf2maf's own raw, tool-agnostic VCF down to the coding-only regions extracted above, and
# reheaders the #CHROM sample columns to guarantee they equal wildcards.tumour_id/normal_id
# literally -- defensively handling both a raw caller VCF (its own original sample names) and a VCF
# already normalized to vcf2maf's own "TUMOR"/"NORMAL" convention (see modules/vcf2maf/1.3's own
# _vcf2maf_run, which does the reverse rename) -- pvacseq's own `sample_name`/
# `--normal-sample-name` arguments require an exact match to a real column name in the VCF's own
# #CHROM line. A raw VCF using some third, unrelated sample-naming convention would not be caught
# by this defensive check -- not yet independently verified against a real vcf2maf deployment (see
# this module's CHANGELOG).
#
# Output NOT permanent -- the CHANGELOG's own temp()-chain-checkpoint lesson from mhc_hammer: this
# step is cheap/fast (small variant count, no expensive computation), so temp()-marking it is safe
# (single consumer, _pvactools_vep_annotate) and avoids leaving a redundant permanent copy of a
# subsetted VCF around indefinitely.
#
# Real production crash fixed here (see this module's CHANGELOG for the full investigation):
# `bcftools annotate -x FORMAT/AF` strips the FORMAT/AF field unconditionally, and `bcftools norm
# -m-` splits multiallelic records into one biallelic record per ALT. Confirmed directly (real
# bcftools + vcfpy, not just reasoning from docs) that a caller can declare FORMAT/AF as
# `Number=1` in its own VCF header while still writing one comma-joined value per ALT allele on a
# multiallelic line (a real, non-compliant-but-real VCF pattern) -- vcfpy then can't parse it as a
# scalar float and keeps it as a raw string, which crashes pvacseq's own coverage_filter input
# conversion (`'>' not supported between instances of 'str' and 'int'`) downstream. `bcftools norm
# -m-` alone does NOT fix this: since it trusts the (wrong) Number=1 declaration, it just copies
# the whole unsplit string into both resulting biallelic records unchanged -- confirmed empirically.
# Dropping AF instead of trying to repair it lets pvacseq fall through to its own already-built-in,
# more robust AD+DP-based VAF calculation (`calculate_vaf()`: a plain `var_count / depth`, no
# scale/encoding ambiguity) -- AD (Number=R) and DP (Number=1) are far more universally,
# consistently declared across variant callers than AF, so this is a general robustness
# improvement for a deliberately caller-agnostic module, not a narrow patch for one caller's bug.
rule _pvactools_subset_to_coding_variants:
    input:
        vcf = _pvactools_get_vcf2maf_raw_vcf,
        vcf_tbi = lambda wildcards: _pvactools_get_vcf2maf_raw_vcf(wildcards) + ".tbi",
        regions = str(rules._pvactools_extract_coding_regions.output.regions)
    output:
        vcf = temp(CFG["dirs"]["coding_variants"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.coding.vcf.gz"),
        vcf_tbi = temp(CFG["dirs"]["coding_variants"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.coding.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["coding_variants"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/subset_to_coding_variants.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        None
    threads:
        CFG["threads"]["subset_coding_variants"]
    resources:
        **CFG["resources"]["subset_coding_variants"]
    # Bundled into one cluster job with _pvactools_vep_annotate/_pvactools_run -- see that rule's
    # own comment for why and the real resource-reservation tradeoff this brings.
    group:
        "vep_and_run"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 20
    shell:
        op.as_one_line("""
        (
        bcftools view -R {input.regions} {input.vcf} |
        bcftools annotate -x FORMAT/AF |
        bcftools norm -m- |
        awk -v t="{wildcards.tumour_id}" -v n="{wildcards.normal_id}"
        'BEGIN{{OFS="\\t"}} /^#CHROM/{{for(i=1;i<=NF;i++){{if($i=="TUMOR"||$i==t)$i=t; else if($i=="NORMAL"||$i==n)$i=n}}}} {{print}}' |
        bgzip -c > {output.vcf} &&
        tabix -p vcf {output.vcf}
        ) > {log.stdout} 2>&1
        """)


# Direct VEP call (not through any wrapper), applying pVACtools' own two required plugins.
# Mirrors modules/mhc_hammer/1.0's own _mhc_hammer_detect_muts convention exactly (options.vep_path
# exported onto PATH, container: None) -- VEP is user-supplied, not conda-managed, to avoid the
# bioconda/Perl solver conflicts already hit with ensembl-vep in the vcf2maf module. Required flags
# confirmed against pVACtools' own docs (version 7.x).
#
# temp()-marked for the same reason as _pvactools_subset_to_coding_variants above: single consumer
# (_pvactools_run), cheap/fast to regenerate (VEP annotation of an already-small, coding-only
# variant set is not the expensive step in this pipeline -- pvacseq's own multi-algorithm binding
# prediction is).
rule _pvactools_vep_annotate:
    input:
        vcf = str(rules._pvactools_subset_to_coding_variants.output.vcf),
        vcf_tbi = str(rules._pvactools_subset_to_coding_variants.output.vcf_tbi),
        # Deliberately the two real plugin FILES, not the bare plugins_dir directory: oncopipe's
        # own op.setup_module() pre-creates every declared subdirectory empty upfront, so a rule
        # depending on the bare directory alone would see it as "already existing" and never
        # actually schedule _pvactools_install_vep_plugins to populate it -- confirmed as a real
        # bug via a real dry-run before this was caught (see this module's CHANGELOG).
        vep_frameshift = str(rules._pvactools_install_vep_plugins.output.frameshift),
        vep_wildtype = str(rules._pvactools_install_vep_plugins.output.wildtype),
        vep_cache = ancient(CFG["inputs"]["vep_cache"]),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = temp(CFG["dirs"]["vep"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.pvacseq_input.vep.vcf.gz"),
        vcf_tbi = temp(CFG["dirs"]["vep"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.pvacseq_input.vep.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["vep"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/vep_annotate.log"
    params:
        vep_path = CFG["options"]["vep_path"],
        pick_flag = "--pick" if CFG["options"]["vep_pick"] else "",
        vep_plugins_dir = lambda wildcards, input: os.path.dirname(input.vep_frameshift),
        assembly = lambda wildcards: VEP_ASSEMBLY_MAP[wildcards.genome_build]
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        None
    threads:
        CFG["threads"]["vep_annotate"]
    resources:
        **CFG["resources"]["vep_annotate"]
    # Bundled into one cluster job with _pvactools_subset_to_coding_variants/_pvactools_run.
    # Real, requested change: VEP annotation of an already-small, coding-only variant set finishes
    # in well under a minute, and this cluster penalizes total job COUNT (short jobs still count
    # fully against that), not just aggregate compute time. group: bundles all three connected
    # per-pair rules into a single cluster submission -- confirmed this repo already uses this
    # convention (modules/freebayes/1.0's own "input_and_run"/"normalize_and_bgzip" groups).
    # Real tradeoff, same one already documented for modules/mhc_hammer/1.0's own grouping
    # consideration: Snakemake's group resource aggregation SUMS every member rule's mem_mb/threads
    # rather than assuming sequential execution, so the combined job reserves
    # 2000+8000+16000=26000 mem_mb and 1+4+8=13 threads for the group's whole duration (a few
    # minutes, dominated by _pvactools_run itself) -- not just _pvactools_run's own 16000/8. Revisit
    # this grouping (or split it back out) if that over-reservation becomes a real problem instead
    # of the job-count penalty it's meant to solve.
    group:
        "vep_and_run"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 20
    shell:
        op.as_one_line("""
        (
        export PATH={params.vep_path}:${{PATH}} &&
        vep -i {input.vcf}
        --format vcf --vcf -o STDOUT
        --symbol --terms SO --mane_select --canonical --tsl --biotype --hgvs
        --fasta {input.fasta}
        --offline --cache --dir_cache {input.vep_cache} --assembly {params.assembly}
        --plugin Frameshift --plugin Wildtype --dir_plugins {params.vep_plugins_dir}
        {params.pick_flag}
        --fork {threads}
        --force_overwrite |
        bgzip -c > {output.vcf} &&
        tabix -p vcf {output.vcf}
        ) > {log.stdout} 2>&1
        """)


# Builds pvacseq's own comma-separated HLA allele-list argument entirely from mhc_hammer's own
# HLA-HD-based typing: Class I (required) and Class II (optional). See
# src/prepare_allele_list.py's own header comment for the real, confirmed allele-naming conventions
# this implements, and for the caveats on Class I's own file format (not independently verified
# against a real HLA-HD final.result.txt in this session).
rule _pvactools_prepare_allele_list:
    input:
        hla_final_result = _pvactools_get_hla_final_result,
        hla2_alleles = _pvactools_get_hla2_alleles,
        script = CFG["options"]["prepare_allele_list_script"]
    output:
        allele_list = CFG["dirs"]["alleles"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.allele_list.txt",
        audit = CFG["dirs"]["outputs"] + "allele_summary/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.allele_summary.tsv"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority (this rule is a parallel branch at the same DAG depth as that one).
    priority: 10
    log:
        stdout = CFG["logs"]["alleles"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/prepare_allele_list.log"
    params:
        hla2_flag = lambda wildcards, input: f"--hla2-alleles {input.hla2_alleles[0]}" if input.hla2_alleles else ""
    shell:
        op.as_one_line("""
        python3 {input.script}
        --hla-final-result {input.hla_final_result}
        {params.hla2_flag}
        --output-allele-list {output.allele_list}
        --output-audit {output.audit}
        > {log.stdout} 2>&1
        """)


# Core rule: runs pVACseq's own binding-prediction pipeline. A single invocation mixes Class I +
# II alleles/algorithms together -- pvacseq automatically separates results into its own
# MHC_Class_I/MHC_Class_II output subdirectories (confirmed via a real installed pvactools=7.1.3).
#
# Real, confirmed CLI (pvacseq run --help, pvactools=7.1.3) -- in particular, there is NO short
# `-n` flag for the normal sample name, only `--normal-sample-name` (an earlier draft of this rule
# assumed one existed; caught before implementation via a real CLI check).
rule _pvactools_run:
    input:
        vep_vcf = str(rules._pvactools_vep_annotate.output.vcf),
        vep_vcf_tbi = str(rules._pvactools_vep_annotate.output.vcf_tbi),
        allele_list = str(rules._pvactools_prepare_allele_list.output.allele_list),
        mhcflurry_sentinel = str(rules._pvactools_download_mhcflurry_models.output.complete)
    output:
        outdir = directory(CFG["dirs"]["pvacseq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/"),
        complete = touch(CFG["dirs"]["pvacseq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pvacseq.complete")
    log:
        stdout = CFG["logs"]["pvacseq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pvacseq_run.log"
    params:
        alleles = lambda wildcards, input: open(input.allele_list).read().strip(),
        algorithms = " ".join(CFG["options"]["prediction_algorithms"]),
        e1 = CFG["options"]["class_i_epitope_lengths"],
        e2 = CFG["options"]["class_ii_epitope_lengths"],
        binding_threshold = CFG["options"]["binding_threshold"],
        top_score_metric = CFG["options"]["top_score_metric"],
        iedb_flag = f"--iedb-install-directory {CFG['options']['iedb_install_directory']}" if CFG["options"]["iedb_install_directory"] else "",
        # EXPERIMENTAL diagnostic knob, not a default-on fix -- see this module's CHANGELOG for the
        # full investigation. `--tumor-purity` only feeds pvacseq's own `vaf_clonal` estimate
        # (confirmed from pvactools/lib/aggregate_all_epitopes.py's calculate_clonal_vaf(): when
        # unset, pvacseq falls back to a crude, single-outlier-sensitive heuristic -- the single
        # largest observed DNA VAF under 0.6, capped at 0.5 -- rather than any real central-tendency
        # estimate of the sample's own VAF distribution), used only for the aggregated report's own
        # Tier/Subclonal *labels*. It is deliberately NOT used to derive --tdna-vaf below -- see that
        # param's own comment for why coupling them was wrong.
        purity_flag = f"--tumor-purity {CFG['options']['assumed_purity']}" if CFG["options"].get("assumed_purity") is not None else "",
        # Independent of purity/vaf_clonal on purpose (an earlier version of this rule derived this
        # from assumed_purity via pvacseq's own vaf_clonal/2 formula -- reverted, real, requested
        # correction): `--tdna-vaf` is coverage_filter's own hard coverage/filtering threshold on
        # whether a variant survives into filtered.tsv at all, completely independent of
        # `--tumor-purity`/vaf_clonal's own, separate role (an aggregate-report *tiering* label,
        # confirmed above). Raw tumour DNA VAF is not a reliable clonality estimator on its own --
        # it also depends on tumour purity, copy number, mutant multiplicity, LOH, and sampling
        # variance -- so using pvacseq's own clonality heuristic as a hard neoantigen-exclusion
        # criterion was never appropriate for our lymphoma/WGS use case, where subclonal mutations
        # should remain available for downstream neoantigen analysis rather than being silently
        # dropped from filtered.tsv. options.tdna_vaf (default 0.05, well below pvacseq's own stock
        # default of 0.25) is a fixed technical minimum -- low enough to retain credible
        # lower-VAF/subclonal candidates -- applied unconditionally, regardless of whether
        # assumed_purity is set. Tier labels (including "Subclonal") are never used anywhere in this
        # module to discard a candidate outright -- neither _pvactools_apply_additional_filter nor
        # _pvactools_join_neoantigen_maf filter on Tier at all, they only ever carry it through as
        # metadata (see this module's CHANGELOG).
        tdna_vaf_flag = f"--tdna-vaf {CFG['options'].get('tdna_vaf', 0.05)}"
    conda:
        CFG["conda_envs"]["pvactools"]
    container:
        None
    threads:
        CFG["threads"]["pvacseq_run"]
    resources:
        **CFG["resources"]["pvacseq_run"]
    # See _pvactools_vep_annotate's own comment for why this rule is grouped together with
    # _pvactools_subset_to_coding_variants/_pvactools_vep_annotate into a single cluster job, and
    # the resource-reservation tradeoff that comes with it.
    group:
        "vep_and_run"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 20
    shell:
        op.as_one_line("""
        (
        mkdir -p {output.outdir} &&
        pvacseq run {input.vep_vcf} {wildcards.tumour_id} "{params.alleles}" {params.algorithms}
        {output.outdir}
        -e1 {params.e1} -e2 {params.e2}
        -b {params.binding_threshold} -m {params.top_score_metric}
        --normal-sample-name {wildcards.normal_id}
        {params.iedb_flag}
        {params.purity_flag}
        {params.tdna_vaf_flag}
        -t {threads}
        ) > {log.stdout} 2>&1
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/').
#
# Filenames confirmed against a real completed pvacseq run (pvactools=7.1.3, 2026-09) -- an earlier
# draft guessed "{tumour_id}.filtered.tsv"/"{tumour_id}.all_epitopes.tsv" with no class infix; the
# real files are "{tumour_id}.MHC_I.filtered.tsv" etc. (and there's also a third, previously
# unaccounted-for "combined/" directory merging both classes -- see
# _pvactools_output_combined_filtered/_pvactools_output_combined_all_epitopes below). Still written
# defensively (symlink if the expected file exists, otherwise touch an empty file and warn) since a
# Class-I-only pair (no Class II alleles) never gets an MHC_Class_II/ or combined/ directory at all.
rule _pvactools_output_class1_filtered:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "mhc_class_i_filtered/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.tsv"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 30
    run:
        src = os.path.join(input.outdir, "MHC_Class_I", f"{wildcards.tumour_id}.MHC_I.filtered.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq Class I filtered output not found at {src} -- "
                  f"writing an empty placeholder (this pair likely has no Class I alleles typed).")
            open(output.tsv, "w").close()

rule _pvactools_output_class1_all_epitopes:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "mhc_class_i_all_epitopes/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.all_epitopes.tsv"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 30
    run:
        src = os.path.join(input.outdir, "MHC_Class_I", f"{wildcards.tumour_id}.MHC_I.all_epitopes.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq Class I all_epitopes output not found at {src} -- "
                  f"writing an empty placeholder (this pair likely has no Class I alleles typed).")
            open(output.tsv, "w").close()

# Class II outputs are only meaningful for a pair whose allele_summary.tsv actually contains a
# Class II allele -- for any other pair, pvacseq never writes an MHC_Class_II/ directory at all, so
# these just fall through to the defensive empty-placeholder-plus-warning path above.
rule _pvactools_output_class2_filtered:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "mhc_class_ii_filtered/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.tsv"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 30
    run:
        src = os.path.join(input.outdir, "MHC_Class_II", f"{wildcards.tumour_id}.MHC_II.filtered.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq Class II filtered output not found at {src} -- "
                  f"writing an empty placeholder (this pair likely has no Class II alleles typed).")
            open(output.tsv, "w").close()

rule _pvactools_output_class2_all_epitopes:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "mhc_class_ii_all_epitopes/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.all_epitopes.tsv"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 30
    run:
        src = os.path.join(input.outdir, "MHC_Class_II", f"{wildcards.tumour_id}.MHC_II.all_epitopes.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq Class II all_epitopes output not found at {src} -- "
                  f"writing an empty placeholder (this pair likely has no Class II alleles typed).")
            open(output.tsv, "w").close()

# The "combined/" directory (confirmed on a real run, not anticipated when this module was first
# designed) merges both classes into a single report -- likely the most useful single deliverable
# for most downstream analysis, since it saves having to concatenate MHC_Class_I/MHC_Class_II
# yourself. Same defensive pattern: only absent for a pair where neither class produced any alleles
# at all (shouldn't happen given _pvactools_prepare_allele_list's own both-empty assertion, but kept
# consistent with the other four output rules regardless).
rule _pvactools_output_combined_filtered:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "combined_filtered/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.tsv"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 30
    run:
        src = os.path.join(input.outdir, "combined", f"{wildcards.tumour_id}.Combined.filtered.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq combined filtered output not found at {src} -- "
                  f"writing an empty placeholder.")
            open(output.tsv, "w").close()

rule _pvactools_output_combined_all_epitopes:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "combined_all_epitopes/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.all_epitopes.tsv"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 30
    run:
        src = os.path.join(input.outdir, "combined", f"{wildcards.tumour_id}.Combined.all_epitopes.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq combined all_epitopes output not found at {src} -- "
                  f"writing an empty placeholder.")
            open(output.tsv, "w").close()

# Previously unexposed (noted as a known gap in this module's own CHANGELOG): pvacseq's own tiered,
# already-one-row-per-mutation aggregated report. Needed as the join source for
# _pvactools_join_neoantigen_maf below -- its own "ID" column encodes genomic coordinates as
# "{Chromosome}-{Start}-{Stop}-{Reference}-{Variant}" (confirmed against real production files;
# there are no separate Chromosome/Start/etc columns in this report, unlike the uncollapsed
# all_epitopes.tsv above). Same defensive symlink-or-placeholder pattern as the six rules above.
rule _pvactools_output_combined_aggregated:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "combined_aggregated/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.all_epitopes.aggregated.tsv"
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 30
    run:
        src = os.path.join(input.outdir, "combined", f"{wildcards.tumour_id}.Combined.all_epitopes.aggregated.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq combined aggregated output not found at {src} -- "
                  f"writing an empty placeholder.")
            open(output.tsv, "w").close()


# Additional standard filterings of the SAME already-computed combined all_epitopes.tsv, via
# `pvacseq binding_filter` -- pure post-processing on existing scores/percentiles (confirmed by
# reading pvactools' own lib/binding_filter.py directly: no model/IEDB calls, just re-applies
# thresholds to an existing report), so these cost essentially nothing beyond the primary pvacseq
# run itself -- no need to rerun anything to get a differently-thresholded view. See
# options.additional_filters for the preset definitions and the real observation that motivated
# this (pvacseq's own default "conservative" filtered.tsv left several samples with very few or
# zero surviving candidates). Localrule: this is fast, pure TSV post-processing, not a real compute
# job.
rule _pvactools_apply_additional_filter:
    input:
        all_epitopes = str(rules._pvactools_output_combined_all_epitopes.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "additional_filters/{filter_preset}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.tsv",
        # binding_filter alone does not collapse per mutation -- a single mutation with several
        # alleles/transcripts/peptide-windows surviving the loosened threshold would otherwise show
        # up as several rows here, inflating counts in a way that scales with how many HLA alleles
        # a sample has (a real, confounding bias for any cross-sample comparison). This uncollapsed
        # intermediate is filtered first, then collapsed to one row per mutation below -- kept as a
        # declared temp() output (not an ad hoc scratch file) so Snakemake tracks/cleans it.
        uncollapsed_tsv = temp(CFG["dirs"]["outputs"] + "additional_filters/{filter_preset}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.uncollapsed.tsv")
    log:
        stdout = CFG["logs"]["outputs"] + "additional_filters/{filter_preset}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/apply_additional_filter.log"
    params:
        preset = lambda wildcards: config["lcr-modules"]["pvactools"]["options"]["additional_filters"][wildcards.filter_preset],
        top_score_metric = CFG["options"]["top_score_metric"]
    wildcard_constraints:
        filter_preset = "|".join(CFG["options"]["additional_filters"].keys())
    conda:
        CFG["conda_envs"]["pvactools"]
    container:
        None
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority -- this rule is the last step, so it gets the highest value.
    priority: 40
    shell:
        op.as_one_line("""
        pvacseq binding_filter {input.all_epitopes} {output.uncollapsed_tsv}
        -b {params.preset[binding_threshold]}
        --binding-percentile-threshold {params.preset[binding_percentile_threshold]}
        --presentation-percentile-threshold {params.preset[presentation_percentile_threshold]}
        --percentile-threshold-strategy {params.preset[percentile_threshold_strategy]}
        -m {params.top_score_metric}
        > {log.stdout} 2>&1
        &&
        pvacseq top_score_filter {output.uncollapsed_tsv} {output.tsv}
        -b {params.preset[binding_threshold]}
        -m {params.top_score_metric}
        >> {log.stdout} 2>&1
        """)


# Stage 1 of the neoantigen MAF feature -- see this module's CHANGELOG for the full design
# reasoning, including the approaches tried and rejected before landing on this one. Joins
# pvactools' own aggregated report onto inputs.vcf2maf_maf -- the SAME native-build MAF this module
# already reads (_pvactools_extract_coding_regions already derives its own coding-region VCF subset
# from this exact file), so every pvactools-predicted mutation is guaranteed a corresponding row in
# it already. No second vcf2maf deployment needed: an earlier draft of this rule introduced one,
# left over from brainstorming Stage 2's own harmonized-build problem (a different problem -- see
# CHANGELOG) and never cleaned back out; removed as a pure simplification once noticed. Localrule:
# pure stdlib Python (csv/re), no conda env needed, matching _pvactools_prepare_allele_list's own
# convention.
#
# Output NOT permanent -- single consumer (_pvactools_crossmap_neoantigen_maf), cheap to
# regenerate, same temp()-chain-checkpoint reasoning already used elsewhere in this module.
rule _pvactools_join_neoantigen_maf:
    input:
        aggregated_report = str(rules._pvactools_output_combined_aggregated.output.tsv),
        native_maf = _pvactools_get_vcf2maf_maf,
        relaxed_filtered = lambda wildcards: _pvactools_get_filtered_for_preset(wildcards, "relaxed"),
        very_relaxed_filtered = lambda wildcards: _pvactools_get_filtered_for_preset(wildcards, "very_relaxed"),
        script = CFG["options"]["join_neoantigen_maf_script"]
    output:
        maf = temp(CFG["dirs"]["outputs"] + "neoantigen_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.native.maf"),
        audit = CFG["dirs"]["outputs"] + "neoantigen_maf_audit/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.neoantigen_maf_audit.tsv"
    log:
        stdout = CFG["logs"]["outputs"] + "neoantigen_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/join_neoantigen_maf.log"
    params:
        minimal_maf_columns = CFG["options"]["minimal_maf_columns"],
        neoantigen_columns = " ".join(shlex.quote(c) for c in CFG["options"]["neoantigen_maf_columns"]),
        # Column name -> threshold spec, as JSON (a nested dict doesn't fit this module's usual
        # space-joined-list param convention). One new pass/fail column per entry, named exactly by
        # its own config key -- see options.neoantigen_pass_filters' own comment in default.yaml for
        # the criteria supported and why NA auto-passes each one.
        pass_filters_json = shlex.quote(json.dumps(CFG["options"]["neoantigen_pass_filters"]))
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority.
    priority: 30
    shell:
        op.as_one_line("""
        python3 {input.script}
        --aggregated-report {input.aggregated_report}
        --native-maf {input.native_maf}
        --relaxed-filtered {input.relaxed_filtered}
        --very-relaxed-filtered {input.very_relaxed_filtered}
        --minimal-maf-columns {params.minimal_maf_columns}
        --neoantigen-columns {params.neoantigen_columns}
        --pass-filters-json {params.pass_filters_json}
        --output-maf {output.maf}
        --output-audit {output.audit}
        > {log.stdout} 2>&1
        """)


# Stage 2 of the neoantigen MAF feature: lifts Stage 1's native-build joined MAF to
# NEOANTIGEN_TARGET_BUILD, using ONLY the CrossMap coordinate-lift half of modules/vcf2maf/1.3's
# own pipeline (lcr-scripts' own convert_maf_coords.sh, same chain files/conda env that module's
# own _vcf2maf_crossmap rule uses) -- deliberately skipping vcf2maf's own reannotation half
# (_vcf2maf_reannotate/maf2maf.pl). Confirmed by reading convert_maf_coords.sh directly: it's
# column-agnostic, overwriting only Chromosome/Start_Position/End_Position in place and passing
# every other column through unchanged -- exactly what's needed, since Stage 1's own
# Hugo_Symbol/Variant_Classification/pvactools columns should ride through untouched rather than
# being re-derived (re-annotating at the target build is both unnecessary -- gene/protein identity
# doesn't depend on which build's coordinates it's expressed in -- and the specific step that would
# reintroduce the transcript-instability risk this module's CHANGELOG already documents rejecting).
# When this pair's native build already matches NEOANTIGEN_TARGET_BUILD, no chain file is needed at
# all -- Stage 1's own output is used as-is.
#
# This rule's own output already lands directly under 99-outputs/ (CFG["dirs"]["outputs"]) -- no
# separate "_pvactools_output_*" symlink wrapper needed, matching _pvactools_apply_additional_filter's
# own precedent (those wrapper rules exist specifically to extract files out of _pvactools_run's own
# nested, pvacseq-internal outdir; this rule has no such nesting to extract from).
#
# Real bug found while verifying this rule (not introduced by it): convert_maf_coords.sh's own
# multi-stage `awk | awk | cut | awk | sed | perl` reconstruction pipeline is intermittently flaky
# under its own `set -o pipefail` -- confirmed directly by running it 10 times against identical
# input, 1 of 10 runs failed with exit 141 (SIGPIPE) partway through that pipe chain, the other 9
# succeeded and lifted coordinates correctly. This is a real, pre-existing issue in that SHARED
# lcr-scripts file (also used by modules/vcf2maf/1.3's own _vcf2maf_crossmap) -- not something
# specific to this module, and not yet confirmed whether it reproduces on Linux (this was only
# observed on macOS) or is more of a pipe-scheduling artifact there. Rather than patch a shared
# script outside this module's own scope, this rule retries its own call up to 3 times (the
# `cp`-only branch never invokes the flaky script at all, so doesn't need this) -- a self-contained
# workaround, not a fix; worth reporting upstream if this turns out to matter on a real cluster too.
rule _pvactools_crossmap_neoantigen_maf:
    input:
        maf = str(rules._pvactools_join_neoantigen_maf.output.maf),
        convert_coord = CFG["options"]["convert_coord_script"]
    output:
        maf = CFG["dirs"]["outputs"] + "neoantigen_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.neoantigen.maf"
    log:
        stdout = CFG["logs"]["outputs"] + "neoantigen_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/crossmap_neoantigen_maf.log"
    params:
        chain = lambda wildcards: _pvactools_get_crossmap_chain(wildcards.genome_build) or ""
    conda:
        CFG["conda_envs"]["crossmap"]
    container:
        None
    threads:
        CFG["threads"]["crossmap_neoantigen_maf"]
    resources:
        **CFG["resources"]["crossmap_neoantigen_maf"]
    # See _pvactools_extract_coding_regions's own comment for why this whole chain has increasing
    # priority -- this rule is the last step for this feature, so it gets the highest value.
    priority: 40
    shell:
        op.as_one_line("""
        (
        if [[ -z "{params.chain}" ]];
        then cp {input.maf} {output.maf};
        else
        for attempt in 1 2 3;
        do {input.convert_coord} {input.maf} {params.chain} {output.maf} crossmap && break || true;
        sleep 1;
        done;
        test -s {output.maf};
        fi
        ) > {log.stdout} 2>&1
        """)


# Generates the target sentinels for each run, which generate the symlinks. Uses
# CFG["paired_runs"] (narrowed to pair_status == "matched" above) so tumour samples without a
# matched germline sample are never requested as targets -- HLA typing needs the patient's own
# germline sample. Class II/combined outputs are requested unconditionally for every pair -- whether
# a given pair actually has a Class II allele can't be known until _pvactools_prepare_allele_list has
# already run (a per-pair expand() can't depend on that rule's own output at DAG-build time), so a
# pair with no real Class II typing just gets the defensive empty-placeholder-plus-warning behaviour
# already built into the relevant output rules above, rather than being excluded from the target
# list entirely.
rule _pvactools_all:
    input:
        expand(
            [
                str(rules._pvactools_output_class1_filtered.output.tsv),
                str(rules._pvactools_output_class1_all_epitopes.output.tsv),
                str(rules._pvactools_output_class2_filtered.output.tsv),
                str(rules._pvactools_output_class2_all_epitopes.output.tsv),
                str(rules._pvactools_output_combined_filtered.output.tsv),
                str(rules._pvactools_output_combined_all_epitopes.output.tsv),
                str(rules._pvactools_crossmap_neoantigen_maf.output.maf)
            ],
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            tumour_id = CFG["paired_runs"]["tumour_sample_id"],
            normal_id = CFG["paired_runs"]["normal_sample_id"],
            pair_status = CFG["paired_runs"]["pair_status"]
        ),
        # Cross product (every pair x every preset), not zip -- every pair gets every preset.
        # allow_missing=True on the inner expand leaves {filter_preset} unresolved so the outer
        # expand can fill it in across all presets.
        expand(
            expand(
                [str(rules._pvactools_apply_additional_filter.output.tsv)],
                zip,
                seq_type = CFG["paired_runs"]["tumour_seq_type"],
                genome_build = CFG["paired_runs"]["tumour_genome_build"],
                tumour_id = CFG["paired_runs"]["tumour_sample_id"],
                normal_id = CFG["paired_runs"]["normal_sample_id"],
                pair_status = CFG["paired_runs"]["pair_status"],
                allow_missing = True
            ),
            filter_preset = list(CFG["options"]["additional_filters"].keys())
        )


##### CLEANUP #####


op.cleanup_module(CFG)
