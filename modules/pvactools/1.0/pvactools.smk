#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Griffith Lab (griffithlab/pVACtools)
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
# Real, confirmed gotcha (pvactools=7.1.3, 2026-09): pointing MHCFLURRY_DOWNLOADS_DIR at a fresh/
# custom directory disables mhcflurry's own default-release auto-detection ("No release defined"),
# so --release must be passed explicitly -- see options.mhcflurry_release's own comment.
rule _pvactools_download_mhcflurry_models:
    output:
        complete = touch(CFG["dirs"]["mhcflurry_reference"] + "mhcflurry_downloads.complete")
    log:
        stdout = CFG["logs"]["mhcflurry_reference"] + "download_mhcflurry_models.log"
    params:
        dir = CFG["dirs"]["mhcflurry_reference"],
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
        export MHCFLURRY_DOWNLOADS_DIR={params.dir} &&
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
    shell:
        op.as_one_line("""
        (
        bcftools view -R {input.regions} {input.vcf} |
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
        vep_plugins_dir = lambda wildcards, input: os.path.dirname(input.vep_frameshift)
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        None
    threads:
        CFG["threads"]["vep_annotate"]
    resources:
        **CFG["resources"]["vep_annotate"]
    shell:
        op.as_one_line("""
        (
        export PATH={params.vep_path}:${{PATH}} &&
        vep -i {input.vcf}
        --format vcf --vcf -o STDOUT
        --symbol --terms SO --mane_select --canonical --tsl --biotype --hgvs
        --fasta {input.fasta}
        --offline --cache --dir_cache {input.vep_cache}
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
        mhcflurry_dir = CFG["dirs"]["mhcflurry_reference"]
    conda:
        CFG["conda_envs"]["pvactools"]
    container:
        None
    threads:
        CFG["threads"]["pvacseq_run"]
    resources:
        **CFG["resources"]["pvacseq_run"]
    shell:
        op.as_one_line("""
        (
        export MHCFLURRY_DOWNLOADS_DIR={params.mhcflurry_dir} &&
        mkdir -p {output.outdir} &&
        pvacseq run {input.vep_vcf} {wildcards.tumour_id} "{params.alleles}" {params.algorithms}
        {output.outdir}
        -e1 {params.e1} -e2 {params.e2}
        -b {params.binding_threshold} -m {params.top_score_metric}
        --normal-sample-name {wildcards.normal_id}
        {params.iedb_flag}
        -t {threads}
        ) > {log.stdout} 2>&1
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/').
#
# Filenames below are PLACEHOLDERS based on general pVACtools documentation knowledge, NOT
# confirmed against a real pvacseq run's actual output directory in this session -- written
# defensively (symlink if the expected file exists, otherwise touch an empty file and warn) so a
# filename mismatch on the first real run doesn't hard-fail the whole pipeline. Correct the path
# here once confirmed via a real `ls -R` (see this module's CHANGELOG), rather than guessing twice.
rule _pvactools_output_class1_filtered:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "mhc_class_i_filtered/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.tsv"
    run:
        src = os.path.join(input.outdir, "MHC_Class_I", f"{wildcards.tumour_id}.filtered.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq Class I filtered output not found at {src} -- "
                  f"writing an empty placeholder. Confirm the real filename with `ls -R {input.outdir}` "
                  f"and fix this rule.")
            open(output.tsv, "w").close()

rule _pvactools_output_class1_all_epitopes:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "mhc_class_i_all_epitopes/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.all_epitopes.tsv"
    run:
        src = os.path.join(input.outdir, "MHC_Class_I", f"{wildcards.tumour_id}.all_epitopes.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq Class I all_epitopes output not found at {src} -- "
                  f"writing an empty placeholder. Confirm the real filename with `ls -R {input.outdir}` "
                  f"and fix this rule.")
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
    run:
        src = os.path.join(input.outdir, "MHC_Class_II", f"{wildcards.tumour_id}.filtered.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq Class II filtered output not found at {src} -- "
                  f"writing an empty placeholder. Confirm the real filename with `ls -R {input.outdir}` "
                  f"and fix this rule.")
            open(output.tsv, "w").close()

rule _pvactools_output_class2_all_epitopes:
    input:
        outdir = str(rules._pvactools_run.output.outdir),
        complete = str(rules._pvactools_run.output.complete)
    output:
        tsv = CFG["dirs"]["outputs"] + "mhc_class_ii_all_epitopes/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.all_epitopes.tsv"
    run:
        src = os.path.join(input.outdir, "MHC_Class_II", f"{wildcards.tumour_id}.all_epitopes.tsv")
        if os.path.isfile(src):
            op.relative_symlink(src, output.tsv, in_module = True)
        else:
            print(f"WARNING: expected pvacseq Class II all_epitopes output not found at {src} -- "
                  f"writing an empty placeholder. Confirm the real filename with `ls -R {input.outdir}` "
                  f"and fix this rule.")
            open(output.tsv, "w").close()


# Generates the target sentinels for each run, which generate the symlinks. Uses
# CFG["paired_runs"] (narrowed to pair_status == "matched" above) so tumour samples without a
# matched germline sample are never requested as targets -- HLA typing needs the patient's own
# germline sample. Class II outputs are requested unconditionally for every pair -- whether a given
# pair actually has a Class II allele can't be known until _pvactools_prepare_allele_list has
# already run (a per-pair expand() can't depend on that rule's own output at DAG-build time), so a
# pair with no real Class II typing just gets the defensive empty-placeholder-plus-warning behaviour
# already built into _pvactools_output_class2_filtered/_pvactools_output_class2_all_epitopes above,
# rather than being excluded from the target list entirely.
rule _pvactools_all:
    input:
        expand(
            [
                str(rules._pvactools_output_class1_filtered.output.tsv),
                str(rules._pvactools_output_class1_all_epitopes.output.tsv),
                str(rules._pvactools_output_class2_filtered.output.tsv),
                str(rules._pvactools_output_class2_all_epitopes.output.tsv)
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
