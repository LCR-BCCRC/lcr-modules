#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  sgillis
# Module Author:    sgillis
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["panel_of_normals"]`
CFG = op.setup_module(
    name = "panel_of_normals",
    version = "1.0",
    subdirectories = ["inputs", "cnvkit_targets_bed", "cnvkit_autobin", "cnvkit_coverage", "cnvkit_pon_cnn", "cnvkit_flat_ref_cnn",
                        "purecn_intervals", "purecn_mutect2", "purecn_merge_vcfs",
                        "purecn_coverage", "purecn_NormalDB", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _panel_of_normals_input_bam,
    _panel_of_normals_input_capspace,
    _panel_of_normals_canonical_capspace,
    _panel_of_normals_cnvkit_accessible_main_chrs,
    _panel_of_normals_cnvkit_output_beds,
    _panel_of_normals_cnvkit_output_cnn,
    _panel_of_normals_cnvkit_output_tsv,
    _panel_of_normals_cnvkit_output_flat_ref_beds,
    _panel_of_normals_cnvkit_output_flat_ref,
    _panel_of_normals_purecn_get_mappability,
    _panel_of_normals_purecn_gatk_interval_list,
    _panel_of_normals_purecn_gatk_interval_list_chrom,
    _panel_of_normals_purecn_gatk_interval_list_targets,
    _panel_of_normals_purecn_samples_map,
    _panel_of_normals_purecn_coverage_list,
    _panel_of_normals_purecn_output_intervals,
    _panel_of_normals_purecn_output_gatk_intervals,
    _panel_of_normals_purecn_output_gatk_intervals_targets,
    _panel_of_normals_purecn_output_mutect2_pon,
    _panel_of_normals_purecn_output_normal_panel_vcf,
    _panel_of_normals_purecn_output_database_cnvkit,
    _panel_of_normals_purecn_output_database_denovo,
    _panel_of_normals_output_samples_tsv,
    _panel_of_normals_all


# Get unique combos in the samples table for rule all
META = CFG["samples"][['seq_type', 'genome_build', 'capture_space']].drop_duplicates()

##### RULES #####

# ---------------------------------------------------------- #
# Shared files section (used by both CNVkit, PureCN, etc)
# ---------------------------------------------------------- #

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _panel_of_normals_input_bam:
    input:
        bam = ancient(CFG["inputs"]["sample_bam"])
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)

# Recreate index so the timestamp will always be later then the bam (for CNVkit)
rule _panel_of_normals_index_bam:
    input:
        bam = str(rules._panel_of_normals_input_bam.output.bam)
    output:
        bai = temp(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam.bai")
    log:
        log = CFG["logs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{sample_id}_index.log"
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["samtools"]
    shell:
        op.as_one_line("""
        samtools index -@ {threads} {input.bam} 2> {log.log} &&
        cd $(dirname {input.bam});
        if [[ -e {wildcards.sample_id}.bam.crai ]];
        then
            ln -s {wildcards.sample_id}.bam.crai {wildcards.sample_id}.bam.bai;
        fi
        """)

# Symlinks the capture bed file
rule _panel_of_normals_input_capspace:
    input:
        bed = ancient(reference_files("genomes/{genome_build}/capture_space/{capture_space}.padded.bed"))
    output:
        bed = CFG["dirs"]["inputs"] + "bed/{seq_type}--{genome_build}/{capture_space}.padded.bed"
    run:
        op.absolute_symlink(input.bed, output.bed)

# Removes non-canonical chroms from capture space bed
rule _panel_of_normals_canonical_capspace:
    input:
        bed = str(rules._panel_of_normals_input_capspace.output.bed)
    output:
        bed = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}.padded.canonical.bed"
    shell:
        op.as_one_line("""
        awk -F"\t" -v OFS="\t" '$1 !~ /(_|M|EBV|HIV)/' {input.bed} > {output.bed}
        """)

# Symlink fasta and it's index (for compatibility with $PURECN/IntervalFile.R)
rule _panel_of_normals_symlink_fasta:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        fasta = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa",
        fai = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa.fai"
    run:
        op.absolute_symlink(input.fasta, output.fasta)
        op.absolute_symlink(input.fai, output.fai)

checkpoint _panel_of_normals_input_chroms_withY:
    input:
        txt = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt")
    output:
        txt = CFG["dirs"]["inputs"] + "references/{genome_build}/main_chromosomes_withY.txt"
    run:
        op.absolute_symlink(input.txt, output.txt)


# Collects all normals per genome_build--capture_space combo
def _get_normal_bams_per_combo(wildcards):
    CFG = config["lcr-modules"]["panel_of_normals"]
    tbl = CFG["samples"]
    samples = tbl[(tbl.genome_build == wildcards.genome_build) & (tbl.capture_space == wildcards.capture_space)]["sample_id"].tolist()
    normals = expand(
        str(rules._panel_of_normals_input_bam.output.bam),
        sample_id = samples,
        allow_missing = True)
    return normals

# Collects all normal indexes per genome_build--capture_space combo
# need to do this separately so that they have a different name
def _get_indexes_per_combo(wildcards):
    CFG = config["lcr-modules"]["panel_of_normals"]
    tbl = CFG["samples"]
    samples = tbl[(tbl.genome_build == wildcards.genome_build) & (tbl.capture_space == wildcards.capture_space)]["sample_id"].tolist()
    normals = expand(
        str(rules._panel_of_normals_index_bam.output.bai),
        sample_id = samples,
        allow_missing = True)
    return normals

# Outputs the sample_ids per capture space combo
rule _panel_of_normals_record_samples:
    input:
        _get_normals_per_combo
    output:
        tsv = CFG["dirs"]["inputs"] + "sample_metadata_tsvs/{seq_type}--{genome_build}/{capture_space}_samples_metadata.tsv"
    params:
        metadata = CFG["samples"]
    threads: 1
    run:
        metadata_for_combo = params.metadata[(params.metadata.genome_build == wildcards.genome_build) & (params.metadata.capture_space == wildcards.capture_space)]
        metadata_for_combo.to_csv(output.tsv, sep="\t", index=False, na_rep='NA')

# ---------------------------------------------------------- #
# CNVkit section
# ---------------------------------------------------------- #

# Download gene annotation files
rule _panel_of_normals_cnvkit_get_refFlat:
    output:
        refFlat = CFG["dirs"]["inputs"] + "references/{genome_build}/refFlat.final.txt"
    params:
        url = "http://hgdownload.soe.ucsc.edu/goldenPath/",
        build = lambda w: "hg38" if "38" in str({w.genome_build}) else "hg19",
        txt = CFG["dirs"]["inputs"] + "references/{genome_build}/refFlat.txt",
        prefix = lambda w: "chr" if "hg" in str({w.genome_build}) else ""
    threads: 1
    shell:
        op.as_one_line("""
        wget {params.url}{params.build}/database/refFlat.txt.gz -O - | gzip -d
         > {params.txt} &&
         sed 's/chr/{params.prefix}/g' {params.txt} > {output.refFlat}
        """)

# divides regions into appropriately sized bins and adds gene annotations
rule _panel_of_normals_cnvkit_targets_bed:
    input:
        bed = str(rules._panel_of_normals_canonical_capspace.output.bed),
        refFlat = str(rules._panel_of_normals_get_refFlat.output.refFlat)
    output:
        targets = CFG["dirs"]["cnvkit_targets_bed"] + "{seq_type}--{genome_build}/{capture_space}/targets.bed"
    log:
        log = CFG["logs"]["cnvkit_targets_bed"] + "{seq_type}--{genome_build}/{capture_space}/annotate_targets.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but needs to be submitted
    resources:
        **CFG["resources"]["cnvkit"]["coverage"]
    shell:
        op.as_one_line("""
        cnvkit.py target {input.bed} --annotate {input.refFlat} --split -o {output.targets} &> {log.log}
        """)

rule _panel_of_normals_cnvkit_accessible_regions:
    input:
        fasta = str(rules._panel_of_normals_symlink_fasta.output.fasta)
    output:
        access = CFG["dirs"]["inputs"] + "references/access.{genome_build}.bed"
    log:
        log = CFG["logs"]["inputs"] + "{genome_build}_access.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["cnvkit"]["autobin"]
    shell:
        op.as_one_line("""
        cnvkit.py access {input.fasta} -o {output.access} &> {log.log}
        """)

# Filters to only the main chromosomes
rule _panel_of_normals_cnvkit_accessible_main_chrs:
    input:
        access = CFG["dirs"]["inputs"] + "references/access.{genome_build}.bed",
        main_txt = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt")
    output:
        access_main = CFG["dirs"]["inputs"] + "references/access_main.{genome_build}.bed"
    shell:
        op.as_one_line("""
        grep -F -f {input.main_txt} {input.access} > {output.access_main}
        """)


rule _panel_of_normals_cnvkit_autobin:
    input:
        access = str(rules._panel_of_normal_cnvkit_accessible_main_chrs.output.access_main),
        bam = _get_normals_per_combo,
        bai = _get_indexes_per_combo,
        targets = str(rules._panel_of_normals_cnvkit_targets_bed.output.targets),
        refFlat = str(rules._panel_of_normals_get_refFlat.output.refFlat)
    output:
        target = CFG["dirs"]["cnvkit_autobin"] + "{seq_type}--{genome_build}/{capture_space}/target_sites.target.bed",
        antitarget = CFG["dirs"]["cnvkit_autobin"] + "{seq_type}--{genome_build}/{capture_space}/target_sites.antitarget.bed"
    log:
        log = CFG["logs"]["cnvkit_autobin"] + "{seq_type}--{genome_build}/{capture_space}_autobin.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["cnvkit"]["autobin"]
    shell:
        op.as_one_line("""
        cnvkit.py autobin {input.bam} -t {input.targets} -g {input.access} --annotate
         {input.refFlat} --short-names --target-output-bed
         {output.target} --antitarget-output-bed {output.antitarget} &> {log.log}
        """)

# Get coverage for each sample
rule _panel_of_normals_cnvkit_coverage_target:
    input:
        bam = str(rules._panel_of_normals_input_bam.output.bam),
        bai = str(rules._panel_of_normals_index_bam.output.bai),
        target = str(rules._panel_of_normals_cnvkit_autobin.output.target),
    output:
        cov = CFG["dirs"]["cnvkit_coverage"] + "target/{seq_type}--{genome_build}/{capture_space}/{sample_id}.targetcoverage.cnn"
    log:
        log = CFG["logs"]["cnvkit_coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}_target.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["cnvkit"]["coverage"]
    resources:
        **CFG["resources"]["cnvkit"]["coverage"]
    shell:
        op.as_one_line("""
        cnvkit.py coverage {input.bam} {input.target} -o {output.cov} -p {threads}
         &> {log.log}
        """)

# Get coverage of anti-target sample
rule _panel_of_normals_cnvkit_coverage_antitarget:
    input:
        bam = str(rules._panel_of_normals_input_bam.output.bam),
        bai = str(rules._panel_of_normals_index_bam.output.bai),
        antitarget = str(rules._panel_of_normals_cnvkit_autobin.output.antitarget),
    output:
        cov = CFG["dirs"]["cnvkit_coverage"] + "antitarget/{seq_type}--{genome_build}/{capture_space}/{sample_id}.antitargetcoverage.cnn"
    log:
        log = CFG["logs"]["cnvkit_coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}_antitarget.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["cnvkit"]["coverage"]
    resources:
        **CFG["resources"]["cnvkit"]["coverage"]
    shell:
        op.as_one_line("""
        cnvkit.py coverage {input.bam} {input.antitarget} -o {output.cov} -p {threads}
         &> {log.log}
        """)

# Create reference coverage cnn file using all normal cnn files per combo
def _get_cnvkit_coverage_per_combo(wildcards):
    CFG = config["lcr-modules"]["panel_of_normals"]
    tbl = CFG["samples"]
    samples = tbl[(tbl.genome_build == wildcards.genome_build) & (tbl.capture_space == wildcards.capture_space)]["sample_id"].tolist()
    target_cov = expand(
        str(rules._panel_of_normals_cnvkit_coverage_target.output.cov),
        sample_id = samples,
        allow_missing = True)
    antitarget_cov = expand(
        str(rules._panel_of_normals_cnvkit_coverage_antitarget.output.cov),
        sample_id = samples,
        allow_missing = True)
    return{
        "control_target": target_cov,
        "control_antitarget": antitarget_cov
    }

rule _panel_of_normals_cnvkit_create_pon_ref:
    input:
        unpack(_get_cnvkit_coverage_per_combo),
        fasta = str(rules._panel_of_normals_symlink_fasta.output.fasta)
    output:
        pon = CFG["dirs"]["cnvkit_pon_cnn"] + "cnn/{seq_type}--{genome_build}/{capture_space}_normal_reference.cnn"
    log:
        log = CFG["logs"]["cnvkit_pon_cnn"] + "{seq_type}--{genome_build}/{capture_space}/pon_cnn.log"
    params:
        male_reference = CFG["options"]["cnvkit"]["male_ref"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["cnvkit"]["reference"]
    shell:
        op.as_one_line("""
        cnvkit.py reference {input.control_target} {input.control_antitarget}
         --fasta {input.fasta} -o {output.pon} {params.male_reference}
         &> {log.log}
        """)

#### The following rules create a "flat" reference file for edge cases where there are no normals and no equivalents
# Annotates target sites with refFlat file
rule _panel_of_normals_cnvkit_flat_ref_annotate_targets:
    input:
        bed = str(rules._panel_of_normals_canonical_capspace.output.bed),
        refFlat = str(rules._panel_of_normals_cnvkit_get_refFlat.output.refFlat)
    output:
        targets = CFG["dirs"]["cnvkit_flat_ref_cnn"] + "{seq_type}--{genome_build}/{capture_space}/target_sites.bed"
    log:
        log = CFG["logs"]["cnvkit_flat_ref_cnn"] + "{seq_type}--{genome_build}/{capture_space}/annotate_targets.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization
    resources:
        **CFG["resources"]["cnvkit"]["coverage"]
    shell:
        op.as_one_line("""
        cnvkit.py target {input.bed} --annotate {input.refFlat} --split -o {output.targets} &> {log.log}
        """)

# Create anti target regions bed
rule _panel_of_normals_cnvkit_flat_ref_antitargets:
    input:
        targets = str(rules._panel_of_normals_cnvkit_flat_ref_annotate_targets.output.targets),
        access_main = str(rules._panel_of_normals_cnvkit_filter_main_chrs.output.access_main)
    output:
        antitargets = CFG["dirs"]["cnvkit_flat_ref_cnn"] + "{seq_type}--{genome_build}/{capture_space}/antitarget_sites.bed"
    log:
        log = CFG["logs"]["cnvkit_flat_ref_cnn"] + "{seq_type}--{genome_build}/{capture_space}/antitargets.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization
    resources:
        **CFG["resources"]["cnvkit"]["coverage"]
    shell:
        op.as_one_line("""
        cnvkit.py antitarget {input.targets} -g {input.access_main} -o {output.antitargets} &> {log.log}
        """)

# Makes the flat reference file
rule _panel_of_normals_cnvkit_flat_ref:
    input:
        targets = str(rules._panel_of_normals_cnvkit_flat_ref_annotate_targets.output.targets),
        antitargets = str(rules._panel_of_normals_cnvkit_flat_ref_antitargets.output.antitargets),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        cnn = CFG["dirs"]["cnvkit_flat_ref_cnn"] + "{seq_type}--{genome_build}/{capture_space}/flat_reference.cnn"
    log:
        log = CFG["logs"]["cnvkit_flat_ref_cnn"] + "{seq_type}--{genome_build}/{capture_space}/flat_ref_cnn.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["cnvkit"]["reference"]
    resources:
        **CFG["resources"]["cnvkit"]["reference"]
    shell:
        op.as_one_line("""
        cnvkit.py reference -o {output.cnn} -f {input.fasta} -t {input.targets} -a {input.antitargets} &> {log.log}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _panel_of_normals_cnvkit_output_beds:
    input:
        target = str(rules._panel_of_normals_cnvkit_autobin.output.target),
        antitarget = str(rules._panel_of_normals_cnvkit_autobin.output.antitarget)
    output:
        target = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{capture_space}_target_sites.bed",
        antitarget = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{capture_space}_antitarget_sites.bed"
    wildcard_constraints: # needed in order not to clash with _panel_of_normals_flat_ref_output_beds
        seq_type="capture",
        genome_build='|'.join(META["genome_build"]),
        capture_space='|'.join(META["capture_space"])
    run:
        op.relative_symlink(input.target, output.target, in_module= True)
        op.relative_symlink(input.antitarget, output.antitarget, in_module= True)

rule _panel_of_normals_cnvkit_output_cnn:
    input:
        cnn = str(rules._panel_of_normals_cnvkit_create_pon_ref.output.pon)
    output:
        cnn = CFG["dirs"]["outputs"] + "cnn/{seq_type}--{genome_build}/{capture_space}_normal_reference.cnn"
    wildcard_constraints: # needed in order not to clash with _panel_of_normals_output_flat_ref
        seq_type="capture",
        genome_build='|'.join(META["genome_build"]),
        capture_space='|'.join(META["capture_space"])
    run:
        op.relative_symlink(input.cnn, output.cnn, in_module= True)

rule _panel_of_normals_cnvkit_output_flat_ref_beds:
    input:
        target = str(rules._panel_of_normals_cnvkit_flat_ref_annotate_targets.output.targets),
        antitarget = str(rules._panel_of_normals_cnvkit_flat_ref_antitargets.output.antitargets)
    output:
        target = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{capture_space}_target_sites.bed",
        antitarget = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{capture_space}_antitarget_sites.bed"
    wildcard_constraints: # needed in order not to clash with _panel_of_normals_output_beds
        seq_type="capture",
        genome_build='|'.join(CFG["options"]["cnvkit"]["flat_ref_combos"]["genome_builds"]),
        capture_space='|'.join(CFG["options"]["cnvkit"]["flat_ref_combos"]["capture_spaces"])
    run:
        op.relative_symlink(input.target, output.target, in_module= True)
        op.relative_symlink(input.antitarget, output.antitarget, in_module= True)

rule _panel_of_normals_cnvkit_output_flat_ref:
    input:
        cnn = str(rules._panel_of_normals_cnvkit_flat_ref.output.cnn)
    output:
        cnn = CFG["dirs"]["outputs"] + "cnn/{seq_type}--{genome_build}/{capture_space}_normal_reference.cnn"
    wildcard_constraints: # needed in order not to clash with _panel_of_normals_output_cnn
        seq_type="capture",
        genome_build='|'.join(CFG["options"]["cnvkit"]["flat_ref_combos"]["genome_builds"]),
        capture_space='|'.join(CFG["options"]["cnvkit"]["flat_ref_combos"]["capture_spaces"])
    run:
        op.relative_symlink(input.cnn, output.cnn, in_module= True)


# ---------------------------------------------------------- #
# PureCN section
# ---------------------------------------------------------- #
# Download PureCN-recommended mappability files
rule _panel_of_normals_purecn_get_mappability:
    output:
        bw = CFG["dirs"]["inputs"] + "references/{genome_build}/mappability.bw"
    params:
        url = lambda w: {
            "hg19": "https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig",
            "hg38": "https://s3.amazonaws.com/purecn/GCA_000001405.15_GRCh38_no_alt_analysis_set_100.bw"
        }[config["lcr-modules"]["purecn"]["options"]["genome_builds_map"].get(w.genome_build, w.genome_build)]
    conda:
        CFG["conda_envs"]["wget"]
    shell:
        op.as_one_line("""
            wget {params.url} -O {output.bw}
        """)
# TODO: remove or re-implement after testing
# PureCN script to create intervals from capture space bed
rule _panel_of_normals_purecn_setinterval:
    input:
        genome = str(rules._panel_of_normals_symlink_fasta.output.fasta),
        bed = str(rules._panel_of_normals_canonical_capspace.output.bed),
        bw = str(rules._panel_of_normals_purecn_get_mappability.output.bw)
    output:
        intervals = CFG["dirs"]["purecn_intervals"] + "{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_intervals.txt"
    params:
        genome_build = lambda w: CFG["options"]["purecn"]["genome_builds_map"][w.genome_build],
        # intervalfile_script = CFG["software"]["intervalfile_script"],
        force = CFG["options"]["purecn"]["setinterval"]["force"],
        opts = CFG["options"]["purecn"]["setinterval"]["opts"]
    conda:
        CFG["conda_envs"]["purecn"]
    resources:
        **CFG["resources"]["purecn"]["setinterval"]
    threads:
        CFG["threads"]["purecn"]["setinterval"]
    log:
        CFG["logs"]["purecn_intervals"] + "{seq_type}--{genome_build}/{capture_space}/purecn_setinterval.log"
    shell:
        op.as_one_line("""
            PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
            Rscript --vanilla $PURECN/IntervalFile.R --in-file {input.bed}
            --fasta {input.genome} --out-file {output.intervals}
            --genome {params.genome_build}
            --mappability {input.bw} {params.force} {params.opts} > {log} 2>&1
        """)

#### Calculating coverage with GATK4 since PureCN cannot take crams directly

# Remove header lines and columns, so only regions remain in format chr:##-##
rule _panel_of_normals_purecn_gatk_interval_list:
    input:
        intervals = str(rules._panel_of_normals_purecn_setinterval.output.intervals)
    output:
        gatk_intervals = CFG["dirs"]["purecn_intervals"] + "{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_intervals_gatk.list"
    shell:
        op.as_one_line("""
            egrep -i  '^.*:.*-.*' {input.intervals} | awk '{{print $1}}' > {output.gatk_intervals}
        """)

# Split by chrom, used later in GATK depthOfCoverage
rule _panel_of_normals_purecn_gatk_interval_list_chrom:
    input:
        gatk_intervals = str(rules._panel_of_normals_purecn_gatk_interval_list.output.gatk_intervals)
    output:
        chrom_int = CFG["dirs"]["purecn_intervals"] + "{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_{chrom}.intervals_gatk.list"
    log:
        CFG["logs"]["purecn_intervals"] + "{seq_type}--{genome_build}/{capture_space}/purecn_gatk_intervals_{chrom}.log"
    shell:
        op.as_one_line(
        """
            num_intervals=$( {{ egrep -i '^{wildcards.chrom}:.*-.*' {input.gatk_intervals} || true; }} | wc -l );
            if [[ $num_intervals -eq 0 ]]; then
                echo "No intervals found for chromosome {wildcards.chrom} in {input.gatk_intervals}" | tee {log};
                echo "{wildcards.chrom}" > {output.chrom_int};
            else
                echo "Found $num_intervals intervals for chromosome {wildcards.chrom} in {input.gatk_intervals}" | tee {log};
                egrep -i '^{wildcards.chrom}:.*-.*' {input.gatk_intervals} > {output.chrom_int};
            fi
        """
        )

# Drop off-target regions for MuTect2
rule _panel_of_normals_purecn_gatk_interval_list_targets:
    input:
        intervals = str(rules._panel_of_normals_purecn_setinterval.output.intervals)
    output:
        gatk_targets = CFG["dirs"]["purecn_intervals"] + "{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_intervals_gatk_targets.list"
    shell:
        op.as_one_line("""
            egrep -i  '^.*:.*-.*' {input.intervals} | egrep "TRUE" | awk '{{print $1}}' > {output.gatk_targets}
        """)

# Run Mutect2 to get germline variants
rule _panel_of_normals_purecn_mutect2_germline:
    input:
        bam = str(rules._panel_of_normals_input_bam.output.bam),
        dbsnp = ancient(reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")),
        fasta = str(rules._panel_of_normals_symlink_fasta.output.fasta),
        gnomad = ancient(reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")),
        target_regions = str(rules._panel_of_normals_purecn_gatk_interval_list_targets.output.gatk_targets)
    output:
        vcf = temp(CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.vcf.gz.tbi"),
        stats = temp(CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.vcf.gz.stats"),
        f1r2 = temp(CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.f1r2.tar.gz")
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["purecn"]["mutect"]["mutect2_opts"]
    log:
        CFG["logs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/mutect2_{chrom}.log"
    conda:
        CFG["conda_envs"]["mutect"]
    group:
        "mutect2_per_chrom"
    shell:
        op.as_one_line(
        """
            if [[ $(egrep "^{wildcards.chrom}:" {input.target_regions} | wc -l) -eq 0 ]]; then
                echo "No intervals found for chromosome {wildcards.chrom} in {input.target_regions}" | tee {log};
                gatk Mutect2
                --java-options "-Xmx{params.mem_mb}m" {params.opts}
                --genotype-germline-sites true
                --genotype-pon-sites true
                --max-mnp-distance 0
                --germline-resource {input.gnomad}
                -R {input.fasta}
                -L {wildcards.chrom}:1-100
                -I {input.bam}
                -O {output.vcf}
                --f1r2-tar-gz {output.f1r2}
                > {log} 2>&1;
            else
                echo "Found intervals for chromosome {wildcards.chrom} in {input.target_regions}" | tee {log};
                gatk Mutect2
                --java-options "-Xmx{params.mem_mb}m" {params.opts}
                --genotype-germline-sites true
                --genotype-pon-sites true
                --max-mnp-distance 0
                --germline-resource {input.gnomad}
                -R {input.fasta}
                -L {wildcards.chrom}
                -L {input.target_regions}
                -isr INTERSECTION
                -I {input.bam}
                -O {output.vcf}
                --f1r2-tar-gz {output.f1r2}
                > {log} 2>&1;
            fi
        """
        )

# Concat vcfs and indexes across chroms per sample
def _get_mutect2_chr_vcfs(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._panel_of_normals_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    vcfs = expand(
        CFG["dirs"]["purecn_mutect2"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}/{{sample_id}}.{chrom}.vcf.gz",
        chrom = chrs
    )
    return(vcfs)


def _get_mutect2_chr_tbis(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._panel_of_normals_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    tbis = expand(
        CFG["dirs"]["purecn_mutect2"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}/{{sample_id}}.{chrom}.vcf.gz.tbi",
        chrom = chrs
    )
    return(tbis)


rule _panel_of_normals_purecn_concat_vcf_per_sample:
    input:
        vcf = _get_mutect2_chr_vcfs,
        tbi = _get_mutect2_chr_tbis
    output:
        vcf = temp(CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.vcf.gz"),
        tbi = temp(CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.vcf.gz.tbi")
    log:
        stderr = CFG["logs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.concat_vcf.stderr.log"
    resources:
        **CFG["resources"]["purecn"]["mutect"]
    threads:
        CFG["threads"]["purecn"]["mutect"]
    conda:
        CFG["conda_envs"]["bcftools"]
    group:
        "mutect2_per_chrom"
    shell:
        op.as_one_line("""
        bcftools concat {input.vcf} -Oz -o {output.vcf} &> {log.stderr};
        tabix -p vcf {output.vcf} &>> {log.stderr}
        """)

def _get_mutect2_chr_stats(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._panel_of_normals_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    stats = expand(
        CFG["dirs"]["purecn_mutect2"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}/{{sample_id}}.{chrom}.vcf.gz.stats",
        chrom = chrs
    )
    return(stats)


# Merge mutect2 stats per chrom, for FilterMutectCalls rule
rule _panel_of_normals_purecn_merge_stats_per_sample:
    input:
        stats = _get_mutect2_chr_stats
    output:
        stats = CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.vcf.gz.stats"
    log:
        CFG["logs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/merge_stats.log"
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["purecn"]["mutect"]
    threads:
        CFG["theads"]["purecn"]["mutect"]
    group:
        "mutect2_per_chrom"
    shell:
        op.as_one_line("""
        gatk MergeMutectStats $(for i in {input.stats}; do echo -n "-stats $i "; done)
        -O {output.stats} > {log} 2>&1
        """)

# Get pileup summaries, for GATK CalculateContamination
rule _panel_of_normals_purecn_pileup_summaries:
    input:
        bam = str(rules._panel_of_normals_input_bam.output.bam),
        snps = ancient(reference_files("genomes/{genome_build}/gatk/mutect2_small_exac.{genome_build}.vcf.gz")),
        fasta = str(rules._panel_of_normals_symlink_fasta.output.fasta)
    output:
        pileup = CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/pileupSummary.table"
    log:
        CFG["logs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/pileupSummary.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["purecn"]["mutect"]
    threads:
        CFG["threads"]["purecn"]["mutect"]
    group:
        "mutect2_per_chrom"
    shell:
        op.as_one_line("""
        gatk GetPileupSummaries
            --java-options "-Xmx{params.mem_mb}m"
            -I {input.normal_bam}
            -R {input.fasta}
            -V {input.snps}
            -L {input.snps}
            -O {output.pileup}
            > {log} 2>&1
        """)

# Calculate contamination
rule _panel_of_normals_purecn_calc_contamination:
    input:
        pileup = str(rules._panel_of_normals_purecn_pileup_summaries.output.pileup)
    output:
        segments =  CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/segments.table",
        contamination =  CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/contamination.table"
    log:
        CFG["logs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/calc_contam.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["purecn"]["mutect"]
    threads:
        CFG["threads"]["purecn"]["mutect"]
    shell:
        op.as_one_line("""
        gatk CalculateContamination
            --java-options "-Xmx{params.mem_mb}m"
            -I {input.pileup}
            -tumor-segmentation {output.segments}
            -O {output.contamination}
            > {log} 2>&1
        """)

# Learn read orientation model
def _get_mutect2_chr_f1r2(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._panel_of_normals_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    f1r2 = expand(
        CFG["dirs"]["purecn_mutect2"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}/{{sample_id}}.{chrom}.f1r2.tar.gz",
        chrom = chrs
    )
    return(f1r2)

rule _panel_of_normals_purecn_calc_contamination:
    input:
        f1r2 = _get_mutect2_chr_f1r2
    output:
        model =  CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/read-orientation-model.tar.gz"
    log:
        CFG["logs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/learn_orient_model.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["purecn"]["mutect"]
    threads:
        CFG["threads"]["purecn"]["mutect"]
    shell:
        op.as_one_line("""
        inputs=$(for input in {input.f1r2}; do printf -- "-I $input "; done);
        gatk LearnReadOrientationModel
        --java-options "-Xmx{params.mem_mb}m"
        $inputs -O {output.model}
        > {log} 2>&1
        """)

# Marks variants filtered or PASS annotations
rule _panel_of_normals_purecn_annotate_vcf:
    input:
        vcf = str(rules._panel_of_normals_purecn_concat_vcf_per_sample.output.vcf),
        tbi = str(rules._panel_of_normals_purecn_concat_vcf_per_sample.output.tbi),
        stat = str(rules._panel_of_normals_purecn_merge_stats_per_sample.output.stats),
        segments = str(rules._panel_of_normals_purecn_calc_contamination.output.segments),
        contamination = str(rules._panel_of_normals_purecn_calc_contamination.output.contamination),
        model = str(rules._panel_of_normals_purecn_calc_contamination.output.model),
        fasta = str(rules._panel_of_normals_symlink_fasta.output.fasta)
    output:
        vcf = temp(CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/annotated.vcf.gz")
    log:
        CFG["logs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/annotate_vcf.log",
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["purecn"]["mutect"]["annotate"]
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["purecn"]["mutect"]
    threads: 1
    shell:
        op.as_one_line("""
        gatk FilterMutectCalls --java-options "-Xmx{params.mem_mb}m"
            {params.opts}
            -V {input.vcf}
            -R {input.fasta}
            --tumor-segmentation {input.segments}
            --contamination-table {input.contamination}
            --ob-priors {input.model}
            -O {output.vcf}
            > {log} 2>&1
        """)

# Filters for PASS and germline variants
# This will only take somatic ones if filtering just for PASSED (need to still maintain germline ones)
rule _panel_of_normals_purecn_mutect2_filter_vcf:
    input:
        vcf = str(rules._panel_of_normals_purecn_annotate_vcf.output.vcf)
    output:
        vcf = CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}_passed.vcf.gz",
        tbi = CFG["dirs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}_passed.vcf.gz.tbi"
    params:
        filter_for_opts = CFG["options"]["purecn"]["mutect"]["filter_for"],
        filter_out_opts = CFG["options"]["purecn"]["mutect"]["filter_out"]
    log:
        stderr = CFG["logs"]["purecn_mutect2"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{sample_id}/filter_vcf.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        **CFG["resources"]["purecn"]["mutect"]
    threads: 1
    shell:
        op.as_one_line("""
        bcftools view {params.filter_for_opts} -e "{params.filter_out_opts}" {input.vcf} |
            bcftools norm -m - -Oz -o {output.vcf} 2> {log.stderr};
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

###### Creating a panel of normals vcf for MuTect2 variant calling in tumours

def _get_normal_vcfs_per_combo(wildcards):
    CFG = config["lcr-modules"]["panel_of_normals"]
    tbl = CFG["samples"]
    samples = tbl[(tbl.genome_build == wildcards.genome_build) & (tbl.capture_space == wildcards.capture_space)]["sample_id"].tolist()
    normals = expand(
        str(rules._panel_of_normals_purecn_mutect2_filter_vcf.output.vcf),
        sample_id = samples,
        allow_missing = True)
    return normals

def _get_normal_tbis_per_combo(wildcards):
    CFG = config["lcr-modules"]["panel_of_normals"]
    tbl = CFG["samples"]
    samples = tbl[(tbl.genome_build == wildcards.genome_build) & (tbl.capture_space == wildcards.capture_space)]["sample_id"].tolist()
    normals = expand(
        str(rules._panel_of_normals_purecn_mutect2_filter_vcf.output.tbi),
        sample_id = samples,
        allow_missing = True)
    return normals

# Create sample map file
rule _panel_of_normals_purecn_samples_map:
    input:
        normal = _get_normal_vcfs_per_combo
    output:
        samples_map = CFG["dirs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/samples_map.txt",
        done = CFG["dirs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/samples_map.done"
    shell:
        """
        for samples in {input.normal}
        do
            name=$(basename $samples)
            name=${{name/.vcf.gz}}
            echo -e "$name\t$samples" >> {output.map_sample}
        done &&
        touch {output.done}
        """

# Create genomicsDb used later in creating mutect2 pon vcf
rule _panel_of_normals_purecn_gatk_genomicsDbimport:
    input:
        normal = _get_normal_vcfs_per_combo,
        normal_tbi = _get_normal_tbis_per_combo,
        map_sample = str(rules._panel_of_normals_purecn_samples_map.output.samples_map),
        done = str(rules._panel_of_normals_purecn_samples_map.output.done),
        target_regions = str(rules._panel_of_normals_canonical_capspace.output.bed)
    output:
        touch(CFG["dirs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/genomicsdb.done")
    log:
        CFG["logs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/genomicsdb.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        db_path = CFG["dirs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/genomicsdb/"
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["purecn"]["mutect"]
    threads: 1
    shell:
        op.as_one_line("""
        gatk GenomicsDBImport --java-options "-Xmx{params.mem_mb}m"
        -L {input.target_regions}
        --sample-name-map {input.map_sample}
        --genomicsdb-workspace-path {params.db_path} --lenient --merge-input-intervals TRUE
        --overwrite-existing-genomicsdb-workspace TRUE
        > {log} 2>&1
        """)

# Create mutect2 pon vcf
rule _panel_of_normals_purecn_mutect2_pon:
    input:
        done = str(rules._panel_of_normals_purecn_gatk_genomicsDbimport.output),
        fasta = str(rules._panel_of_normals_symlink_fasta.output.fasta)
    output:
        pon = CFG["dirs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/{capture_space}_mutect2_pon.vcf.gz"
    log:
        CFG["logs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/mutect2_pon_vcf.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = "gendb://" + CFG["dirs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/genomicsdb/"
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["purecn"]["mutect"]
    threads: 1
    shell:
        op.as_one_line("""
        gatk CreateSomaticPanelOfNormals --java-options "-Xmx{params.mem_mb}m"
        --reference {input.fasta}
        --variant {params.opts}
        -O {output.pon}
        > {log} 2>&1
        """)

###### Create panel of normals database files for PureCN

# Creating a panel of normals vcf for PureCN
rule _panel_of_normals_purecn_merge_vcfs:
    input:
        normal = _get_normal_vcfs_per_combo,
        normal_tbi = _get_normal_tbis_per_combo,
    output:
        normal_panel = CFG["dirs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/{capture_space}_normalpanel.vcf.gz",
        normal_panel_tbi = CFG["dirs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/{capture_space}_normalpanel.vcf.gz.tbi"
    log:
        dirs["logs"]["purecn_merge_vcfs"] + "{seq_type}--{genome_build}/{capture_space}/normalpanel.vcf.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads: 1
    shell:
        op.as_one_line("""
            bcftools merge {input.normal} -Oz -o {output.normal_panel} --force-samples &> {log};
            tabix -p vcf {output.normal_panel} &>> {log}
        """)

# PureCN - by extension rsamtools, does not have CRAM compatibility, even with R v4.5.3
# https://github.com/Bioconductor/Rsamtools/issues/21
# https://github.com/Bioconductor/Rsamtools/issues/56
# Work around is to use GATK to also calculate coverage and then feed it into pureCN for GC normalization

rule _panel_of_normals_purecn_gatk_depthOfCoverage:
    input:
        bam = str(rules._panel_of_normals_input_bam.output.bam),
        bai = str(rules._panel_of_normals_index_bam.output.bai),
        intervals =  str(rules._panel_of_normals_purecn_gatk_interval_list_chrom.output.chrom_int),
        fasta = str(rules._panel_of_normals_symlink_fasta.output.fasta)
    output:
        coverage = CFG["dirs"]["purecn_coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.sample_interval_summary",
        statistics = temp(CFG["dirs"]["purecn_coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.sample_interval_statistics")
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["purecn"]["coverage"]["depth_coverage"],
        base_name = CFG["dirs"]["purecn_coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}"
    log:
        CFG["logs"]["purecn_coverage"] + "{seq_type}--{genome_build}/{capture_space}/gatk_coverage/{sample_id}.{chrom}.log"
    conda:
        CFG["conda_envs"]["mutect"]
    group:
        "purecn_coverage"
    shell:
        op.as_one_line("""
        gatk DepthOfCoverage
            --java-options "-Xmx{params.mem_mb}m"
            {params.opts}
            --omit-depth-output-at-each-base
            --omit-locus-table
            --omit-per-sample-statistics
            --interval-merging-rule OVERLAPPING_ONLY
            -R {params.genome_fasta}
            -I {input.bam}
            -O {params.base_name}
            -L {input.intervals}
            > {log} 2>&1
        """)


def _get_gatk_chr_cov_depth(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._panel_of_normals_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    coverage = expand(
        CFG["dirs"]["purecn_coverage"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}/{{sample_id}}.{chrom}.sample_interval_summary",
        chrom = chrs
    )
    return(coverage)

def _get_gatk_chr_cov_stats(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._panel_of_normals_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    statistics = expand(
        CFG["dirs"]["purecn_coverage"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}/{{sample_id}}.{chrom}.sample_interval_statistics",
        chrom = chrs
    )
    return(statistics)

rule _panel_of_normals_purecn_gatk_coverage_concatenate_depths:
    input:
        depth = _get_gatk_chr_cov_depth,
        statistics = _get_gatk_chr_cov_stats,
    output:
        depth = CFG["dirs"]["purecn_coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.sample_interval_summary.gz"
    resources:
        **CFG["resources"]["purecn"]["coverage"]
    threads:
        CFG["threads"]["coverage"]
    group:
        "purecn_coverage"
    shell:
        op.as_one_line("""
        file1=$(echo {input.depth} | cut -d " " -f1 ) ;
        head -n 1 $file1 | gzip > {output.depth} ;
        for sample in {input.depth};
        do
            awk '(NR > 1)' $sample | gzip >> {output.depth} ;
        done
        """)


# pureCN Coverage.R used to normalize by GC, since depth is calculated above by GATK
# These pureCN coverage files are used in making the PureCN DB that is used in the denovo case
rule _panel_of_normals_purecn_coverage:
    input:
        bam = str(rules._panel_of_normals_input_bam.output.bam),
        bai = str(rules._panel_of_normals_index_bam.output.bai),
        intervals =  str(rules._panel_of_normals_purecn_setinterval.output.intervals),
        coverage = str(rules._panel_of_normals_purecn_gatk_coverage_concatenate_depths.output.depth)
    output:
        coverage = CFG["dirs"]["purecn_coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}_coverage_loess.txt.gz"
    params:
        name = "{sample_id}",
        coverage_script = CFG["software"]["purecn"]["coverage_script"],
        outdir = CFG["dirs"]["purecn_coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}",
        force =  CFG["options"]["purecn"]["coverage"]["force"],
        opt =  CFG["options"]["purecn"]["coverage"]["opts"],
        genome_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    conda:
        CFG["conda_envs"]["purecn"]
    resources:
        **CFG["resources"]["purecn"]["coverage"]
    threads:
        CFG["threads"]["purecn"]["coverage"]
    log:
        CFG["logs"]["purecn_coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}_coverage_loess.log"
    shell:
        op.as_one_line("""
            PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
            export R_LIBS=$CONDA_DEFAULT_ENV/lib/R/library/ ;
            echo -e "Using {params.coverage_script} instead of default $PURECN/Coverage.R..." ;
            Rscript --vanilla {params.coverage_script}  --out-dir {params.outdir}
            --bam {input.bam}
            --name {params.name}
            --reference {params.genome_fasta}
            --coverage {input.coverage}
            --intervals {input.intervals} {params.force} {params.opt} > {log} 2>&1
        """)

def _get_normals_coverage_per_combo(wildcards):
    CFG = config["lcr-modules"]["panel_of_normals"]
    tbl = CFG["samples"]
    samples = tbl[(tbl.genome_build == wildcards.genome_build) & (tbl.capture_space == wildcards.capture_space)]["sample_id"].tolist()
    normals = expand(
        str(rules._panel_of_normals_purecn_coverage.output.coverage),
        sample_id = samples,
        allow_missing = True)
    return normals

# Write out a list of coverage files
rule _panel_of_normals_purecn_coverage_list:
    input:
        coverage = _get_normals_coverage_per_combo
    output:
        cov_list = CFG["dirs"]["panel_of_normals"] + "{seq_type}--{genome_build}/{capture_space}/cov_list.txt"
    shell:
        op.as_one_line("""
            ls -a {input.coverage} | cat > {output.cov_list}
        """)

# Setting a mapping bias database - used for cnvkit segs
rule _panel_of_normals_purecn_database_cnvkit:
    input:
        normal_panel = str(rules._panel_of_normals_purecn_merge_vcfs.output.normal_panel),
        normal_panel_tbi = str(rules._panel_of_normals_purecn_merge_vcfs.output.normal_panel_tbi)
    output:
        mapping_bias = CFG["dirs"]["purecn_NormalDB"] + "{seq_type}--{genome_build}/{capture_space}/purecn_cnvkit_normal/mapping_bias_{genome_build}_{capture_space}.rds"
    log:
        CFG["logs"]["purecn_NormalDB"] + "{seq_type}--{genome_build}/{capture_space}/purecn_cnvkit_normaldb.log"
    params:
        dirOut = CFG["dirs"]["purecn_NormalDB"] + "{seq_type}--{genome_build}/{capture_space}/purecn_cnvkit_normal/",
        genome = "{genome_build}",
        capture_space = "{capture_space}"
    conda:
        CFG["conda_envs"]["purecn"]
    resources:
        **CFG["resources"]["purecn"]["normalDB"]
    threads:
        CFG["threads"]["purecn"]["normalDB"]
    shell:
        op.as_one_line("""
            echo $CONDA_DEFAULT_ENV ;
            PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
            export R_LIBS=$CONDA_DEFAULT_ENV/lib/R/library/ ;
            Rscript --vanilla $PURECN/NormalDB.R --out-dir {params.dirOut} --normal-panel {input.normal_panel}
            --assay {params.capture_space} --genome {params.genome} --force > {log} 2>&1
        """)

# For pureCN de novo PSCBS seg method, using its own coverage files
    rule _panel_of_normals_purecn_database_denovo:
        input:
            normal_panel = str(rules._panel_of_normals_purecn_merge_vcfs.output.normal_panel),
            normal_panel_tbi = str(rules._panel_of_normals_purecn_merge_vcfs.output.normal_panel_tbi),
            cov_list = str(rules._panel_of_normals_purecn_coverage_list.output.cov_list)
        output:
            mapping_bias = CFG["dirs"]["purecn_NormalDB"] + "{seq_type}--{genome_build}/{capture_space}/purecn_denovo_normal/mapping_bias_{genome_build}_{capture_space}.rds",
            normal_db = CFG["dirs"]["purecn_NormalDB"] + "{seq_type}--{genome_build}/{capture_space}/purecn_denovo_normal/normalDB_{genome_build}_{capture_space}.rds"
        log:
            CFG["logs"]["purecn_NormalDB"] + "{seq_type}--{genome_build}/{capture_space}/purecn_denovo_normaldb.log"
        params:
            dirOut = CFG["dirs"]["purecn_NormalDB"] + "{seq_type}--{genome_build}/{capture_space}/purecn_denovo_normal/",
            genome = "{genome_build}",
            platform = "{capture_space}"
        conda:
            CFG["conda_envs"]["purecn"]
        resources:
            **CFG["resources"]["purecn"]["normalDB"]
        threads:
            CFG["threads"]["purecn"]["normalDB"]
        shell:
            op.as_one_line("""
                echo $CONDA_DEFAULT_ENV ;
                PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
                export R_LIBS=$CONDA_DEFAULT_ENV/lib/R/library/ ;
                Rscript --vanilla $PURECN/NormalDB.R --out-dir {params.dirOut} --normal-panel {input.normal_panel}
                --coverage-files {input.cov_list}
                --assay {params.platform} --genome {params.genome} --force > {log} 2>&1
            """)

rule _panel_of_normals_purecn_output_intervals:
    input:
        intervals = str(rules._panel_of_normals_purecn_setinterval.output.intervals)
    output:
        intervals =  CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_intervals.txt"
    run:
        op.relative_symlink(input.intervals, output.intervals, in_module=True)

rule _panel_of_normals_purecn_output_gatk_intervals:
    input:
        gatk_intervals = str(rules._panel_of_normals_purecn_gatk_interval_list.output.gatk_intervals)
    output:
        gatk_intervals = CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_intervals_gatk.list"
    run:
        op.relative_symlink(input.gatk_intervals, output.gatk_intervals, in_module=True)

rule _panel_of_normals_purecn_output_gatk_intervals_targets:
    input:
        gatk_targets = str(rules._panel_of_normals_purecn_gatk_interval_list_targets.output.gatk_targets)
    output:
        gatk_targets = CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_intervals_gatk_targets.list"
    run:
        op.relative_symlink(input.gatk_targets, output.gatk_targets, in_module=True)

rule _panel_of_normals_purecn_output_mutect2_pon:
    input:
        pon = str(rules._panel_of_normals_purecn_mutect2_pon.output.pon),
    output:
        pon =  CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{capture_space}/{capture_space}_mutect2_pon.vcf.gz"
    run:
        op.relative_symlink(input.pon, output.pon, in_module=True)

rule _panel_of_normals_purecn_output_normal_panel_vcf:
    input:
        normal_panel = str(rules._panel_of_normals_purecn_merge_vcfs.output.normal_panel),
        normal_panel_tbi = str(rules._panel_of_normals_purecn_merge_vcfs.output.normal_panel_tbi)
    output:
        normal_panel =  CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{capture_space}/{capture_space}_normalpanel.vcf.gz",
        normal_panel_tbi = CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{capture_space}/{capture_space}_normalpanel.vcf.gz.tbi"
    run:
        op.relative_symlink(input.normal_panel, output.normal_panel, in_module=True)
        op.relative_symlink(input.normal_panel_tbi, output.normal_panel_tbi, in_module=True)

rule _panel_of_normals_purecn_output_database_cnvkit:
    input:
        mapping_bias = str(rules._panel_of_normals_purecn_database_cnvkit.output.mapping_bias)
    output:
        mapping_bias =  CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{capture_space}/purecn_cnvkit_normal/mapping_bias_{genome_build}_{capture_space}.rds"
    run:
        op.relative_symlink(input.mapping_bias, output.mapping_bias, in_module=True)

rule _panel_of_normals_purecn_output_database_denovo:
    input:
        mapping_bias = str(rules._panel_of_normals_purecn_database_denovo.output.mapping_bias)
    output:
        mapping_bias =  CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{capture_space}/purecn_denovo_normal/mapping_bias_{genome_build}_{capture_space}.rds",
        normal_db =  CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{capture_space}/purecn_denovo_normal/normalDB_{genome_build}_{capture_space}.rds"
    run:
        op.relative_symlink(input.mapping_bias, output.mapping_bias, in_module=True)
        op.relative_symlink(input.normal_db, output.normal_db, in_module=True)

# ---------------------------------------------------------- #
# rule all section
# ---------------------------------------------------------- #
rule _panel_of_normals_output_samples_tsv:
    input:
        tsv = str(rules._panel_of_normals_record_samples.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{capture_space}_samples_metadata.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module=True)

# Generates the target sentinels for each run, which generate the symlinks
rule _panel_of_normals_all:
    input:
        # record samples used
        expand(
            [str(rules._panel_of_normals_output_samples_tsv.output.tsv)],
            zip,
            seq_type=META["seq_type"],
            genome_build=META["genome_build"],
            capture_space=META["capture_space"]
        )
        # cnvkit
        expand(
            [
                str(rules._panel_of_normals_cnvkit_output_cnn.output.cnn),
                str(rules._panel_of_normals_cnvkit_output_tsv.output.tsv),
                str(rules._panel_of_normals_cnvkit_output_beds.output.target),
                str(rules._panel_of_normals_cnvkit_output_beds.output.antitarget)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=META["seq_type"],
            genome_build=META["genome_build"],
            capture_space=META["capture_space"]
        ),
        # cnvkit flat reference file
        expand(
            [
                str(rules._panel_of_normals_cnvkit_output_flat_ref.output.cnn),
                str(rules._panel_of_normals_cnvkit_output_flat_ref_beds.output.target),
                str(rules._panel_of_normals_cnvkit_output_flat_ref_beds.output.antitarget)
            ],
            zip,
            seq_type="capture",
            genome_build=CFG["options"]["cnvkit"]["flat_ref_combos"]["genome_builds"],
            capture_space=CFG["options"]["cnvkit"]["flat_ref_combos"]["capture_spaces"]
        ),
        # pureCN
        expand(
            [
                str(rules._panel_of_normals_purecn_output_intervals.output.intervals),
                str(rules._panel_of_normals_purecn_output_gatk_intervals.output.gatk_intervals),
                str(rules._panel_of_normals_purecn_output_gatk_intervals_targets.output.gatk_targets),
                str(rules._panel_of_normals_purecn_output_mutect2_pon.output.pn),
                str(rules._panel_of_normals_purecn_output_normal_panel_vcf.output.normal_panel),
                str(rules._panel_of_normals_purecn_output_normal_panel_vcf.output.normal_panel_tbi),
                str(rules._panel_of_normals_purecn_output_database_cnvkit.output.mapping_bias),
                str(rules._panel_of_normals_purecn_output_database_denovo.output.mapping_bias),
                str(rules._panel_of_normals_purecn_output_database_denovo.output.normal_db)
            ],
            zip
            seq_type="capture",
            genome_build=CFG["options"]["cnvkit"]["flat_ref_combos"]["genome_builds"],
            capture_space=CFG["options"]["cnvkit"]["flat_ref_combos"]["capture_spaces"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
