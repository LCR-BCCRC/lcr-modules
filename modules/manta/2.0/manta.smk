#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Snakefile Author:    Bruno Grande
# Module Author:                Bruno Grande
# Additional Contributors:      N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules.
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`.
CFG = op.setup_module(
    name = "manta", 
    version = "1.0",
    subdirectories = ["inputs", "chrom_bed", "manta", "calc_vaf", "bedpe", "outputs"]
)

# Define rules to be run locally when using a compute cluster.
localrules: 
    _manta_input_bam,
    _manta_input_bam_none,
    _manta_index_bed,
    _manta_configure,
    _manta_output_bedpe,
    _manta_output_vcf, 
    _manta_all_dispatch,
    _manta_all


##### RULES #####


# Symlinks the input BAM files into the module output directory (under '00-inputs/').
rule _manta_input_bam:
    input:
        sample_bam = CFG["inputs"]["sample_bam"],
        sample_bai = CFG["inputs"]["sample_bai"]
    output:
        sample_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        sample_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.sample_bam, output.sample_bam)
        op.relative_symlink(input.sample_bai, output.sample_bai)


# Create empty file for "no normal" runs (but this is ultimately omitted from downstream rules)
rule _manta_input_bam_none:
    output:
        touch(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/None.bam")


# bgzip-compress and tabix-index the BED file to meet Manta requirement
rule _manta_index_bed:
    input:
        bed = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.bed")
    output:
        bedz = CFG["dirs"]["chrom_bed"] + "{genome_build}.main_chroms.bed.gz"
    conda:
        CFG["conda_envs"]["tabix"]
    shell:
        op.as_one_line("""
        bgzip -c {input.bed} > {output.bedz}
            &&
        tabix {output.bedz}
        """)


# Configures the manta workflow with the input BAM files and reference FASTA file.
rule _manta_configure:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome.fa"),
        config = op.switch_on_wildcard("seq_type", CFG["switches"]["manta_config"]),
        bedz = rules._manta_index_bed.output.bedz
    output:
        runwf = CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        stdout = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_configure.stdout.log",
        stderr = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_configure.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["configure"]),
        normal_bam_arg = op.switch_on_wildcard("pair_status", CFG["switches"]["normal_bam_arg"]),
        tumour_bam_arg = op.switch_on_wildcard("seq_type", CFG["switches"]["tumour_bam_arg"])
    conda:
        CFG["conda_envs"]["manta"]
    shell:
        op.as_one_line("""
        configManta.py {params.opts} --referenceFasta {input.fasta} --callRegions {input.bedz}
        --runDir "$(dirname {output.runwf})" {params.tumour_bam_arg} {params.normal_bam_arg}
        --config {input.config} > {log.stdout} 2> {log.stderr}
        """)


# Launches manta workflow in parallel mode and deletes unnecessary files upon success.
checkpoint _manta_run:
    input:
        runwf = CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    output:
        variants_dir = directory(CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/"),
    log:
        stdout = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_run.stdout.log",
        stderr = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_run.stderr.log"
    params:
        opts = CFG["options"]["manta"]
    conda:
        CFG["conda_envs"]["manta"]
    threads:
        CFG["threads"]["manta"]
    resources: 
        mem_mb = CFG["mem_mb"]["manta"]
    shell:
        op.as_one_line("""
        {input.runwf} {params.opts} --jobs {threads} > {log.stdout} 2> {log.stderr}
            &&
        rm -rf "$(dirname {input.runwf})/workspace/"
        """)


# Fixes sample IDs in VCF header for compatibility with svtools vcftobedpe Otherwise, 
# manta uses the sample name from the BAM read groups, which may not be useful.
rule _manta_fix_vcf_ids:
    input:
        variants_dir = rules._manta_run.output.variants_dir
    output:
        vcf = CFG["dirs"]["calc_vaf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.with_ids.vcf"
    log:
        stderr = CFG["logs"]["calc_vaf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_fix_vcf_ids.{vcf_name}.stderr.log"
    threads:
        CFG["threads"]["fix_vcf_ids"]
    resources: 
        mem_mb = CFG["mem_mb"]["fix_vcf_ids"]
    shell:
        op.as_one_line("""
        gzip -dc {input.variants_dir}/{wildcards.vcf_name}.vcf.gz"
            |
        awk 'BEGIN {{FS=OFS="\\t"}}
        $1 == "#CHROM" && $10 != "" && $11 != "" {{$10="{wildcards.normal_id}"}}
        $1 == "#CHROM" && $10 != "" && $11 == "" {{$10="{wildcards.tumour_id}"}}
        $1 == "#CHROM" && $11 != "" {{$11="{wildcards.tumour_id}"}}
        {{print $0}}' > {output.vcf} 2> {log.stderr}
        """)


# Calculates the tumour and normal variant allele fraction (VAF) from the allele counts
# and creates new fields in the INFO column for convenience.
rule _manta_calc_vaf:
    input:
        vcf  = rules._manta_fix_vcf_ids.output.vcf,
        cvaf = CFG["inputs"]["calc_manta_vaf"]
    output:
        vcf = CFG["dirs"]["calc_vaf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.with_ids.with_vaf.vcf"
    log:
        stderr = CFG["logs"]["calc_vaf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_calc_vaf.{vcf_name}.stderr.log"
    conda:
        CFG["conda_envs"]["calc_manta_vaf"]
    threads:
        CFG["threads"]["calc_vaf"]
    resources: 
        mem_mb = CFG["mem_mb"]["calc_vaf"]
    shell:
        "{input.cvaf} {input.vcf} > {output.vcf} 2> {log.stderr}"


# Converts the VCF file into a more tabular BEDPE file, which is easier to handle in R
# and automatically pairs up breakpoints for interchromosomal events.
rule _manta_vcf_to_bedpe:
    input:
        vcf  = rules._manta_calc_vaf.output.vcf
    output:
        bedpe = CFG["dirs"]["bedpe"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.bedpe"
    log:
        stderr = CFG["logs"]["bedpe"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_vcf_to_bedpe.{vcf_name}.stderr.log"
    conda:
        CFG["conda_envs"]["svtools"]
    threads:
        CFG["threads"]["vcf_to_bedpe"]
    resources: 
        mem_mb = CFG["mem_mb"]["vcf_to_bedpe"]
    shell:
        "svtools vcftobedpe -i {input.vcf} > {output.bedpe} 2> {log.stderr}"


# Symlinks the VCF files with FORMAT columns
rule _manta_output_vcf_with_vafs:
    input:
        vcf = rules._manta_calc_vaf.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{vcf_name}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)


# Symlinks the VCF files without FORMAT columns
rule _manta_output_vcf_without_vafs:
    input:
        vcf = rules._manta_fix_vcf_ids.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{vcf_name}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)


# Symlinks the final BEDPE files
rule _manta_output_bedpe:
    input:
        bedpe = rules._manta_vcf_to_bedpe.output.bedpe
    output:
        bedpe = CFG["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/{vcf_name}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.bedpe"
    run:
        op.relative_symlink(input.bedpe, output.bedpe)


def _manta_predict_manta_output(wildcards):
    """Request symlinks for all Manta VCF/BEDPE files.
    
    This function is required in conjunction with a Snakemake
    checkpoint because Manta produces different files based
    on whether it's run in paired mode or not and based on
    some parameters (like `--rna`). This function dynamically
    generates target symlinks for the raw VCF files and the
    processed BEDPE files based on what was actually produced.
    """

    # Get module config (since `CFG` is temporary)
    CFG = config["lcr-modules"]["manta"]

    # Standard VCF outputs
    vcf_names = ["candidateSV", "candidateSmallIndels"]

    # Create switch to emulate what is done in `_manta_configure`
    configure_switch = op.switch_on_wildcard(
        "seq_type", CFG["options"]["configure"]
    )
    # And use the switch right away to get the command-line parameters used
    configure_params = configure_switch(wildcards)
    # Check if run with '--rna'
    if "--rna" in configure_params:
        vcf_names.append("rnaSV")

    # Check it run in `no_normal` mode
    if wildcards.pair_status == "no_normal":
        vcf_names.append("tumorSV")
    else:
        vcf_names.append("diploidSV")
        vcf_names.append("somaticSV")

    # VCF files without VAF information (due to lack of FORMAT columns)
    no_vafs = ["candidateSV", "candidateSmallIndels"]
    vcf_names_with_vafs = set(vcf_names) - set(no_vafs)
    vcf_names_without_vafs = set(vcf_names) & set(no_vafs)

    # Request the output files based on whether VAF info is available
    outputs_with_vafs = expand(
        [
            rules._manta_output_vcf_with_vafs.output.vcf,
            rules._manta_output_bedpe.output.bedpe,
        ],
        vcf_name=vcf_names_with_vafs,
        **wildcards
    )
    outputs_without_vafs = expand(
        rules._manta_output_vcf_without_vafs.output.vcf,
        vcf_name=vcf_names_without_vafs,
        **wildcards
    )

    return outputs_with_vafs + outputs_without_vafs


# Generates the target symlinks for each run depending on the Manta output VCF files
rule _manta_dispatch:
    input:
        _manta_predict_output
    output:
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _manta_all:
    input:
        expand(
            [
                rules._manta_dispatch.output.dispatched, 
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
