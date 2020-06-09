#!/usr/bin/env snakemake


##### SETUP #####


# Import standard packages
import os
import gzip

# Import package with useful functions for developing analysis modules.
import modutils as md

# Make sure the `CFG` variable doesn't exist yet
assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

# Setup module and store module-specific configuration in `CFG`.
CFG = md.setup_module(
    config = config, 
    name = "manta", 
    version = "1.0",
    subdirs = ["inputs", "chrom_bed", "manta", "bedpe", "outputs"],
    req_references = ["genome_fasta", "genome_fasta_index", "main_chroms"]
)

# Define rules to be run locally when using a compute cluster.
localrules: 
    _manta_input_bam,
    _manta_generate_bed,
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
        sample_bam = CFG["inputs"]["sample_bam"]
    output:
        sample_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.sample_bam, output.sample_bam)
        op.relative_symlink(input.sample_bai, output.sample_bai)


# Generate BED file for main chromosomes to exclude small contigs from Manta run
rule _manta_generate_bed:
    input:
        fai = md.get_reference(CFG, "genome_fasta_index")
    output:
        bed = CFG["dirs"]["chrom_bed"] + "{genome_build}.main_chroms.bed"
    params:
        chroms = md.get_reference(CFG, "main_chroms")
    run:
        with open(input.fai) as fai, open(output.bed, "w") as bed:
            for line in fai:
                chrom, length, _, _, _ = line.rstrip("\n").split("\t")
                if chrom not in params.chroms:
                    continue
                bed_line = f"{chrom}\t0\t{length}\n"
                bed.write(bed_line)


# bgzip-compress and tabix-index the BED file to meet Manta requirement
rule _manta_index_bed:
    input:
        bed = rules._manta_generate_bed.output.bed
    output:
        bedz = CFG["dirs"]["chrom_bed"] + "{genome_build}.main_chroms.bed.gz"
    conda:
        CFG["conda_envs"].get("tabix") or "envs/tabix-0.2.6.yaml"
    shell:
        "bgzip {input.bed} && tabix {output.bedz}"


# Configures the manta workflow with the input BAM files and reference FASTA file.
rule _manta_configure:
    input:
        # Do not have a normal_bam as input in 'no_normal' mode
        unpack(md.switch_on_wildcard("pair_status", {
            "_default" : {"normal_bam": CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"},
            "no_normal" : {}
        })),
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        config = op.switch_on_wildcard("seq_type", CFG["switches"]["manta_config"]),
        bedz = rules._manta_index_bed.output.bedz
    output:
        runwf = CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        stdout = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_configure.stdout.log",
        stderr = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_configure.stderr.log"
    params:
        opts   = md.switch_on_column("seq_type", CFG["samples"], CFG["options"]["configure"]),
        fasta  = md.get_reference(CFG, "genome_fasta"),
        # Omit the normal BAM CLI argument if there is no normal
        normal_bam = md.switch_on_wildcard("pair_status", {
            "_default" : "--normalBam {input.normal_bam}",
            "no_normal" : ""
        }),
        # Use --bam for mrna data
        tumour_bam = md.switch_on_wildcard("seq_type", {
            "_default" : "--tumourBam {input.tumour_bam}",
            "mrna" : "--bam {input.tumour_bam}"
        })
    conda:
        CFG["conda_envs"].get("manta") or "envs/manta-1.6.0.yaml"
    shell:
        md.as_one_line("""
        configManta.py {params.opts} --referenceFasta {params.fasta} --callRegions {input.bedz}
        --runDir "$(dirname {output.runwf})" {params.tumour_bam} {params.normal_bam}
        --config {input.config} > {log.stdout} 2> {log.stderr}
        """)


# Launches manta workflow in parallel mode and deletes unnecessary files upon success.
checkpoint _manta_run:
    input:
        runwf = CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    output:
        vcf = CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/candidateSV.vcf.gz"
    log:
        stdout = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_run.stdout.log",
        stderr = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_run.stderr.log"
    params:
        variants_dir = CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/",
        opts   = CFG["options"]["manta"]
    conda:
        CFG["conda_envs"].get("manta") or "envs/manta-1.6.0.yaml"
    threads:
        CFG["threads"].get("manta") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("manta") or 1000
    shell:
        md.as_one_line("""
        {input.runwf} {params.opts} --jobs {threads} > {log.stdout} 2> {log.stderr}
            &&
        rm -rf "$(dirname {input.runwf})/workspace/"
        """)


# Fixes sample IDs in VCF header for compatibility with svtools vcftobedpe Otherwise, 
# manta uses the sample name from the BAM read groups, which may not be useful.
rule _manta_fix_vcf_ids:
    input:
        vcf = rules._manta_run.params.variants_dir + "{vcf_name}.vcf.gz"
    output:
        vcf = pipe(CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/{vcf_name}.with_ids.vcf")
    log:
        stderr = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_fix_vcf_ids.{vcf_name}.stderr.log"
    shell:
        md.as_one_line("""
        gzip -dc {input.vcf}
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
        vcf = pipe(CFG["dirs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/{vcf_name}.with_ids.with_vaf.vcf")
    log:
        stderr = CFG["logs"]["manta"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/manta_calc_vaf.{vcf_name}.stderr.log"
    conda:
        CFG["conda_envs"].get("calc_manta_vaf")
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
        CFG["conda_envs"].get("svtools") or "envs/svtools-0.5.1.yaml"
    threads:
        CFG["threads"].get("vcf_to_bedpe") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("vcf_to_bedpe") or 1000
    shell:
        "svtools vcftobedpe -i {input.vcf} > {output.bedpe} 2> {log.stderr}"


# Symlinks the VCF files
rule _manta_output_vcf:
    input:
        vcf = rules._manta_run.params.variants_dir + "{vcf_name}.vcf.gz"
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{vcf_name}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz"
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


def _get_manta_files(wildcards):
    """Request symlinks for all Manta VCF/BEDPE files.
    
    This function is required in conjunction with a Snakemake
    checkpoint because Manta produces different files based
    on whether it's run in paired mode or not and based on
    some parameters (like `--rna`). This function dynamically
    generates target symlinks for the raw VCF files and the
    processed BEDPE files based on what was actually produced.
    """
    no_bedpe = ["candidateSV"]
    manta_vcf = checkpoints._manta_run.get(**wildcards).output.vcf
    variants_dir = os.path.dirname(manta_vcf)
    all_files = os.listdir(variants_dir)
    vcf_files = [f for f in all_files if f.endswith(".vcf.gz")]
    vcf_names = [f.replace(".vcf.gz", "") for f in vcf_files]
    
    # Remove any empty VCF files from bedpe_targets
    vcf_filepaths = [os.path.join(variants_dir, f) for f in vcf_files]
    for vcf_name, vcf_filepath in zip(vcf_names, vcf_filepaths):
        with gzip.open(vcf_filepath, "rt") as vcf:
            content = [l for l in vcf.readlines() if not l.startswith("#")]
            if len(content) == 0:
                no_bedpe.append(vcf_name)
    
    vcf_targets = expand(rules._manta_output_vcf.output.vcf,
                         vcf_name=vcf_names, **wildcards)
    bedpe_targets = expand(rules._manta_output_bedpe.output.bedpe,
                           vcf_name=(set(vcf_names) - set(no_bedpe)), 
                           **wildcards)
    return vcf_targets + bedpe_targets


# Generates the target symlinks for each run
rule _manta_all_dispatch:
    input:
        _get_manta_files
    output:
        touch(CFG["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/.{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _manta_all:
    input:
        expand(rules._manta_all_dispatch.output, 
               zip,  # Run expand() with zip(), not product()
               seq_type=CFG["runs"]["tumour_seq_type"],
               genome_build=CFG["runs"]["tumour_genome_build"],
               tumour_id=CFG["runs"]["tumour_sample_id"],
               normal_id=CFG["runs"]["normal_sample_id"],
               pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk (including the samples and runs)
md.cleanup_module(CFG)

# Delete any local variables to avoid interfering with other code
del CFG
del _get_manta_files
