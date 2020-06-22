#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Bruno Grande
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["vcf2maf"]`
CFG = op.setup_module(
    name = "vcf2maf",
    version = "1.0",
    subdirectories = ["inputs", "vcf2maf", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _vcf2maf_input_vcf,
    _vcf2maf_output_maf,
    _vcf2maf_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _vcf2maf_input_vcf:
    input:
        vcf = op.switch_on_wildcards("caller", CFG["inputs"]["sample_vcf"])
    output:
        vcf = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}/{vartype}.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _vcf2maf_run:
    input:
        vcf = rules._vcf2maf_input_vcf.output.vcf,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        vep_cache = CFG["inputs"]["vep_cache"]
    output:
        maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}/{vartype}.maf"
    log:
        stdout = CFG["logs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}_{vartype}.stdout.log",
        stderr = CFG["logs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}_{vartype}.stderr.log"
    params:
        opts = CFG["options"]["vcf2maf"]
    conda:
        CFG["conda_envs"]["vcf2maf"]
    threads:
        CFG["threads"]["vcf2maf"]
    resources:
        mem_mb = CFG["mem_mb"]["vcf2maf"]
    shell:
        op.as_one_line("""
        vepPATH=$(dirname $(which vep))/../share/variant-effect-predictor*;
        vcf2maf.pl 
        --input-vcf {input.vcf} 
        --output-maf {output.maf} 
        --tumor-id {wildcards.tumour_id} --normal_id {wildcards.normal_id}
        --vcf-tumor-id TUMOR --vcf-normal_id NORMAL
        --ref-fasta {params.fasta}
        --vep-data {params.vep} --vep-path $vepPATH {params.opts} 2> {log}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _vcf2maf_output_maf:
    input:
        maf = rules._vcf2maf_run.output.maf
    output:
        maf = CFG["dirs"]["outputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}/{vartype}.maf"
    run:
        op.relative_symlink(input, output)


def _vcf2maf_get_output(wildcards):
    CFG = config["lcr-modules"]["vcf2maf"]



# Generates the target sentinels for each run, which generate the symlinks
rule _vcf2maf_all:
    input:
        expand(rules._vcf2maf_output_maf.output.maf,
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
