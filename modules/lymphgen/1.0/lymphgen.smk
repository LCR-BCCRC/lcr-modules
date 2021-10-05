#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  NA
# Module Author:    Chris Rushton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["lymphgen"]`
CFG = op.setup_module(
    name = "lymphgen",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "lymphgen", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _lymphgen_input_maf,
    _lymphgen_step_2,
    _lymphgen_output_txt,
    _lymphgen_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _lymphgen_input_maf:
    input:
        maf = CFG["inputs"]["sample_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{sample_id}.maf"
    run:
        op.relative_symlink(input.maf, output.maf)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _lymphgen_step_1:
    input:
        tumour_maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}.maf",
        normal_maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{normal_id}.maf",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        txt = CFG["dirs"]["lymphgen"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.txt"
    log:
        stdout = CFG["logs"]["lymphgen"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
        stderr = CFG["logs"]["lymphgen"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["step_1"]
    resources:
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        op.as_one_line("""
        <TODO> {params.opts} --tumour {input.tumour_maf} --normal {input.normal_maf}
        --ref-fasta {input.fasta} --output {output.txt} --threads {threads}
        > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _lymphgen_step_2:
    input:
        txt = str(rules._lymphgen_step_1.output.txt)
    output:
        txt = CFG["dirs"]["lymphgen"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.txt"
    log:
        stderr = CFG["logs"]["lymphgen"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.txt} > {output.txt} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _lymphgen_output_txt:
    input:
        txt = str(rules._lymphgen_step_2.output.txt)
    output:
        txt = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.txt"
    run:
        op.relative_symlink(input.txt, output.txt)


# Generates the target sentinels for each run, which generate the symlinks
rule _lymphgen_all:
    input:
        expand(
            [
                str(rules._lymphgen_output_txt.output.txt),
                # TODO: If applicable, add other output rules here
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
