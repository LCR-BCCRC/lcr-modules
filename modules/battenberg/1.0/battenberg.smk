#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["battenberg"]`
CFG = op.setup_module(
    name = "battenberg",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "battenberg", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _battenberg_input_bam,
    _battenberg_step_2,
    _battenberg_output_seg,
    _battenberg_all,


##### RULES #####
print(CFG)

# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _battenberg_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)

#Currently I think this rule has to be run twice for it to work properly because the conda environment is created here. 
rule _install_battenberg:
    output:
        complete = "config/envs/battenberg_dependencies_installed.success"
    conda:
        CFG["conda_envs"]["battenberg"]
    shell:
        """
        R -q -e 'BiocManager::install(c("devtools", "splines", "readr", "doParallel", "ggplot2", "RColorBrewer", "gridExtra", "gtools", "parallel"));' && 
        R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")' &&
        R -q -e 'devtools::install_github("morinlab/battenberg")' &&
        touch {output.complete}"""


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _run_battenberg:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
#        _install_battenberg.outputs.complete
    output:
        refit=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_refit_suggestion.txt",
        sub=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.txt",
        sub_1=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones_1.txt",
        ac=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleCounts.tab",
        mb=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantBAF.tab",
        mlrg=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR_gcCorrected.tab",
        mlr=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR.tab",
        nlr=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalLogR.tab",
        nb=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalBAF.tab"
    log:
        stdout = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stdout.log",
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stderr.log"
    params:
        reference_path = CFG["reference_path"],
        script = CFG["inputs"]["battenberg_script"],
        calc_sex_status = CFG["inputs"]["calc_sex_status"],
        x_chrom = CFG["options"]["x_chrom"],
        y_chrom = CFG["options"]["y_chrom"],
        chr_prefixed = CFG["options"]["chr_prefixed_reference"],
        out_dir = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}"
    conda:
        CFG["conda_envs"]["battenberg"]
    resources:
        mem_mb = CFG["mem_mb"]["battenberg"]
    threads:
        CFG["threads"]["battenberg"]
    resources: mem_mb = 200000
    shell:
        """sex=$({params.calc_sex_status} {input.normal_bam} {params.x_chrom} {params.y_chrom} | tail -1 | awk '{{if( $3/$2 > 0.1) print "male"; else print "female"}}') ;"""
        "Rscript {params.script} -t {wildcards.tumour_id} "
        "-n {wildcards.normal_id} --tb {input.tumour_bam} --nb {input.normal_bam} --cpu {threads} "
        "-o {params.out_dir} --sex $sex --reference {params.reference_path} {params.chr_prefixed}  > {log.stdout} 2> {log.stderr} "


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
#rule _battenberg_step_2:
#    input:
#        seg = rules._battenberg_step_1.output.seg
#    output:
#        seg = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.seg"
#    log:
#        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
#    params:
#        opts = CFG["options"]["step_2"]
#    shell:
#        "grep {params.opts} {input.seg} > {output.seg} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
#rule _battenberg_output_seg:
#    input:
#        seg = rules._battenberg_step_2.output.seg
#    output:
#        seg = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.seg"
#    run:
#        op.relative_symlink(input, output)


# Generates the target sentinels for each run, which generate the symlinks
rule _battenberg_all:
    input:
        expand(
            [
                rules._run_battenberg.output.refit,
                rules._install_battenberg.output.complete,
                #rules._battenberg_input_bam.output.bam,
                #rules._battenberg_output_seg.output.refit,
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
