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

#this preserves the variable when using lambda functions
_battenberg_CFG = CFG

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _battenberg_input_bam,
    _battenberg_step_2,
    _battenberg_output_seg,
    _battenberg_all,


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _battenberg_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)

# Installs the Battenberg R dependencies and associated software (impute2, alleleCounter)
# Currently I think this rule has to be run twice for it to work properly because the conda environment is created here. 
rule _install_battenberg:
    output:
        complete = "config/envs/battenberg_dependencies_installed.success"
    conda:
        CFG["conda_envs"]["battenberg"]
    shell:
        """
        R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")' && ##move some of this to config?
        R -q -e 'devtools::install_github("morinlab/battenberg")' &&              ##move some of this to config?
        touch {output.complete}"""


#This rule runs the entire Battenberg pipeline
rule _run_battenberg:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        tumour_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam.bai",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        normal_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam.bai",
        installed = "config/envs/battenberg_dependencies_installed.success"
    output:
        refit=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_refit_suggestion.txt",
        sub=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.txt",
        sub1=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones_1.txt",
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
        reference_path = lambda w: _battenberg_CFG["reference_path"][w.genome_build],
        script = CFG["inputs"]["battenberg_script"],
        calc_sex_status = CFG["inputs"]["calc_sex_status"],
        x_chrom = lambda w: _battenberg_CFG["options"]["x_chrom"][w.genome_build],
        y_chrom = lambda w: _battenberg_CFG["options"]["y_chrom"][w.genome_build],
        chr_prefixed = lambda w: _battenberg_CFG["options"]["chr_prefixed_reference"][w.genome_build],
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
        "-n {wildcards.normal_id} --tb {input.tumour_bam} --nb {input.normal_bam}  "
        "-o {params.out_dir} --sex $sex --reference {params.reference_path} {params.chr_prefixed} --cpu {threads} > {log.stdout} 2> {log.stderr} "

# Converts the subclones.txt (best fit) and subclones_1.txt (second best fit) to igv-friendly SEG files. 
rule _battenberg_to_igv_seg:
    input:
        sub = rules._run_battenberg.output.sub,
        sub1 = rules._run_battenberg.output.sub1,
        cnv2igv = CFG["inputs"]["cnv2igv"]
    output:
        seg = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.igv.seg",
        seg1 = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones_1.igv.seg"
    log:
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_seg2igv.stderr.log"
    threads: 1
    shell:
        op.as_one_line("""
        python {input.cnv2igv} --mode battenberg --sample {wildcards.tumour_id} 
        {input.sub} > {output.seg} 2> {log.stderr} &&
        python {input.cnv2igv} --mode battenberg --sample {wildcards.tumour_id} 
        {input.sub1} > {output.seg1} 2>> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _battenberg_output_seg:
    input:
        seg = rules._battenberg_to_igv_seg.output.seg,
        seg1 = rules._battenberg_to_igv_seg.output.seg1
    output:
        seg = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_subclones.igv.seg",
        seg1 = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_subclones_1.igv.seg"
    run:
        op.relative_symlink(input.seg, output.seg)  #note: I think there's a bug in the cookicutter. This line is op.relative_symlink(input, output) but presumably should include the names
        op.relative_symlink(input.seg1, output.seg1)


# Generates the target sentinels for each run, which generate the symlinks
rule _battenberg_all:
    input:
        expand(
            [
                rules._run_battenberg.output.refit,
                rules._install_battenberg.output.complete,
                rules._battenberg_to_igv_seg.output.seg,
                rules._battenberg_to_igv_seg.output.seg1,
                rules._battenberg_output_seg.output.seg,
                rules._battenberg_output_seg.output.seg1,
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
