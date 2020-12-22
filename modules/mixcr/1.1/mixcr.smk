#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Anita dos Santos
# Module Author:    Anita dos Santos
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mixcr"]`
CFG = op.setup_module(
    name = "mixcr",
    version = "1.0",
    subdirectories = ["inputs", "mixcr", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _install_mixcr,
    _mixcr_input_fastq,
    _mixcr_output_txt,
    _mixcr_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mixcr_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"],
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R2.fastq.gz",
    run:
        op.relative_symlink(input.fastq_1, output.fastq_1)
        op.relative_symlink(input.fastq_2, output.fastq_2)

# Installs latest MiXCR release from github if the mixcr folder is not present yet
rule _install_mixcr:
    params:
        mixcr = CFG["inputs"]["mixcr_exec"]
    output: 
        complete = CFG["inputs"]["mixcr_exec"] + "/mixcr_dependencies_installed.success"
    shell:
        '''
        download_url=$(curl --silent "https://api.github.com/repos/milaboratory/mixcr/releases/latest" | grep '"browser_download_url":' | sed -E 's/.*\"([^\"]+)\".*/\\1/');
        mkdir -p {params.mixcr};

        if [ ! -f {params.mixcr}/mixcr ]; then
            wget -cO - $download_url > {params.mixcr}/mixcr.zip && unzip {params.mixcr}/mixcr.zip -d {params.mixcr}/ && rm {params.mixcr}/mixcr.zip;
            mv {params.mixcr}/mixcr*/* {params.mixcr}/ && rm -r {params.mixcr}/mixcr*/;
        fi

        touch  {output.complete};
        '''

# Run MiXCR rule
rule _mixcr_run:
    input:
        fastq_1 = rules._mixcr_input_fastq.output.fastq_1,
        fastq_2 = rules._mixcr_input_fastq.output.fastq_2,
        fastq_1_real = CFG["inputs"]["sample_fastq_1"], # Prevent premature deletion of fastqs marked as temp
        fastq_2_real = CFG["inputs"]["sample_fastq_2"],
        installed = rules._install_mixcr.output.complete
    output:
         txt = CFG["dirs"]["mixcr"] + "{seq_type}--{genome_build}/{sample_id}/mixcr.{sample_id}.clonotypes.ALL.txt",
         report = CFG["dirs"]["mixcr"] + "{seq_type}--{genome_build}/{sample_id}/mixcr.{sample_id}.report"
    log:
        stdout = CFG["logs"]["mixcr"] + "{seq_type}--{genome_build}/{sample_id}/mixcr_run.stdout.log",
        stderr = CFG["logs"]["mixcr"] + "{seq_type}--{genome_build}/{sample_id}/mixcr_run.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["mixcr_run"]),
        prefix = CFG["dirs"]["mixcr"] + "{seq_type}--{genome_build}/{sample_id}/mixcr.{sample_id}", 
        mixcr = CFG["inputs"]["mixcr_exec"] + "/mixcr"
    conda: CFG["conda_envs"]["java"]
    threads:
        CFG["threads"]["mixcr_run"]
    resources:
        **CFG["resources"]["mixcr_run"],
    message:
        "{params.mixcr}"    
    shell:
        op.as_one_line("""
        {params.mixcr} analyze shotgun -s hsa -t {threads} {params.opts} {input.fastq_1} {input.fastq_2} {params.prefix} > {log.stdout} 2> {log.stderr};
        touch "{output.txt}";
        """)
        

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mixcr_output_txt:
    input:
        txt = rules._mixcr_run.output.txt,
        report = rules._mixcr_run.output.report
    output:
        txt = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/mixcr.{sample_id}.clonotypes.ALL.txt",
        report = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/mixcr.{sample_id}.report"
    run:
        op.relative_symlink(input.txt, output.txt)
        op.relative_symlink(input.report, output.report)


# Generates the target sentinels for each run, which generate the symlinks
rule _mixcr_all:
    input:
        expand(
            [
                rules._install_mixcr.output.complete,
                rules._mixcr_output_txt.output.txt
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            sample_id=CFG["runs"]["tumour_sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
