#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Daniel Nicorici
# Module Author:    Ryan Morin
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
# `CFG` is a shortcut to `config["lcr-modules"]["fusioncatcher"]`
CFG = op.setup_module(
    name = "fusioncatcher",
    version = "1.0",
    subdirectories = ["inputs", "fusioncatcher", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _fusioncatcher_input_fastq,
    _fusioncatcher_output_all,
    _fusioncatcher_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _fusioncatcher_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"],
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R2.fastq.gz",
    run:
        op.absolute_symlink(input.fastq_1, output.fastq_1)
        op.absolute_symlink(input.fastq_2, output.fastq_2)



rule _fusioncatcher_run:
    input:
        fastq_1 = str(rules._fusioncatcher_input_fastq.output.fastq_1),
        fastq_2 = str(rules._fusioncatcher_input_fastq.output.fastq_2),
        fastq_1_real = CFG["inputs"]["sample_fastq_1"], # Prevent premature deletion of fastqs marked as temp
        fastq_2_real = CFG["inputs"]["sample_fastq_2"],
        ref_path = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/fusioncatcher/fusioncatcher/data/human_v102/"
    output:
        complete = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/fusioncatcher.complete",
        summary = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/summary_candidate_fusions.txt",
        junk = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/junk-chimeras.txt",
        hg19 = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/final-list_candidate-fusion-genes.hg19.txt",
        hg38 = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/final-list_candidate-fusion-genes.txt",
        markdown = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/final-list_candidate-fusion-genes.caption.md.txt",
        virus = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/viruses_bacteria_phages.txt"
    log:
        stdout = CFG["logs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/fusioncatcher_run.stdout.log",
        stderr = CFG["logs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/fusioncatcher_run.stderr.log"
    conda:
        CFG["conda_envs"]["fusioncatcher"]
    threads:
        CFG["threads"]["fusioncatcher_run"]
    resources:
        **CFG["resources"]["fusioncatcher_run"]
    params:
        reference = CFG["options"]["reference_path"],
        out_dir = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{sample_id}/",
        # Base dir for fusioncatcher's ~150 GB of transient intermediates. "${TMPDIR:-/tmp}"
        # keeps them on node-local disk (no shared-scratch contention); a scratch path is
        # the spillover for libraries too deep for node-local. Only the small final outputs
        # are copied back to out_dir. The EXIT trap removes the work dir even on failure
        # (ENOSPC) so a partial 150 GB blob never orphans on the node.
        work_base = CFG["options"].get("work_base", "${TMPDIR:-/tmp}"),
        # Self-imposed wall-clock cap: fusioncatcher norm ~2h, but some libraries deadlock
        # early and pin the node for days. `timeout` SIGKILLs it at this limit so the job
        # exits (trap cleans up) and the node frees. (Not srun --time: that broke on this
        # partition - the UNLIMITED default exceeds the 50d MaxTime -> PartitionTimeLimit.)
        timeout = CFG["options"].get("timeout", "600m")
    shell:
        op.as_one_line("""
        work="{params.work_base}/fusioncatcher.{wildcards.seq_type}--{wildcards.genome_build}.{wildcards.sample_id}";
        cleanup() {{ rm -rf "$work"; }}; trap cleanup EXIT INT TERM;
        rm -rf "$work"; mkdir -p "$work";
        timeout -s KILL {params.timeout} fusioncatcher --no-update-check -d {params.reference} -i {input.fastq_1},{input.fastq_2}
        -o "$work" > {log.stdout} 2> {log.stderr};
        mkdir -p {params.out_dir};
        cp "$work/summary_candidate_fusions.txt" "$work/junk-chimeras.txt"
        "$work/final-list_candidate-fusion-genes.hg19.txt" "$work/final-list_candidate-fusion-genes.txt"
        "$work/final-list_candidate-fusion-genes.caption.md.txt" "$work/viruses_bacteria_phages.txt"
        {params.out_dir};
        touch {output.complete}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')

rule _fusioncatcher_output_all:
    input:
        complete = str(rules._fusioncatcher_run.output.complete),
        summary = str(rules._fusioncatcher_run.output.summary),
        junk = str(rules._fusioncatcher_run.output.junk),
        hg19 = str(rules._fusioncatcher_run.output.hg19),
        hg38 = str(rules._fusioncatcher_run.output.hg38),
        markdown = str(rules._fusioncatcher_run.output.markdown),
        virus = str(rules._fusioncatcher_run.output.virus)
    output:
        summary = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.summary.txt",
        junk = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.junk-chimeras.txt",
        hg19 = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.final-list_candidate-fusion-genes.hg19.txt",
        hg38 = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.final-list_candidate-fusion-genes.hg38.txt",
        markdown = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.final-list_column_explanation.md",
        virus = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.viruses_bacteria_phages.txt"
    run:
        op.relative_symlink(input.summary, output.summary, in_module= True)
        op.relative_symlink(input.junk, output.junk, in_module= True)
        op.relative_symlink(input.hg19, output.hg19, in_module= True)
        op.relative_symlink(input.hg38, output.hg38, in_module= True)
        op.relative_symlink(input.virus, output.virus, in_module= True)
        op.relative_symlink(input.markdown, output.markdown, in_module= True)



# Generates the target sentinels for each run, which generate the symlinks
rule _fusioncatcher_all:
    input:
        expand(
            [
                str(rules._fusioncatcher_output_all.output.summary)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            sample_id=CFG["runs"]["tumour_sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
