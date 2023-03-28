#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  n/a
# Module Author:    Krysta M Coyle PhD
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
# `CFG` is a shortcut to `config["lcr-modules"]["sv_repair"]`
CFG = op.setup_module(
    name = "sv_repair",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "window_bedpe","bed","fasta","tumor_reads","clustal", "sv_repair", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _sv_repair_input_bedpe,
    _sv_repair_step_2,
    _sv_repair_output_unk,
    _sv_repair_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
# TODO: If applicable, create second symlink to .crai file in the input function, to accomplish cram support
rule _sv_repair_input_bedpe:
    input:
        bedpe = CFG["inputs"]["sample_bedpe"]
        tumour_bam = CFG["inputs"]["tumour_bam"]
        tumour_bai = CFG["inputs"]["tumour_bai"]
    output:
        bedpe = CFG["dirs"]["inputs"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.bedpe"
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.{genome_build}.bam"
        tumour_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.{genome_build}.bam.bai
    group: 
        "input_and_step_1"
    run:
        op.absolute_symlink(input.bedpe, output.bedpe),
        op.absolute_symlink(input.tumour_bam, output.tumour_bam),
        op.absolute_symlink(input.tumour_bai, output.tumour_bai)


## Add windows to bedpe
rule _sv_repair_step_1:rule add_windows_bedpe:
	input:
		bedpe = str(rules.input_bedpe.output.bedpe)
	output:
		bedpe = CFG["dirs"]["window_bedpe"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.bedpe"
	log:
		stdout = CFG["logs"]["window_bedpe"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/add_windows.stdout.log"
	params:
        CFG["options"]["window_size"]
	conda:
		CFG["conda_envs"]["readr"]
	script:
		"../../src/R/bedpe.windows.R"

## Convert bedpe to bed
rule bedpe_to_bed:
	input:
		bedpe = str(rules.add_windows_bedpe.output.bedpe)
	output:
		bed = OUTPUT_DIR + "03-bed/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.bed"
	log:
		stdout = OUTPUT_DIR + "00-logs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/bedpe.to.bed.stdout.log"
	conda:
		CFG["conda_envs"]["readr"]
	script:
		"../../src/R/bedpe.to.bed.R"


rule get_fasta:
	conda:
        CFG["conda_envs"]["bedtools"]
	input:
		bed = str(rules.bedpe_to_bed.output.bed),
		genome_fasta = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/ref/lcr-modules-references-STABLE/genomes/{genome_build}/genome_fasta/genome.fa"
	log:
		stdout = OUTPUT_DIR + "00-logs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/get_fasta.stdout.log"
	output:
		reference_fasta = OUTPUT_DIR + "04-fasta/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/reference.fa"
	shell:
		op.as_one_line("""
			bedtools getfasta -s -fo {output.reference_fasta} -fi {input.genome_fasta} -bed {input.bed} 
		""")

# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _sv_repair_output_unk:
    input:
        unk = str(rules._sv_repair_step_2.output.unk)
    output:
        unk = CFG["dirs"]["outputs"] + "unk/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.unk"
    run:
        op.relative_symlink(input.unk, output.unk, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _sv_repair_all:
    input:
        expand(
            [
                str(rules._sv_repair_output_unk.output.unk),
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
