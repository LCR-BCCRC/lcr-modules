# Original Author:  Nicole Thomas
# Module Author:    Nicole Thomas
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import inspect

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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["cutesv"]`
CFG = op.setup_module(
    name = "cutesv",
    version = "1.0",
    subdirectories = ["inputs", "cutesv", "cutesv_working", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _promethion_bam,
    _cutesv_output,
    _cutesv_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _promethion_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        

rule _cutesv:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    threads:
         CFG["threads"]["cutesv"]
    resources: 
        mem_mb = CFG["mem_mb"]["cutesv"]        
    output:
        dir = temp(directory(CFG["dirs"]["cutesv_working"] + "{seq_type}--{genome_build}/{sample_id}")),
        vcf = CFG["dirs"]["cutesv"] + "{seq_type}--{genome_build}/{sample_id}.vcf" 
    conda :
        CFG["conda_envs"]["cutesv"] 
    log:
        stderr = CFG["logs"]["cutesv"] + "{seq_type}--{genome_build}/{sample_id}/cutesv.stderr.log"             
    shell:
         op.as_one_line("""
            mkdir {output.dir} &&
            cuteSV {input.bam} {input.fasta} {output.vcf} {output.dir}
            --threads {threads} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 
            --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3""")    


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _cutesv_output:
    input:
        cutesv = str(rules._cutesv.output.vcf)
    output:
        cutesv = CFG["dirs"]["outputs"] + "cutesv/{seq_type}--{genome_build}/{sample_id}.vcf"

    run:
        op.relative_symlink(input.cutesv, output.cutesv, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _cutesv_all:
    input:
        expand(
            [
                str(rules._cutesv_output.output.cutesv)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)                  








