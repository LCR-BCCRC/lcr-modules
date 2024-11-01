
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
# `CFG` is a shortcut to `config["lcr-modules"]["sniffles"]`
CFG = op.setup_module(
    name = "sniffles",
    version = "1.0",
    subdirectories = ["inputs", "sniffles", "bedpe", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _promethion_bam,
    _sniffles_output,
    _sniffles_all


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
        


rule _sniffles:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["sniffles"] + "{seq_type}--{genome_build}/{sample_id}.vcf"
    conda: 
        CFG["conda_envs"]["sniffles"]
    resources: 
       mem_mb = CFG["mem_mb"]["sniffles"]
    threads:    
        CFG["threads"]["sniffles"]  
    conda :
        CFG["conda_envs"]["sniffles"]
    log:
        CFG["logs"]["sniffles"] + "{seq_type}--{genome_build}/{sample_id}/sniffles.log"               
    shell:
        op.as_one_line('''
        sniffles -t {threads} --input {input.bam} 
        --vcf {output.vcf} --reference {input.fasta} 
        --non-germline 
        2>&1 | tee -a {log}
        ''') 

rule _sniffles_vcf_to_bedpe:
    input:
        vcf = str(rules._sniffles.output.vcf)
    output:
        bedpe = CFG["dirs"]["bedpe"] + "{seq_type}--{genome_build}/{sample_id}.bedpe"
    log:
        stderr = CFG["logs"]["bedpe"] + "{seq_type}--{genome_build}/{sample_id}/sniffles_vcf_to_bedpe.stderr.log"
    conda:
        CFG["conda_envs"]["svtools"]
    threads:
        CFG["threads"]["vcf_to_bedpe"]
    resources: 
        mem_mb = CFG["mem_mb"]["vcf_to_bedpe"]
    shell:
        "svtools vcftobedpe -i {input.vcf} > {output.bedpe} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sniffles_output:
    input:
        vcf = str(rules._sniffles.output.vcf),
        bedpe = str(rules._sniffles_vcf_to_bedpe.output.bedpe)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.vcf",
        bedpe = CFG["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/{sample_id}.bedpe"

    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True),
        op.relative_symlink(input.bedpe, output.bedpe, in_module= True)




# Generates the target sentinels for each run, which generate the symlinks
rule _sniffles_all:
    input:
        expand(
            [
                str(rules._sniffles_output.output.vcf),
                str(rules._sniffles_output.output.bedpe)

            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)  
