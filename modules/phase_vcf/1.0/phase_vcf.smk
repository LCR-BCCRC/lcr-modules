# Original Author:  Nicole Thomas
# Module Author:    Nicole Thomas
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import inspect
import hashlib
import glob

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
# `CFG` is a shortcut to `config["lcr-modules"]["phase_vcf"]`
CFG = op.setup_module(
    name = "clair3",
    version = "1.0",
    subdirectories = ["inputs", "clair3", "filter_clair3", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _clair3_input,
    _clair3_get_models,
    _clair3_output,
    _clair3_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _clair3_input:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai", 
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)

# Download the models trained for clair3 variant calling
rule _clair3_get_models: 
    output: 
        model = directory(CFG["dirs"]["inputs"] + "clair3_models/{model}")
    params: 
        url = "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/{model}.tar.gz" 
    log: 
        stderr = CFG["logs"]["inputs"] + "{model}.log"
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        op.as_one_line("""
            wget -qO- {params.url} | tar -xvzf - -C $(dirname {output.model}) 2> {log.stderr}
        """)
    
def _clair3_which_model(wildcards): 
    CFG = config["lcr-modules"]["clair3"]
    this_sample = op.filter_samples(
        CFG["samples"], 
        sample_id = wildcards.sample_id, 
        seq_type = wildcards.seq_type, 
        genome_build = wildcards.genome_build
        )
    clair3_model = CFG["options"]["model"][this_sample["basecalling_model"].tolist()[0]]
    model = expand(
        CFG["dirs"]["inputs"] + "clair3_models/{model}", 
        model = clair3_model
    )
    return model

rule _clair3_run:
    input:
        bam = str(rules._clair3_input.output.bam),
        bai = str(rules._clair3_input.output.bai),
        cram = str(rules._clair3_input.output.crai), 
        model = _clair3_which_model,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    params: 
        dir = CFG["dirs"]["clair3"] + "{seq_type}--{genome_build}/{sample_id}",
        platform = CFG["options"]["platform"], 
        options = CFG["options"]["clair3"]
    conda:
        CFG["conda_envs"]["clair3"]  
    threads:
        CFG["resources"]["clair3"]["threads"]    
    resources: 
        bam = 1,
        mem_mb = CFG["resources"]["clair3"]["mem_mb"]    
    log:
        stderr = CFG["logs"]["clair3"] + "{seq_type}--{genome_build}/{sample_id}/clair3.log"    
    output:
        vcf = CFG["dirs"]["clair3"] + "{seq_type}--{genome_build}/{sample_id}/phased_merge_output.vcf.gz"          
    shell:
        op.as_one_line("""run_clair3.sh  
            -b {input.bam} -f {input.fasta} -t {threads} -p {params.platform}
            -m {input.model} -o {params.dir} {params.options}
            --remove_intermediate_dir --enable_phasing
            2>&1 {log} """)

rule _clair3_filter:
    input:
        vcf = str(rules._clair3_run.output.vcf)
    params:
        filter = CFG["options"]["filter"]
    output: 
        filtered = CFG["dirs"]["filter_clair3"] + "{seq_type}--{genome_build}/{sample_id}.filtered_phased_merged.vcf.gz",
        index = CFG["dirs"]["filter_clair3"] + "{seq_type}--{genome_build}/{sample_id}.filtered_phased_merged.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["clair3"]
    resources: 
        mem_mb = CFG["resources"]["filter"]["mem_mb"]
    threads: CFG["resources"]["filter"]["threads"]
    log:
        stderr = CFG["logs"]["filter_clair3"] + "{seq_type}--{genome_build}/{sample_id}/filter_clair3.stderr.log"      
    shell:
        "zcat {input.vcf} | awk '$6 > {params.filter} || $1 ~ /^#/' | bgzip > {output.filtered} && tabix -p vcf {output.filtered}  2> {log.stderr}"       


rule _clair3_output:
    input:
        vcf = str(rules._clair3_filter.output.filtered),
        index = str(rules._clair3_filter.output.index)
    output:
        vcf = CFG["dirs"]["outputs"] + "phased_long_read_variants/{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz",
        index = CFG["dirs"]["outputs"] + "phased_long_read_variants/{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True),
        op.relative_symlink(input.index, output.index, in_module= True)        

# Generates the target sentinels for each run, which generate the symlinks
rule _clair3_all:
    input:
        expand(
            [
                str(rules._clair3_output.output.vcf),
                str(rules._clair3_output.output.index)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG) 
