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
    name = "phase_vcf",
    version = "1.0",
    subdirectories = ["inputs", "clair3", "filter_clair3", "whatshap_phase_vcf","outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _promethion_input,
    _filter_clair3,
    _phase_vcf_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _promethion_input:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)

rule _clair3:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    params: 
        model = "/bin/models/" + CFG["clair3"]["model"],
        dir = CFG["dirs"]["clair3"] + "{seq_type}--{genome_build}/{sample_id}",
        platform = CFG["clair3"]["platform"]
    conda :
        CFG["conda_envs"]["clair3"]  
    threads:
        CFG["threads"]["clair3"]    
    resources: 
        mem_mb = CFG["mem_mb"]["clair3"]     
    log:
        stderr = CFG["logs"]["clair3"] + "{seq_type}--{genome_build}/{sample_id}/clair3.stderr.log"    
    output:
        vcf = CFG["dirs"]["clair3"] + "{seq_type}--{genome_build}/{sample_id}/phased_merge_output.vcf.gz"          
    shell:
         op.as_one_line("""run_clair3.sh  
            -b {input.bam} -f {input.fasta} -t {threads} -p {params.platform}
            -m ${{CONDA_PREFIX}}{params.model} -o {params.dir}
            --remove_intermediate_dir --enable_phasing
            2> {log.stderr} """)

rule _filter_clair3:
    input:
        vcf = str(rules._clair3.output.vcf)
    params:
        filter = CFG["clair3"]["filter"]
    output: 
        filtered = CFG["dirs"]["filter_clair3"] + "{seq_type}--{genome_build}/{sample_id}.filtered_phased_merged.vcf.gz",
        index = CFG["dirs"]["filter_clair3"] + "{seq_type}--{genome_build}/{sample_id}.filtered_phased_merged.vcf.gz.tbi"
    conda :
        CFG["conda_envs"]["clair3"]
    log:
        stderr = CFG["logs"]["filter_clair3"] + "{seq_type}--{genome_build}/{sample_id}/filter_clair3.stderr.log"      
    shell :
        "zcat {input.vcf} | awk '$6 > {params.filter} || $1 ~ /^#/' | bgzip > {output.filtered} && tabix -p vcf {output.filtered}  2> {log.stderr}"       

rule _whatshap_phase_vcf:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        vcf = CFG["inputs"]["vcf"] 
    conda :
        CFG["conda_envs"]["whatshap"] 
    resources: 
        mem_mb = CFG["mem_mb"]["whatshap"]        
    log:
        stderr = CFG["logs"]["whatshap_phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}/whatshap_phase_vcf.stderr.log"    
    output:
        vcf = temp(CFG["dirs"]["whatshap_phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.phased.vcf"),
        vcf_gz = CFG["dirs"]["whatshap_phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz",
        index = CFG["dirs"]["whatshap_phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz.tbi"
    shell:
        op.as_one_line(""" whatshap phase -o {output.vcf} --ignore-read-groups --reference={input.fasta} {input.vcf} {input.bam} 
            && bgzip -c {output.vcf} > {output.vcf_gz} && tabix -p vcf {output.vcf_gz} 2> {log.stderr} """)


rule _phase_vcf_output:
    input:
        vcf = str(rules._filter_clair3.output.filtered) if not CFG["inputs"]["vcf"] else str(rules._whatshap_phase_vcf.output.vcf_gz),
        index = str(rules._filter_clair3.output.index) if not CFG["inputs"]["vcf"] else str(rules._whatshap_phase_vcf.output.index)

    output:
        vcf = CFG["dirs"]["outputs"] + "phased_long_read_variants/{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz" if not CFG["inputs"]["vcf"] else CFG["dirs"]["outputs"] + "phased_short_read_variants/{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz",
        index = CFG["dirs"]["outputs"] + "phased_long_read_variants/{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz.tbi" if not CFG["inputs"]["vcf"] else CFG["dirs"]["outputs"] + "phased_short_read_variants/{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True),
        op.relative_symlink(input.index, output.index, in_module= True)        

# Generates the target sentinels for each run, which generate the symlinks
rule _phase_vcf_all:
    input:
        expand(
            [
                str(rules._phase_vcf_output.output.vcf),
                str(rules._phase_vcf_output.output.index)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG) 
