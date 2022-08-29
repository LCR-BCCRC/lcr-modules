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
# `CFG` is a shortcut to `config["lcr-modules"]["nanomethphase"]`
CFG = op.setup_module(
    name = "nanomethphase",
    version = "1.0",
    subdirectories = ["inputs", "methyl_call", "clair3", "filter_clair3", "whatshap_phase_vcf","nanomethphase","outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _promethion_input,
    _nanomethphase,
    _nanomethphase_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _promethion_input:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"],
        mc = CFG["inputs"]["meth_calls"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        mc = CFG["dirs"]["inputs"] + "meth_calls/{seq_type}--{genome_build}/{sample_id}.calls.tsv.gz"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.mc, output.mc)

rule _methyl_call_processor:
    input:
        calls = CFG["dirs"]["inputs"] + "meth_calls/{seq_type}--{genome_build}/{sample_id}.calls.tsv.gz"
    output: 
        bed = CFG["dirs"]["methyl_call"] + "{seq_type}--{genome_build}/{sample_id}.bed.gz"
    threads:
        CFG["threads"]["methyl_call"]
    conda:
        CFG["conda_envs"]["nanomethphase"] 
    resources: 
        mem_mb = CFG["mem_mb"]["methyl_call"]        
    log:
        stderr = CFG["logs"]["methyl_call"] + "{seq_type}--{genome_build}/{sample_id}/methyl_call_processor.stderr.log"             
    shell:
        op.as_one_line("""
        nanomethphase methyl_call_processor -mc {input.calls} --threads {threads} | 
        sort -k1,1 -k2,2n -k3,3n | bgzip > {output.bed} && tabix -p bed {output.bed}
        """)   

rule _clair3:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    params: 
        model = CFG["clair3"]["model"],
        dir = CFG["dirs"]["clair3"] + "{seq_type}--{genome_build}/{sample_id}"
    threads:
        CFG["threads"]["clair3"]    
    conda :
        CFG["conda_envs"]["clair3"]  
    resources: 
        mem_mb = CFG["mem_mb"]["clair3"]     
    log:
        stderr = CFG["logs"]["clair3"] + "{seq_type}--{genome_build}/{sample_id}/clair3.stderr.log"    
    output:
        vcf = CFG["dirs"]["clair3"] + "{seq_type}--{genome_build}/{sample_id}/phased_merge_output.vcf.gz"       
    shell:
         op.as_one_line("""run_clair3.sh  
            -b {input.bam} -f {input.fasta} -t {threads} -p "ont"
            -m {params.model} -o {params.dir}
            --remove_intermediate_dir --enable_phasing """)

rule _filter_clair3:
    input:
        vcf = str(rules._clair3.output.vcf)
    output: 
        filtered = CFG["dirs"]["filter_clair3"] + "{seq_type}--{genome_build}/{sample_id}.filtered_phased_merged.vcf.gz"
    log:
        stderr = CFG["logs"]["filter_clair3"] + "{seq_type}--{genome_build}/{sample_id}/filter_clair3.stderr.log"      
    shell :
        "zcat {input.vcf} | awk '$6 > 20' | gzip > {output.filtered}"       

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
        vcf = temp(CFG["dirs"]["whatshap_phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.output.vcf"),
        vcf_gz = CFG["dirs"]["whatshap_phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.output.vcf.gz"
    shell:
        op.as_one_line(""" whatshap phase --ignore-read-groups --sample={wildcards.sample_id} -o {output.vcf} 
            --reference={input.fasta} {input.vcf} {input.bam} 
            && bgzip -c {output.vcf} > {output.vcf_gz} """)


rule _nanomethphase:
    input:  
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        vcf = str(rules._whatshap_phase_vcf.output.vcf_gz) if CFG["inputs"]["vcf"] else str(rules._filter_clair3.output.filtered),
        mc = str(rules._methyl_call_processor.output.bed)
    threads:
        CFG["threads"]["nanomethphase"]
    params:
        prefix =  CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}"
    conda :
        CFG["conda_envs"]["nanomethphase"] 
    log:
        stderr = CFG["logs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/nanomethphase.stderr.log"         
    output:
        HP1_bis_bam = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP1_Converted2Bisulfite.bam",
        HP2_bis_bam = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP2_Converted2Bisulfite.bam",
        HP1_bis_bai = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP1_Converted2Bisulfite.bam.bai",
        HP2_bis_bai = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP2_Converted2Bisulfite.bam.bai",
        HP1_freq = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP1_MethylFrequency.tsv",
        HP2_freq = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP2_MethylFrequency.tsv"
    shell:
        op.as_one_line("""
        nanomethphase phase -mc {input.mc} -of methylcall,bam2bis -o {params.prefix} 
        -b {input.bam} -r {input.fasta} -v {input.vcf}
        -t {threads} -is && samtools index -@ {threads} {output.HP1_bis_bam} &&
        samtools index -@ {threads} {output.HP2_bis_bam}
        """)  


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _nanomethphase_output:
    input:
        HP1_bam = str(rules._nanomethphase.output.HP1_bis_bam),
        HP2_bam = str(rules._nanomethphase.output.HP2_bis_bam),
        HP1_bai = str(rules._nanomethphase.output.HP1_bis_bai),
        HP2_bai = str(rules._nanomethphase.output.HP2_bis_bai),
        HP1_freq = str(rules._nanomethphase.output.HP1_freq),
        HP2_freq = str(rules._nanomethphase.output.HP2_freq)

    output:
        HP1_bam = CFG["dirs"]["outputs"] + "bam2bis/{seq_type}--{genome_build}/{sample_id}_NanoMethPhase_HP1_Converted2Bisulfite.bam",
        HP2_bam = CFG["dirs"]["outputs"] + "bam2bis/{seq_type}--{genome_build}/{sample_id}_NanoMethPhase_HP2_Converted2Bisulfite.bam",
        HP1_bai = CFG["dirs"]["outputs"] + "bam2bis/{seq_type}--{genome_build}/{sample_id}_NanoMethPhase_HP1_Converted2Bisulfite.bam.bai",
        HP2_bai = CFG["dirs"]["outputs"] + "bam2bis/{seq_type}--{genome_build}/{sample_id}_NanoMethPhase_HP2_Converted2Bisulfite.bam.bai",
        HP1_freq = CFG["dirs"]["outputs"] + "meth_freq/{seq_type}--{genome_build}/{sample_id}_NanoMethPhase_HP1_MethylFrequency.tsv",
        HP2_freq = CFG["dirs"]["outputs"] + "meth_freq/{seq_type}--{genome_build}/{sample_id}_NanoMethPhase_HP2_MethylFrequency.tsv"

    run:
        op.relative_symlink(input.HP1_bam, output.HP1_bam, in_module= True),
        op.relative_symlink(input.HP2_bam, output.HP2_bam, in_module= True),
        op.relative_symlink(input.HP1_bai, output.HP1_bai, in_module= True),
        op.relative_symlink(input.HP2_bai, output.HP2_bai, in_module= True),
        op.relative_symlink(input.HP1_freq, output.HP1_freq, in_module= True),
        op.relative_symlink(input.HP2_freq, output.HP2_freq, in_module= True)

# Generates the target sentinels for each run, which generate the symlinks
rule _nanomethphase_all:
    input:
        expand(
            [
                str(rules._nanomethphase_output.output.HP1_bam),
                str(rules._nanomethphase_output.output.HP2_bam),
                str(rules._nanomethphase_output.output.HP1_bai),
                str(rules._nanomethphase_output.output.HP2_bai),
                str(rules._nanomethphase_output.output.HP1_freq),
                str(rules._nanomethphase_output.output.HP2_freq)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)  

