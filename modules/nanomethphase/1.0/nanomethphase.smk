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
    subdirectories = ["inputs", "methyl_call", "nanomethphase", "dma", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _nanomethphase_input_promethion,
    _nanomethphase_input_methylation, 
    _nanomethphase_input_vcf,
    _nanomethphase_output_bam,
    _nanomethphase_output_dma,
    _nanomethphase_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _nanomethphase_input_promethion:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        
rule _nanomethphase_input_methylation: 
    input:
        mc = CFG["inputs"]["meth_calls"]
    output:
        mc = CFG["dirs"]["inputs"] + "meth_calls/{seq_type}--{genome_build}/{sample_id}.calls.tsv.gz"
    run:
        op.absolute_symlink(input.mc, output.mc)

rule _nanomethphase_input_vcf: 
    input:
        vcf = CFG["inputs"]["vcf"],
        index = CFG["inputs"]["index"]
    output:
        vcf = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.vcf.gz",
        index = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.vcf.gz.tbi"
    run:
        op.absolute_symlink(input.vcf, output.vcf)
        op.absolute_symlink(input.index, output.index)

rule _nanomethphase_methyl_call:
    input:
        calls = CFG["dirs"]["inputs"] + "meth_calls/{seq_type}--{genome_build}/{sample_id}.calls.tsv.gz"
    output: 
        bed = CFG["dirs"]["methyl_call"] + "{seq_type}--{genome_build}/{sample_id}.bed.gz"
    threads:
        CFG["threads"]["methyl_call"]
    conda:
        CFG["conda_envs"]["nanomethphase"] 
    resources: **CFG["resources"]["methyl_call"]        
    log: CFG["logs"]["methyl_call"] + "{seq_type}--{genome_build}/{sample_id}/methyl_call_processor.stderr.log"             
    shell:
        op.as_one_line("""
        nanomethphase methyl_call_processor -mc {input.calls} --threads {threads} | 
        sort -k1,1 -k2,2n -k3,3n | bgzip > {output.bed} 2> {log} && tabix -p bed {output.bed}
        """)   

rule _nanomethphase_run:
    input:  
        bam = str(rules._nanomethphase_input_promethion.output.bam),
        bai = str(rules._nanomethphase_input_promethion.output.bai),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        vcf = str(rules._nanomethphase_input_vcf.output.vcf),
        tbi = str(rules._nanomethphase_input_vcf.output.index),
        mc = str(rules._nanomethphase_methyl_call.output.bed), 
        realbam =CFG["inputs"]["sample_bam"], # Enables this module to make use of temp bam files
        realbai = CFG["inputs"]["sample_bai"]
    threads:
        CFG["threads"]["nanomethphase"]
    resources: **CFG["resources"]["nanomethphase"]
    params:
        prefix =  CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}",
        options = CFG["options"]["nanomethphase"]
    conda:
        CFG["conda_envs"]["nanomethphase"] 
    log: CFG["logs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/nanomethphase.log"         
    output:
        HP1_bis_bam = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP1_Converted2Bisulfite.bam",
        HP2_bis_bam = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP2_Converted2Bisulfite.bam",
        HP1_freq = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP1_MethylFrequency.tsv",
        HP2_freq = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP2_MethylFrequency.tsv"
    shell:
        op.as_one_line("""
        nanomethphase phase -mc {input.mc} -of methylcall,bam2bis -o {params.prefix} 
        -b {input.bam} -r {input.fasta} -v {input.vcf} --overwrite {params.options} 
        -t {threads} -is > {log} 2>&1 
        """)  
        
rule _nanomethphase_cram: 
    input: 
        bam = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP{haplotype}_Converted2Bisulfite.bam",
        bai = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP{haplotype}_Converted2Bisulfite.bam.bai", 
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        cram = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP{haplotype}_Converted2Bisulfite.cram",
        crai = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP{haplotype}_Converted2Bisulfite.cram.crai"
    threads: CFG["threads"]["cram"]
    resources: **CFG["resources"]["cram"]
    conda: 
        CFG["conda_envs"]["nanomethphase"]
    log: CFG["logs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/HP{haplotype}_cram_conversion.log"
    shell:
        op.as_one_line("""
            samtools view -@ {threads} -T {input.fasta} -C -o {output.cram} {input.bam} > {log} 2>&1 &&
            samtools index -@ {threads} {output.cram} >> {log} 2>&1
        """)

rule _nanomethphase_dma: 
    input: 
        HP1_freq = str(rules._nanomethphase_run.output.HP1_freq),
        HP2_freq = str(rules._nanomethphase_run.output.HP2_freq)
    output:
        dml = CFG["dirs"]["dma"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_callDML.txt", 
        dmr = CFG["dirs"]["dma"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_callDMR.txt"
    params: 
        outdir = CFG["dirs"]["dma"] + "{seq_type}--{genome_build}/{sample_id}/", 
        options = CFG["options"]["dma"]
    conda:
        CFG["conda_envs"]["nanomethphase"] 
    threads: CFG["threads"]["dma"]
    resources: **CFG["resources"]["dma"]
    shell: 
        op.as_one_line("""
            nanomethphase dma -c 1,2,4,5,7 
                -rs "$(which Rscript) --vanilla"
                -ca {input.HP1_freq}
                -co {input.HP2_freq}
                -o {params.outdir}
                -op {wildcards.sample_id}
                --overwrite 
                {params.options}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _nanomethphase_output_cram:
    input:
        cram = str(rules._nanomethphase_cram.output.cram),
        crai = str(rules._nanomethphase_cram.output.crai),
        freq = CFG["dirs"]["nanomethphase"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_NanoMethPhase_HP{haplotype}_MethylFrequency.tsv"
    output:
        cram = CFG["dirs"]["outputs"] + "bam2bis/{seq_type}--{genome_build}/{sample_id}_NanoMethPhase_HP{haplotype}_Converted2Bisulfite.cram",
        crai = CFG["dirs"]["outputs"] + "bam2bis/{seq_type}--{genome_build}/{sample_id}_NanoMethPhase_HP{haplotype}_Converted2Bisulfite.cram.crai",
        freq = CFG["dirs"]["outputs"] + "meth_freq/{seq_type}--{genome_build}/{sample_id}_NanoMethPhase_HP{haplotype}_MethylFrequency.tsv"
    run:
        op.relative_symlink(input.cram, output.cram, in_module= True),
        op.relative_symlink(input.crai, output.crai, in_module= True),
        op.relative_symlink(input.freq, output.freq, in_module= True)

rule _nanomethphase_output_dma: 
    input: 
        dml = str(rules._nanomethphase_dma.output.dml),
        dmr = str(rules._nanomethphase_dma.output.dmr)
    output: 
        dml = CFG["dirs"]["outputs"] + "meth_freq/{seq_type}--{genome_build}/{sample_id}_callDML.txt",
        dmr = CFG["dirs"]["outputs"] + "meth_freq/{seq_type}--{genome_build}/{sample_id}_callDMR.txt"
    run:
        op.relative_symlink(input.dml, output.dml, in_module= True),
        op.relative_symlink(input.dmr, output.dmr, in_module= True)        

# Generates the target sentinels for each run, which generate the symlinks
rule _nanomethphase_all:
    input:
        expand(
            expand(
                [
                    str(rules._nanomethphase_output_cram.output.cram),
                    str(rules._nanomethphase_output_cram.output.crai),
                    str(rules._nanomethphase_output_cram.output.freq),
                    str(rules._nanomethphase_output_dma.output.dml),
                    str(rules._nanomethphase_output_dma.output.dmr)
                ],
                haplotype=[1, 2], 
                allow_missing=True
            ),
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)  

