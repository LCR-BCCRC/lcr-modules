#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
# Contributors:     N/A

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


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["pathseq"]`
CFG = op.setup_module(
    name = "pathseq",
    version = "1.0",
    subdirectories = ["inputs", "flagstat", "pathseq", "outputs"],
)


# Define rules to be run locally when using a compute cluster
localrules:
    _pathseq_input_bam,
    _pathseq_input_fasta,
    _pathseq_fasta_dictionary,
    _pathseq_download_microbe_references,
    _pathseq_download_taxonomy,
    _pathseq_calculate_ebv,
    _pathseq_output,
    _pathseq_all

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _pathseq_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bam+ ".bai", output.bai)
        op.absolute_symlink(input.bam+ ".bai", output.crai)


# Setup shared reference files. We need fasta file without EBV chromosome to be used in Pathseq module, so will process
# fasta file here and create index for the new fasta
rule _pathseq_input_fasta: 
    input: 
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: 
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.noEBV.fa", 
        genome_fai = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.noEBV.fa.fai"
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    shell: 
        op.as_one_line("""
        cat {input.genome_fa} | perl -ne 'print if not /^\>(chr)EBV/../^[X]+\s.+/' > {output.genome_fa} &&
        samtools faidx {output.genome_fa} -o {output.genome_fai}
        """)

# create dictionary for the fasta without EBV
rule _pathseq_fasta_dictionary:
    input: 
        genome_fa = str(rules._pathseq_input_fasta.output.genome_fa),
        genome_fai = str(rules._pathseq_input_fasta.output.genome_fai)
    output:
        genome_dict = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.noEBV.fa.dict"
    log:
        stdout = CFG["logs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.CreateSequenceDictionary.stdout.log",
        stderr = CFG["logs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.CreateSequenceDictionary.stderr.log"
    params:
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    shell:
        op.as_one_line("""
          picard CreateSequenceDictionary
          R={input.genome_fa}
          O={output.genome_dict}
          >> {log.stdout}
          2>> {log.stderr}
        """)    

# Create k-mer file following GATK recommendations
rule _pathseq_reference_img:
    input: 
        genome_fa = str(rules._pathseq_input_fasta.output.genome_fa),
        genome_fai = str(rules._pathseq_input_fasta.output.genome_fai),
        genome_dict= str(rules._pathseq_fasta_dictionary.output.genome_dict)
    output:
        genome_img = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.hss"
    log:
        stdout = CFG["logs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.PathSeqBuildKmers.stdout.log",
        stderr = CFG["logs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.PathSeqBuildKmers.stderr.log"
    params:
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    shell:
        op.as_one_line("""
          gatk PathSeqBuildKmers
          --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
          --reference {input.genome_fa}
          --output {output.genome_img}
          >> {log.stdout}
          2>> {log.stderr}
        """)    


# Create BWA-MEM index image of the reference
rule _pathseq_reference_bfi:
    input: 
        genome_fa = str(rules._pathseq_input_fasta.output.genome_fa),
        genome_fai = str(rules._pathseq_input_fasta.output.genome_fai),
        genome_dict= str(rules._pathseq_fasta_dictionary.output.genome_dict)
    output:
        genome_bfi = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.img"
    log:
        stdout = CFG["logs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.BwaMemIndexImageCreator.stdout.log",
        stderr = CFG["logs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.BwaMemIndexImageCreator.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    shell:
        op.as_one_line("""
          gatk BwaMemIndexImageCreator
          -I {input.genome_fa}
          -O {output.genome_bfi}
          >> {log.stdout}
          2>> {log.stderr}
        """)


# Download specific reference files
rule _pathseq_download_microbe_references: 
    output: 
        microbe_fasta = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_microbe.fa",
        microbe_fai = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_microbe.fa.fai",
        microbe_dict = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_microbe.dict",
        microbe_image = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_microbe.fa.img"
    params:
        url = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_microbe.tar.gz',
        folder = CFG["dirs"]["inputs"] + "references/microbe_reference"
    shell: 
        op.as_one_line("""
        wget -qO- {params.url} |
        tar xzf - -C {params.folder}  
        """)


# Download pathogene taxonomy
rule _pathseq_download_taxonomy: 
    output:
        taxonomy = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_taxonomy.db"
    params:
        url = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_taxonomy.tar.gz',
        folder = CFG["dirs"]["inputs"] + "references/microbe_reference"
    shell: 
        op.as_one_line("""
        wget -qO- {params.url} |
        tar xzf - -C {params.folder}
        """)


# Collect flagstat information for the read normalization
rule _pathseq_collect_flagstats:
    input:
        bam = str(rules._pathseq_input_bam.output.bam),
    output:
        flagstat = CFG["dirs"]["flagstat"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.{genome_build}.flagstat"
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    shell:
        op.as_one_line("""
        samtools flagstat {input.bam} > {output.flagstat}
        """)


# Run Pathseq analysis
rule _pathseq_run:
    input:
        bam = str(rules._pathseq_input_bam.output.bam),
        k_mers = str(rules._pathseq_reference_img.output.genome_img),
        genome_img = str(rules._pathseq_reference_bfi.output.genome_bfi),
        taxonomy = str(rules._pathseq_download_taxonomy.output.taxonomy),
        microbe_dict = str(rules._pathseq_download_microbe_references.output.microbe_dict),
        microbe_image = str(rules._pathseq_download_microbe_references.output.microbe_image)
    output:
        bam = temp(CFG["dirs"]["pathseq"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.{genome_build}.pathseq.bam"),
        scores = CFG["dirs"]["pathseq"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.{genome_build}.pathseq.scores.txt"
    log:
        stdout = CFG["logs"]["pathseq"] + "{seq_type}--{genome_build}/{sample_id}/run_pathseq.stdout.log",
        stderr = CFG["logs"]["pathseq"] + "{seq_type}--{genome_build}/{sample_id}/run_pathseq.stderr.log"
    params:
        min_read_length = CFG["options"]["min_read_length"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        flags = CFG["options"]["flags"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["pathseq"]
    resources:
        **CFG["resources"]["pathseq"]
    shell:
        op.as_one_line("""
          if [ -e {output.bam}.parts ]; then rm -R {output.bam}.parts; fi
              &&
          echo "running {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
          gatk PathSeqPipelineSpark --spark-master local[{threads}]
          --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
          {params.flags}
          --input {input.bam}
          --filter-bwa-image {input.genome_img}
          --kmer-file {input.k_mers}
          --min-clipped-read-length {params.min_read_length} 
          --microbe-dict {input.microbe_dict}
          --microbe-bwa-image {input.microbe_image}
          --taxonomy-file {input.taxonomy}
          --output {output.bam}
          --scores-output {output.scores}
          >> {log.stdout}
          2>> {log.stderr} &&  
          echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)


# Calculate the proportion of EBV reads in the given sample
rule _pathseq_calculate_ebv:
    input:
        scores = str(rules._pathseq_run.output.scores),
        flagstats = str(rules._pathseq_collect_flagstats.output.flagstat)
    output:
        ebv_status = CFG["dirs"]["pathseq"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.{genome_build}.pathseq_ebv_results.tsv"
    params:
        opts = CFG["options"]["ebv_cutoff"]
    conda:
        CFG["conda_envs"]["R"]
    threads:
        CFG["threads"]["calculate_ebv"]
    resources:
        **CFG["resources"]["calculate_ebv"]
    script:
        "src/R/calculate_ebv.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _pathseq_output:
    input:
        bam = str(rules._pathseq_run.output.bam),
        scores = str(rules._pathseq_run.output.scores),
        ebv_status = str(rules._pathseq_calculate_ebv.output.ebv_status)
    output:
        scores = CFG["dirs"]["outputs"] + "scores/{seq_type}--{genome_build}/{sample_id}.pathseq_scores.txt",
        ebv_status = CFG["dirs"]["outputs"] + "ebv_status/{seq_type}--{genome_build}/{sample_id}.ebv_status.tsv"
    run:
        op.relative_symlink(input.scores, output.scores, in_module=True)
        op.relative_symlink(input.ebv_status, output.ebv_status, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _pathseq_all:
    input:
        expand(
            [
                str(rules._pathseq_output.output.scores),
                str(rules._pathseq_output.output.ebv_status)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
