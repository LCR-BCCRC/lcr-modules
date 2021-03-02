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
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
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
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "flagstat", "pathseq", "outputs"],
)

# Filter samples to include only tumours and discard normal samples
config["lcr-modules"]["pathseq"]["samples"] = op.discard_samples(SAMPLES, tissue_status="normal")


# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _pathseq_input_bam,
    _pathseq_input_references,
    _pathseq_reference_img,
    _pathseq_reference_bfi,
    _pathseq_download_microbe_references,
    _pathseq_download_taxonomy,
    _pathseq_collect_flagstats,
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


# Setup shared reference files. Symlinking these files to 00-inputs to ensure index and dictionary are present
# before pipeline starts, otherwise on a fresh directory dictionary is not created and worflow exits.
rule _pathseq_input_references: 
    input: 
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        genome_fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        genome_dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output: 
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa", 
        genome_fai = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa.fai", 
        genome_dict = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.dict"
    shell: 
        op.as_one_line("""
        ln -s {input.genome_fa} {output.genome_fa} &&
        ln -s {input.genome_fai} {output.genome_fai} &&
        ln -s {input.genome_dict} {output.genome_dict}
        """)


# Create the k-mer file following GATK recommendations
rule _pathseq_reference_img:
    input: 
        genome_fa = str(rules._pathseq_input_references.output.genome_fa),
        genome_fai = str(rules._pathseq_input_references.output.genome_fai)
    output:
        genome_img = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.fa.img"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["reference"]
    resources:
        **CFG["resources"]["reference"]
    shell:
        op.as_one_line("""
          gatk-launch PathSeqBuildKmers
          --referencePath {input.genome_fa}
          --bloomFalsePositiveProbability 0.001
          --kSize 31
          --kmerMask 15
          -O {output.genome_img}
        """)    


# Create BWA-MEM index image of the reference
rule _pathseq_reference_bfi:
    input: 
        genome_fa = str(rules._pathseq_input_references.output.genome_fa),
        genome_fai = str(rules._pathseq_input_references.output.genome_fai)
    output:
        genome_bfi = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_{genome_build}.bfi"
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
        """)


# Download specific reference files
rule _pathseq_download_microbe_references: 
    input: 
        genome_fa = str(rules._pathseq_input_references.output.genome_fa)
    output: 
        microbe_fasta = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_microbe.{genome_build}.fa",
        microbe_fai = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_microbe.{genome_build}.fa.fai",
        microbe_dict = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_microbe.{genome_build}.dict",
        microbe_image = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_microbe.{genome_build}.fa.img"
    params:
        url = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_microbe.tar.gz',
        folder = CFG["dirs"]["inputs"] + "references/microbe_reference"
    shell: 
        op.as_one_line("""
        wget -qO- {params.url} |
        tar xzf - -C {params.folder}
          &&
        ln -s {params.folder}/pathseq_microbe.fa {output.microbe_fasta}
          &&
        ln -s {params.folder}/pathseq_microbe.fa.fai {output.microbe_fai}
          &&
        ln -s {params.folder}/pathseq_microbe.dict {output.microbe_dict}
          &&
        ln -s {params.folder}/pathseq_microbe.fa.img {output.microbe_image}       
        """)


# Download pathogene taxonomy
rule _pathseq_download_taxonomy: 
    input: 
        genome_fa = str(rules._pathseq_input_references.output.genome_fa)
    output:
        taxonomy = CFG["dirs"]["inputs"] + "references/microbe_reference/pathseq_taxonomy.{genome_build}.db"
    params:
        url = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_taxonomy.tar.gz',
        folder = CFG["dirs"]["inputs"] + "references/microbe_reference"
    shell: 
        op.as_one_line("""
        wget -qO- {params.url} |
        tar xzf - -C {params.folder}
          &&
        ln -s {params.folder}/pathseq_taxonomy.db {output.taxonomy}
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


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _pathseq_run:
    input:
        bam = str(rules._pathseq_input_bam.output.bam),
        genome_img = str(rules._pathseq_reference_img.output.genome_img),
        k_mers = str(rules._pathseq_reference_bfi.output.genome_bfi),
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
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["pathseq"]
    resources:
        **CFG["resources"]["pathseq"]
    shell:
        op.as_one_line("""
          gatk PathSeqPipelineSpark --spark-master local[{threads}]
          --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
          --input {input.bam}
          --filter-bwa-image {input.genome_img}
          --kmer-file {input.k_mers}
          --min-clipped-read-length {params.min_read_length} 
          --microbe-dict {input.microbe_dict}
          --microbe-bwa-image {input.microbe_image}
          --taxonomy-file {input.taxonomy}
          --output {output.bam}
          --scores-output {output.scores}
          > {log.stdout}
          2> {log.stderr}
        """)


# Calculate the proportion of EBV reads in the given sample
rule _pathseq_calculate_ebv:
    input:
        scores = str(rules._pathseq_run.output.scores),
        flagstats = str(rules._pathseq_collect_flagstats.output.flagstat)
    output:
        ebv_status = CFG["dirs"]["pathseq"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.{genome_build}.pathseq_ebv_results.bam"
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
        scores = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{sample_id}.pathseq_scores.txt",
        ebv_status = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{sample_id}.ebv_status.txt"
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
