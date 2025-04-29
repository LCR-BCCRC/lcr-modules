#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Jasper
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import glob
import os
import re
import sys

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
# `CFG` is a shortcut to `config["lcr-modules"]["purecn"]`
CFG = op.setup_module(
    name = "purecn",
    version = "1.0",
    subdirectories = ["inputs", "normals", "pon", "coverage", "mutect2", "pureCN_cnvkit", "pureCN", "convert_coordinates", "fill_regions", "normalize", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _purecn_input_bam,
    _purecn_symlink_cnvkit_seg,
    _purecn_symlink_cnvkit_cnr,
    _purecn_cnvkit_output_seg,
    _purecn_denovo_output_seg,
    _purecn_cnvkit_output_files,
    _purecn_denovo_output_files,
    coverage_list_normals,
    _purecn_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')

# Input only accepts bam files
rule _purecn_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"],
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Import CNVkit seg and coverage information as needed
if CFG["cnvkit_seg"] == True:
    rule _purecn_symlink_cnvkit_seg:
        input:
            seg = CFG["cnvkit_directory"] + "09-seg/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg"
        output:
            seg = CFG["dirs"]["inputs"] + "seg/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg"
        run:
            op.absolute_symlink(input.seg, output.seg)

if CFG["cnvkit_seg"] == True:
    rule _purecn_symlink_cnvkit_cnr:
        input:
            cnr = CFG["cnvkit_directory"] + "02-fix/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr"
        output:
            cnr = CFG["dirs"]["inputs"] + "cnr/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr"
        run:
            op.absolute_symlink(input.cnr, output.cnr)


# -------------------------------------------------------------------------------------------------- #
# Part I - Set up reference files
# -------------------------------------------------------------------------------------------------- #
rule _download_GEM:
    output:
        touch(CFG["dirs"]["inputs"] + "references/GEM/.done")
    params:
        dirOut = CFG["dirs"]["inputs"] + "references/GEM/"
    conda:
        CFG["conda_envs"]["wget"]
    resources: **CFG["resources"]["gem"]
    shell:
        """
            wget https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download -O {params.dirOut}/GEM-lib.tbz2 && bzip2 -dc {params.dirOut}/GEM-lib.tbz2 | tar -xvf - -C {params.dirOut}/
        """


# Generate mappability files
# grch37 and grch38 from ensembl have additional information in header - need to remove
def which_genome_fasta(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    if  "38" in str({wildcards.genome_build}):
        return reference_files("genomes/grch38_masked/genome_fasta/genome.fa")
    else:
        return reference_files("genomes/grch37_masked/genome_fasta/genome.fa")
    

rule _set_up_grch_genomes:
    input:
        reference = which_genome_fasta
    output:
        reference = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/genome_header.fa"
    resources: **CFG["resources"]["gem"]
    shell:
        "cat {input.reference} | perl -ne 's/(^\>\S+).+/$1/;print;' > {output.reference} "


def get_genome_fasta(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    if  "grch" in str({wildcards.genome_build}):
        if "38" in str({wildcards.genome_build}):
            return  CFG["dirs"]["inputs"] + "references/grch38_masked/freec/genome_header.fa"
        elif "37" in str({wildcards.genome_build}):
            return CFG["dirs"]["inputs"] + "references/grch37_masked/freec/genome_header.fa"
        else:
            raise AttributeError(f"The specified grch genome build is not available for use.")
    elif "hs37d5" in str({wildcards.genome_build}):
        return CFG["dirs"]["inputs"] + "references/grch37_masked/freec/genome_header.fa"
    else:
        if "38" in str({wildcards.genome_build}):
            return reference_files("genomes/hg38_masked/genome_fasta/genome.fa")
        elif "19" in str({wildcards.genome_build}):
            return reference_files("genomes/hg19_masked/genome_fasta/genome.fa")
        else:
            raise AttributeError(f"The specified genome build is not available for use.")
        
        
rule _purecn_generate_gem_index:
    input:
        software = CFG["dirs"]["inputs"] + "references/GEM/.done",
        reference = get_genome_fasta,
        idxpref = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all_index"
    output:
        index = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all_index.gem"
    params:
        gemDir = CFG["dirs"]["inputs"] + "references/GEM/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin"
    threads: CFG["threads"]["gem"]
    resources: **CFG["resources"]["gem"]
    log: CFG["logs"]["inputs"] + "gem/{genome_build}/gem_index.stderr.log"
    shell:
        "PATH=$PATH:{params.gemDir}; {params.gemDir}/gem-indexer -T {threads} -c dna -i {input.reference} -o {input.idxpref} > {log} 2>&1 "


rule _purecn_generate_mappability:
    input:
        software = CFG["dirs"]["inputs"] + "references/GEM/.done",
        index = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all_index.gem",
        pref =  CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem"
    output:
        mappability = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem.mappability"
    params:
        gemDir = CFG["dirs"]["inputs"] + "references/GEM/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin",
        kmer = CFG["options"]["kmer"],
        mismatch = CFG["options"]["mismatch"],
        maxEditDistance = CFG["options"]["maxEditDistance"],
        maxBigIndel = CFG["options"]["maxBigIndel"],
        strata = CFG["options"]["strata"]
    threads: CFG["threads"]["gem"]
    resources: **CFG["resources"]["gem"]
    log: CFG["logs"]["inputs"] + "gem/{genome_build}/gem_map.stderr.log"
    shell:
        "PATH=$PATH:{params.gemDir}; {params.gemDir}/gem-mappability -T {threads} -I {input.index} -l {params.kmer} -m {params.mismatch} -t disable --mismatch-alphabet ACGNT -e {params.maxEditDistance} --max-big-indel-length {params.maxBigIndel} -s {params.strata} -o {input.pref} > {log} 2>&1 "


rule _purecn_symlink_map:
    input:
        mappability = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem.mappability"
    output:
        mappability = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/out100m2_{genome_build}.gem"
    resources: **CFG["resources"]["gem"]
    run:
        op.relative_symlink(input.mappability, output.mappability)
            
            
rule _purecn_set_mappability:
    input:
        index = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all_index.gem",
        mappability = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/out100m2_{genome_build}.gem",
        pref =  CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem"
    output:
        wig = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem.wig",
        size = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem.sizes"
    params:
        gemDir = CFG["dirs"]["inputs"] + "references/GEM/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin"
    threads: CFG["threads"]["gem"]
    resources: **CFG["resources"]["gem"]
    shell:
        """
            PATH=$PATH:{params.gemDir}; 
            {params.gemDir}/gem-2-wig -I {input.index} -i {input.mappability} -o {input.pref} 
        """

rule _purecn_gem_wig2bw:
    input:
        wig = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem.wig",
        sizes = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem.sizes"
    output:
        bw = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem.bw",
    conda: CFG["conda_envs"]["ucsc_bigwigtowig"]
    threads: CFG["threads"]["gem"]
    resources: **CFG["resources"]["gem"]
    shell:
        """
            wigToBigWig {input.wig} {input.sizes} {output.bw}
        """

# set up intervals for coverage calculations (used for de novo pureCN CNV calling)
rule _purecn_setup_blacklist_bed:
    input:
        blacklist = reference_files("genomes/{genome_build}/repeatmasker/repeatmasker.{genome_build}.bed")
    output:
        blacklist = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/repeatmasker.{genome_build}.bed"
    shell:
        """
            awk '{{print $1"\t"$2"\t"$3}}' {input.blacklist} > {output.blacklist}
        """

def _purecn_get_capspace(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    # default is just taking values from genome_build and capture_space in metadata
    default_bed = reference_files("genomes/" + wildcards.genome_build + "/capture_space/" + wildcards.capture_space + ".bed")
    # if there is a custom bed in the config that pertains to the capture_space variable, then use that
    if str(wildcards.capture_space) in CFG["options"]["target_bed"][wildcards.genome_build].keys():
        bed = CFG["options"]["target_bed"][wildcards.genome_build][wildcards.capture_space]
    # if there is "none" instead, then we go towards the default, which is the full exon-space
    elif "none" in str({wildcards.capture_space}):
        if "38" in str({wildcards.genome_build}):
            bed = reference_files("genomes/" + wildcards.genome_build + "/capture_space/exome-utr-grch38.bed")
        else:
            bed = reference_files("genomes/" + wildcards.genome_build + "/capture_space/exome-utr-grch37.bed")
    # if it's anything else, then we can assume it will just use the default, which is what can be taken from the metadata and found via reference_files()
    else:
        bed = default_bed
    return bed

rule _purecn_input_bed: 
    input:
        bed = _purecn_get_capspace
    output:
        bed = CFG["dirs"]["inputs"] + "references/{genome_build}/{capture_space}/baits_{genome_build}.bed"
    run:
        op.absolute_symlink(input.bed, output.bed)

rule _purecn_setinterval:
    input:
        genome = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        bed = str(rules._purecn_input_bed.output.bed),
        bw = CFG["dirs"]["inputs"] + "references/{genome_build}_masked/freec/{genome_build}.hardmask.all.gem.bw",
    output:
        intervals = CFG["dirs"]["inputs"] + "references/{genome_build}/{capture_space}/baits_{genome_build}_intervals.txt"
    params:
        genome_build = lambda w: config["lcr-modules"]["purecn"]["options"]["genome_builds"][w.genome_build],
        intervalfile_script = CFG["software"]["intervalfile_script"],
        outdir =  CFG["dirs"]["inputs"] + "references/{genome_build}/{capture_space}/",
        force = CFG["options"]["setinterval"]["force"],
        opts = CFG["options"]["setinterval"]["opts"]
    conda: CFG["conda_envs"]["purecn"]
    resources:
        **CFG["resources"]["setinterval"]
    threads: CFG["threads"]["setinterval"]
    log:
        CFG["logs"]["inputs"] + "purecn_setinterval/{genome_build}--{capture_space}.log"
    shell:
        """
            PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
            echo -e "Using {params.intervalfile_script} instead of default $PURECN/IntervalFile.R..."
            mkdir -p {params.outdir} ;
            Rscript --vanilla {params.intervalfile_script} --in-file {input.bed} \
            --fasta {input.genome} --out-file {output.intervals} \
            --genome {params.genome_build} \
            --mappability {input.bw} {params.force} {params.opts} > {log} 2>&1
        """

# -------------------------------------------------------------------------------------------------- #
# Part II - Set up normals 
# -------------------------------------------------------------------------------------------------- #

# Run Mutect2
# Germline (normal) samples
rule _purecn_mutect2_germline:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        dbsnp = ancient(reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")),
        fasta = ancient(reference_files("genomes/{genome_build}/genome_fasta/genome.fa")),
        gnomad = ancient(reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"))
    output:
        vcf = temp(CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}.{chrom}.vcf.gz.tbi"),
        stats = temp(CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}.{chrom}.vcf.gz.stats"),
        f1r2 = temp(CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}.{chrom}.f1r2.tar.gz"),
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
        opts = CFG["options"]["mutect2_norm"]["mutect2_opts"]
    log: CFG["logs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{chrom}.log"
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    shell:
        """
            gatk Mutect2 --java-options "-Xmx{params.mem_mb}m" {params.opts} --genotype-germline-sites true --genotype-pon-sites true --interval-padding 50 --max-mnp-distance 0 --germline-resource {input.gnomad} -R {input.fasta} -L {wildcards.chrom} -I {input.bam} -O {output.vcf} --f1r2-tar-gz {output.f1r2} > {log} 2>&1
        """

#### set-up mpileups for BAF calling ####
def _mutect2_normal_get_chr_vcf(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    mpileups = expand(
        CFG["dirs"]["normals"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{normal_id}}/{{normal_id}}.{chrom}.vcf.gz", 
        chrom = chrs
    )
    return(mpileups)
    
    
def _mutect2_normal_get_chr_vcf_tbi(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    mpileups = expand(
        CFG["dirs"]["normals"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{normal_id}}/{{normal_id}}.{chrom}.vcf.gz.tbi", 
        chrom = chrs
    )
    return(mpileups)
    

rule _purecn_mutect2_concatenate_vcf:
    input: 
        vcf = _mutect2_normal_get_chr_vcf,
        tbi = _mutect2_normal_get_chr_vcf_tbi,
    output: 
        vcf = temp(CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}.vcf.gz"),
        tbi = temp(CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}.vcf.gz.tbi")
    resources: 
        **CFG["resources"]["concatenate_vcf"]
    conda:
        CFG["conda_envs"]["bcftools"]
    shell: 
        """
            bcftools concat {input.vcf} -Oz -o {output.vcf} && 
            tabix -p vcf {output.vcf}
        """

def _purecn_mutect2_normal_get_chr_stats(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    stats = expand(
        CFG["dirs"]["normals"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{normal_id}}/{{normal_id}}.{chrom}.vcf.gz.stats", 
        chrom = chrs
    )
    return(stats)


# Merge chromosome mutect2 stats for FilterMutectCalls rule
rule _purecn_mutect2_normal_merge_stats:
    input:
        stats = _purecn_mutect2_normal_get_chr_stats
    output:
        stats = CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}.vcf.gz.stats"
    log:
        CFG["logs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/stats.log"
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    shell:
        op.as_one_line("""
        gatk MergeMutectStats $(for i in {input.stats}; do echo -n "-stats $i "; done)
        -O {output.stats} > {log} 2>&1
        """)

# Get pileup summaries
rule _purecn_mutect2_normal_pileup_summaries: 
    input: 
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
        snps = reference_files("genomes/{genome_build}/gatk/mutect2_small_exac.{genome_build}.vcf.gz"), 
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: 
        pileup = CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/pileupSummary.table"
    log: 
        CFG["logs"]["normals"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{normal_id}/pileupSummary.log"
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    threads: CFG["threads"]["mutect2"]
    shell: 
        op.as_one_line("""
        gatk GetPileupSummaries 
            --java-options "-Xmx{params.mem_mb}m"
            -I {input.normal_bam}
            -R {input.fasta} 
            -V {input.snps}
            -L {input.snps}
            -O {output.pileup}
            > {log} 2>&1
        """)

# Calculate contamination  
rule _purecn_mutect2_normal_calc_contamination: 
    input: 
        pileup = str(rules._purecn_mutect2_normal_pileup_summaries.output.pileup)
    output: 
        segments =  CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/segments.table", 
        contamination =  CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/contamination.table"
    log: 
        CFG["logs"]["normals"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{normal_id}/contam.log"
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    threads: CFG["threads"]["mutect2"]
    shell: 
        op.as_one_line("""
        gatk CalculateContamination 
            --java-options "-Xmx{params.mem_mb}m"
            -I {input.pileup}
            -tumor-segmentation {output.segments}
            -O {output.contamination}
            > {log} 2>&1
        """)
    
# Learn read orientation model
def _purecn_mutect2_normal_get_chr_f1r2(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    f1r2 = expand(
        CFG["dirs"]["normals"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{normal_id}}/{{normal_id}}.{chrom}.f1r2.tar.gz", 
        chrom = chrs
    )
    return(f1r2)

rule _purecn_mutect2_normal_learn_orient_model: 
    input: 
        f1r2 = _purecn_mutect2_normal_get_chr_f1r2
    output:
        model =  CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/read-orientation-model.tar.gz"
    log: 
        stdout = CFG["logs"]["normals"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{normal_id}/read-orientation-model.stdout.log",
        stderr = CFG["logs"]["normals"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{normal_id}/read-orientation-model.stderr.log"
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    threads: CFG["threads"]["mutect2"]
    shell: 
        op.as_one_line("""
        inputs=$(for input in {input.f1r2}; do printf -- "-I $input "; done);
        gatk LearnReadOrientationModel 
        --java-options "-Xmx{params.mem_mb}m" 
        $inputs -O {output.model}
        > {log.stdout} 2> {log.stderr}
        """)
    

# Marks variants filtered or PASS annotations
rule _purecn_mutect2_normal_filter:
    input:
        vcf = str(rules._purecn_mutect2_concatenate_vcf.output.vcf),
        tbi = str(rules._purecn_mutect2_concatenate_vcf.output.tbi),
        stat = str(rules._purecn_mutect2_normal_merge_stats.output.stats),
        segments = str(rules._purecn_mutect2_normal_calc_contamination.output.segments), 
        contamination = str(rules._purecn_mutect2_normal_calc_contamination.output.contamination), 
        model = str(rules._purecn_mutect2_normal_learn_orient_model.output.model),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = temp(CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/output.unfilt.vcf.gz")
    log:
        CFG["logs"]["normals"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{normal_id}/mutect2_filter.log",
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["mutect2_norm"]["mutect2_filter"]
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    shell:
        op.as_one_line("""
        gatk FilterMutectCalls --java-options "-Xmx{params.mem_mb}m" 
            {params.opts} 
            -V {input.vcf} 
            -R {input.fasta}
            --tumor-segmentation {input.segments}
            --contamination-table {input.contamination}
            --ob-priors {input.model}
            -O {output.vcf} 
            > {log} 2>&1
        """)


# Filters for PASS and germline variants
# This will only take somatic ones if filtering just for PASSED (need to still maintain germline ones)
rule _purecn_mutect2_normal_filter_passed:
    input:
        vcf = str(rules._purecn_mutect2_normal_filter.output.vcf)
    output:
        vcf = CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}_passed.vcf.gz",
        tbi = CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}_passed.vcf.gz.tbi"
    params:
        filter_for_opts = CFG["options"]["mutect2_norm"]["mutect2_filter_for"],
        filter_out_opts = CFG["options"]["mutect2_norm"]["mutect2_filter_out"],
    log:
        stderr = CFG["logs"]["normals"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{normal_id}/mutect2_passed.log"
    conda: CFG["conda_envs"]["bcftools"]
    resources: **CFG["resources"]["concatenate_vcf"]
    shell:
        op.as_one_line(""" 
        bcftools view "{params.filter_for_opts}" -e "{params.filter_out_opts}" -Oz -o {output.vcf} {input.vcf} 2> {log.stderr}
            &&
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


# Creating a panel of normals vcf
# This is a list of file paths to normal vcfs that you want removed from your panel of normal
# You can QC these vcf files after running _purecn_mutect2_germline
omit_normals_list =  CFG["options"]["normals"]["omit_list"]

def _get_normals_vcfs(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    capture_space = CFG["samples"][CFG["samples"]['capture_space'].isin([wildcards.capture_space])]
    capture_space = capture_space[capture_space["genome_build"].isin([wildcards.genome_build])]
    capture_space = capture_space[capture_space["seq_type"].isin([wildcards.seq_type])]
    capture_space = op.filter_samples(capture_space, tissue_status = "normal")
    normals = expand(
        CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}_passed.vcf.gz", 
        zip,
        seq_type = capture_space["seq_type"],
        genome_build = capture_space["genome_build"],
        normal_id = capture_space["sample_id"], 
        capture_space = capture_space["capture_space"])
    normals = list(dict.fromkeys(normals))
    omit_normals_list =  CFG["options"]["normals"]["omit_list"]
    if os.path.exists(omit_normals_list):
        with open(omit_normals_list) as file:
            remove = file.read().rstrip("\n").split("\n")
            normals = [e for e in normals if e not in remove]
    return normals
    
def _get_normals_tbi(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    capture_space = CFG["samples"][CFG["samples"]['capture_space'].isin([wildcards.capture_space])]
    capture_space = capture_space[capture_space["genome_build"].isin([wildcards.genome_build])]
    capture_space = capture_space[capture_space["seq_type"].isin([wildcards.seq_type])]
    capture_space = op.filter_samples(capture_space, tissue_status = "normal")
    normals = expand(
        CFG["dirs"]["normals"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}_passed.vcf.gz.tbi", 
        zip,
        seq_type = capture_space["seq_type"],
        genome_build = capture_space["genome_build"],
        normal_id = capture_space["sample_id"], 
        capture_space = capture_space["capture_space"])
    normals = list(dict.fromkeys(normals))
    omit_normals_list =  CFG["options"]["normals"]["omit_list"]
    if os.path.exists(omit_normals_list):
        with open(omit_normals_list) as file:
            remove = file.read().rstrip("\n").split("\n")
            remove = [s + ".tbi" for s in list]
            normals = [e for e in normals if e not in remove]
    return normals


rule _purecn_merge_normal_vcfs:
    input:
        normal = _get_normals_vcfs,
        normal_tbi = _get_normals_tbi
    output:
        normal_panel = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/normal_pon_vcf/{capture_space}_normalpanel.vcf.gz",
        normal_panel_tbi = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/normal_pon_vcf/{capture_space}_normalpanel.vcf.gz.tbi",
    conda: CFG["conda_envs"]["bcftools"]
    shell:
        """
            bcftools merge {input.normal} -Oz -o {output.normal_panel} --force-samples && 
            tabix -p vcf {output.normal_panel}
        """


rule _purecn_generate_samples_map:
    input:
        normal = _get_normals_vcfs
    output:
        map_sample = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/map/{capture_space}_samples_map.txt",
        done = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/map/.{capture_space}_samples_map.done"
    shell:
        """
            for samples in {input.normal}
            do
                name=$(basename $samples)
                name=${{name/.vcf.gz}}
                echo -e "$name\t$samples" >> {output.map_sample}
            done &&
            touch {output.done}
        """

# Create panel of normals database
if CFG['options']['new_normals'] == True:
    rule _purecn_gatk_genomicsDbimport:
        input:
            normal = _get_normals_vcfs,
            normal_tbi = _get_normals_tbi,
            map_sample = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/map/{capture_space}_samples_map.txt", 
            target_regions = str(rules._purecn_input_bed.output.bed),
            done = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/map/.{capture_space}_samples_map.done"
        output:
            touch(CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/genomicsdb/{capture_space}_database.done")
        params:
            mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
            pon_path = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/genomicsdb/{capture_space}_database/"
        conda: CFG["conda_envs"]["mutect"]
        resources: **CFG["resources"]["mutect"]
        shell:
            """
                gatk GenomicsDBImport --java-options "-Xmx{params.mem_mb}m" -L {input.target_regions} \
                --sample-name-map {input.map_sample} \
                --genomicsdb-workspace-path {params.pon_path} --lenient --merge-input-intervals TRUE \
                --overwrite-existing-genomicsdb-workspace TRUE 
            """

if CFG['options']['new_normals'] == True:
    rule _purecn_mutect2_create_pon:
        input:
            CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/genomicsdb/{capture_space}_database.done",
            fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
        output:
            pon = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}_mutect2_pon.vcf.gz"
        params:
            mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
            opts = "gendb://" + CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/genomicsdb/{capture_space}_database/"
        conda: CFG["conda_envs"]["mutect"]
        resources: **CFG["resources"]["mutect"]
        shell:
            """
                gatk CreateSomaticPanelOfNormals --java-options "-Xmx{params.mem_mb}m" --reference {params.fasta} --variant {params.opts} -O {output.pon}
            """


# -------------------------------------------------------------------------------------------------- #
# Part III - Set up tumour coverage and vcf files
# -------------------------------------------------------------------------------------------------- #

# Coverage files are dependent on if the file is BAM or CRAM:
# PureCN - by extension rsamtools does not have CRAM compatibility if R version is less than 4.2
# Work around is to use GATK to also calculate coverage and then feed it into pureCN for GC normalization

# Calculates coverage in intervals using GATK
rule _purecn_gatk_interval_list:
    input:
        intervals = CFG["dirs"]["inputs"] + "references/{genome_build}/{capture_space}/baits_{genome_build}_intervals.txt"
    output:
        gatk_intervals = CFG["dirs"]["inputs"] + "references/{genome_build}/{capture_space}/baits_{genome_build}_intervals_gatk.list"
    shell:
        """
            egrep -i  '^.*:.*-.*' {input.intervals} | awk '{{print $1}}' > {output.gatk_intervals}
        """


rule _purecn_gatk_interval_list_chrom:
    input:
        gatk_intervals = CFG["dirs"]["inputs"] + "references/{genome_build}/{capture_space}/baits_{genome_build}_intervals_gatk.list"
    output:
        chrom_int = CFG["dirs"]["inputs"] + "references/{genome_build}/{capture_space}/baits_{genome_build}_{chrom}.intervals_gatk.list"
    log:
        CFG["logs"]["inputs"] + "purecn_gatk_intervals/{genome_build}--{capture_space}/baits_{genome_build}_{chrom}.log"
    shell:
        op.as_one_line(
        """
            num_intervals=$( {{ egrep -i '^{wildcards.chrom}:.*-.*' {input.gatk_intervals} || true; }} | wc -l );
            if [[ $num_intervals -eq 0 ]]; then
                echo "No intervals found for chromosome {wildcards.chrom} in {input.gatk_intervals}" | tee {log}; 
                echo "{wildcards.chrom}:1-100" > {output.chrom_int}; 
            else
                echo "Found $num_intervals intervals for chromosome {wildcards.chrom} in {input.gatk_intervals}" | tee {log}; 
                egrep -i '^{wildcards.chrom}:.*-.*' {input.gatk_intervals} > {output.chrom_int}; 
            fi
        """
        )

rule _purecn_gatk_depthOfCoverage:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.crai",
        intervals =  CFG["dirs"]["inputs"] + "references/{genome_build}/{capture_space}/baits_{genome_build}_{chrom}.intervals_gatk.list"
    output:
        coverage = temp(CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.sample_interval_summary"),
        statistics = temp(CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.sample_interval_statistics")
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
        opts = CFG["options"]["mutect2"]["depth_coverage"],
        genome_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        base_name = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}"
    log: 
        CFG["logs"]["coverage"] + "gatk_coverage/{seq_type}--{genome_build}/{capture_space}/{sample_id}.{chrom}.log"
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    shell:
        """
            gatk DepthOfCoverage --java-options "-Xmx{params.mem_mb}m" {params.opts} --omit-depth-output-at-each-base --omit-locus-table --omit-per-sample-statistics --interval-merging-rule OVERLAPPING_ONLY -R {params.genome_fasta} -I {input.bam} -O {params.base_name} --intervals {input.intervals} > {log} 2>&1
        """

def _purecn_gatk_coverage_get_chr_depth(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    coverage = expand(
        CFG["dirs"]["coverage"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}/{{sample_id}}.{chrom}.sample_interval_summary", 
        chrom = chrs
    )
    return(coverage)
    
def _purecn_gatk_coverage_get_chr_statistics(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    statistics = expand(
        CFG["dirs"]["coverage"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}/{{sample_id}}.{chrom}.sample_interval_statistics", 
        chrom = chrs
    )
    return(statistics)

rule _purecn_gatk_coverage_concatenate_depths:
    input: 
        depth = _purecn_gatk_coverage_get_chr_depth,
        statistics = _purecn_gatk_coverage_get_chr_statistics,
    output: 
        depth = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.sample_interval_summary.gz"
    resources: **CFG["resources"]["concatenate_vcf"]
    shell: 
        op.as_one_line(
            """
            file1=$(echo {input.depth} | cut -d " " -f1 ) ;
            head -n 1 $file1 | gzip > {output.depth} ; 
            for sample in {input.depth};            
            do
                awk '(NR > 1)' $sample | gzip >> {output.depth} ;
            done
            """)


# Generate pureCN coverage files for all - used only in de novo pureCN
# Using GATK depthOfCoverage to generate all coverage files; pureCN used to normalize by GC
rule _purecn_coverage:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        intervals =  CFG["dirs"]["inputs"] + "references/{genome_build}/{capture_space}/baits_{genome_build}_intervals.txt",
        coverage = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.sample_interval_summary.gz"
    output:
        coverage = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}_coverage_loess.txt.gz"
    params:
        name = "{sample_id}",
        coverage_script = CFG["software"]["coverage_script"],
        outdir = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}",
        force = CFG["options"]["coverage"]["force"],
        opt = CFG["options"]["coverage"]["opts"],
        genome_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    conda: CFG["conda_envs"]["purecn"]
    resources:
        **CFG["resources"]["coverage"]
    threads: CFG["threads"]["coverage"]
    log:
        CFG["logs"]["coverage"] + "coverage/{seq_type}--{genome_build}/{capture_space}/{sample_id}.log"
    shell:
        """
            PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
            echo -e "Using {params.coverage_script} instead of default $PURECN/Coverage.R..."
            Rscript --vanilla {params.coverage_script}  --out-dir {params.outdir} \
            --bam {input.bam} \
            --name {params.name} \
            --reference {params.genome_fasta} \
            --coverage {input.coverage} \
            --intervals {input.intervals} {params.force} {params.opt} > {log} 2>&1 
        """

def _get_normals_coverage(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    capture_space = CFG["samples"][CFG["samples"]['capture_space'].isin([wildcards.capture_space])]
    capture_space = capture_space[capture_space["genome_build"].isin([wildcards.genome_build])]
    capture_space = capture_space[capture_space["seq_type"].isin([wildcards.seq_type])]
    capture_space = op.filter_samples(capture_space, tissue_status = "normal")
    normals = expand(
        CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}/{normal_id}_coverage_loess.txt.gz", 
        zip,
        seq_type = capture_space["seq_type"],
        genome_build = capture_space["genome_build"],
        normal_id = capture_space["sample_id"], 
        capture_space = capture_space["capture_space"])
    omit_normals_list =  CFG["options"]["normals_coverage_loess"]["omit_list"]
    if os.path.exists(omit_normals_list):
        with open(omit_normals_list) as file:
            remove = file.read().rstrip("\n").split("\n")
            normals = [e for e in normals if e not in remove]
    return normals
    
    
rule coverage_list_normals:
    input:
        coverage = _get_normals_coverage
    output:
        cov_list = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/cov_list.txt"
    shell:
        """
            ls -a {input.coverage} | cat > {output.cov_list} 
        """

# Run Mutect2 
# Tumour samples
rule _purecn_mutect2_tumour_germline:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        dbsnp = ancient(reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")),
        fasta = ancient(reference_files("genomes/{genome_build}/genome_fasta/genome.fa")),
        gnomad = ancient(reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")),
        pon = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}_mutect2_pon.vcf.gz"
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.vcf.gz.tbi"),
        stats = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.vcf.gz.stats"),
        f1r2 = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.f1r2.tar.gz"),
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
        opts = CFG["options"]["mutect2"]["mutect2_opts"]
    log: CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/mutect2_germline/{capture_space}/{tumour_id}/{chrom}.log"
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    shell:
        """
            gatk Mutect2 --java-options "-Xmx{params.mem_mb}m" {params.opts} --genotype-germline-sites true --genotype-pon-sites true --interval-padding 50 --germline-resource {input.gnomad} -R {input.fasta} -L {wildcards.chrom} -pon {input.pon} -I {input.bam} -O {output.vcf} --f1r2-tar-gz {output.f1r2} > {log} 2>&1
        """
        
        
def _mutect2_tumour_get_chr_vcf(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    mpileups = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{tumour_id}}/{{tumour_id}}.{chrom}.vcf.gz", 
        chrom = chrs
    )
    return(mpileups)
    
def _mutect2_tumour_get_chr_vcf_tbi(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    mpileups = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{tumour_id}}/{{tumour_id}}.{chrom}.vcf.gz.tbi", 
        chrom = chrs
    )
    return(mpileups)
    


# Merge chromosome mutect2 vcfs
rule _purecn_mutect2_tumour_concatenate_vcf:
    input: 
        vcf = _mutect2_tumour_get_chr_vcf,
        tbi = _mutect2_tumour_get_chr_vcf_tbi,
    output: 
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.vcf.gz.tbi"),
    resources: 
        **CFG["resources"]["concatenate_vcf"]
    conda:
        CFG["conda_envs"]["bcftools"]
    shell: 
        """
            bcftools concat {input.vcf} -Oz -o {output.vcf} && 
            tabix -p vcf {output.vcf}
        """


def _purecn_mutect2_get_chr_stats(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    stats = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{tumour_id}}/{{tumour_id}}.{chrom}.vcf.gz.stats", 
        chrom = chrs
    )
    return(stats)


# Merge chromosome mutect2 stats for FilterMutectCalls rule
rule _purecn_mutect2_merge_stats:
    input:
        stats = _purecn_mutect2_get_chr_stats
    output:
        stats = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.vcf.gz.stats"
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{tumour_id}/stats.log"
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    shell:
        op.as_one_line("""
        gatk MergeMutectStats $(for i in {input.stats}; do echo -n "-stats $i "; done)
        -O {output.stats} > {log} 2>&1
        """)

# Get pileup summaries
rule _purecn_mutect2_pileup_summaries: 
    input: 
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam", 
        snps = reference_files("genomes/{genome_build}/gatk/mutect2_small_exac.{genome_build}.vcf.gz"), 
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: 
        pileup = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/pileupSummary.table"
    log: 
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{tumour_id}/pileupSummary.log"
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    threads: CFG["threads"]["mutect2"]
    shell: 
        op.as_one_line("""
        gatk GetPileupSummaries 
            --java-options "-Xmx{params.mem_mb}m"
            -I {input.tumour_bam}
            -R {input.fasta} 
            -V {input.snps}
            -L {input.snps}
            -O {output.pileup}
            > {log} 2>&1
        """)

# Calculate contamination  
rule _purecn_mutect2_calc_contamination: 
    input: 
        pileup = str(rules._purecn_mutect2_pileup_summaries.output.pileup)
    output: 
        segments =  CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/segments.table", 
        contamination =  CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/contamination.table"
    log: 
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{tumour_id}/contam.log"
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    threads: CFG["threads"]["mutect2"]
    shell: 
        op.as_one_line("""
        gatk CalculateContamination 
            --java-options "-Xmx{params.mem_mb}m"
            -I {input.pileup}
            -tumor-segmentation {output.segments}
            -O {output.contamination}
            > {log} 2>&1
        """)
    
# Learn read orientation model
def _purecn_mutect2_get_chr_f1r2(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    f1r2 = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{tumour_id}}/{{tumour_id}}.{chrom}.f1r2.tar.gz", 
        chrom = chrs
    )
    return(f1r2)

rule _purecn_mutect2_learn_orient_model: 
    input: 
        f1r2 = _purecn_mutect2_get_chr_f1r2
    output:
        model =  CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/read-orientation-model.tar.gz"
    log: 
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{tumour_id}/read-orientation-model.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{tumour_id}/read-orientation-model.stderr.log"
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    threads: CFG["threads"]["mutect2"]
    shell: 
        op.as_one_line("""
        inputs=$(for input in {input.f1r2}; do printf -- "-I $input "; done);
        gatk LearnReadOrientationModel 
        --java-options "-Xmx{params.mem_mb}m" 
        $inputs -O {output.model}
        > {log.stdout} 2> {log.stderr}
        """)
    

# Marks variants filtered or PASS annotations
rule _purecn_mutect2_filter:
    input:
        vcf = str(rules._purecn_mutect2_tumour_concatenate_vcf.output.vcf),
        tbi = str(rules._purecn_mutect2_tumour_concatenate_vcf.output.tbi),
        stat = str(rules._purecn_mutect2_merge_stats.output.stats),
        segments = str(rules._purecn_mutect2_calc_contamination.output.segments), 
        contamination = str(rules._purecn_mutect2_calc_contamination.output.contamination), 
        model = str(rules._purecn_mutect2_learn_orient_model.output.model),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/output.unfilt.vcf.gz")
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{tumour_id}/mutect2_filter.log",
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["mutect2"]["mutect2_filter"]
    conda: CFG["conda_envs"]["mutect"]
    resources: **CFG["resources"]["mutect"]
    shell:
        op.as_one_line("""
        gatk FilterMutectCalls --java-options "-Xmx{params.mem_mb}m" 
            {params.opts} 
            -V {input.vcf} 
            -R {input.fasta}
            --tumor-segmentation {input.segments}
            --contamination-table {input.contamination}
            --ob-priors {input.model}
            -O {output.vcf} 
            > {log} 2>&1
        """)


# Filters for PASS variants
# This will only take somatic ones if filtering for PASSED (need to still maintain germline ones)
rule _purecn_mutect2_filter_passed:
    input:
        vcf = str(rules._purecn_mutect2_filter.output.vcf)
    output:
        vcf = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_passed.vcf.gz",
        tbi = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_passed.vcf.gz.tbi"
    params:
        filter_for_opts = CFG["options"]["mutect2"]["mutect2_filter_for"],
        filter_out_opts = CFG["options"]["mutect2"]["mutect2_filter_out"],
    log:
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/mutect2/{capture_space}/{tumour_id}/mutect2_passed.log"
    conda: CFG["conda_envs"]["bcftools"]
    resources: **CFG["resources"]["concatenate_vcf"]
    shell:
        op.as_one_line(""" 
        bcftools view "{params.filter_for_opts}" -e "{params.filter_out_opts}" -Oz -o {output.vcf} {input.vcf} 2> {log.stderr}
            &&
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

# -------------------------------------------------------------------------------------------------- #
# Part IV  - run pureCN
# -------------------------------------------------------------------------------------------------- #
# For CNVkit seg method

# setting a mapping bias database - used for cnvkit segs
if CFG['options']['new_normals'] == True:
    rule _purecn_set_normal:
        input:
            normal_panel = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/normal_pon_vcf/{capture_space}_normalpanel.vcf.gz",
            normal_panel_tbi = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/normal_pon_vcf/{capture_space}_normalpanel.vcf.gz.tbi"
        output:
            normal_db = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/purecn_cnvkit_normal/mapping_bias_{capture_space}_{genome_build}.rds"
        params:
            dirOut = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/purecn_cnvkit_normal/",
            genome = "{genome_build}",
            capture_space = "{capture_space}"
        conda: CFG["conda_envs"]["purecn"]
        resources:
            **CFG["resources"]["pureCN"]
        threads: CFG["threads"]["pureCN"]
        log:
            CFG["logs"]["pon"] + "{seq_type}--{genome_build}/purecn_cnvkit_normaldb/{capture_space}.log"
        shell:
            """
                echo $CONDA_DEFAULT_ENV ;
                PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
                mkdir -p {params.dirOut} ;
                Rscript --vanilla $PURECN/NormalDB.R --out-dir {params.dirOut} --normal-panel {input.normal_panel} \
                --assay {params.capture_space} --genome {params.genome} --force > {log} 2>&1 || true
            """

def _adjust_segMethod(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    if any(sample in str({wildcards.tumour_id}) for sample in CFG["options"]["purecn_cnvkit"]["flagged_samples"]):
        return 'Hclust'
    elif any(sample in str({wildcards.tumour_id}) for sample in CFG["options"]["purecn_cnvkit"]["failed_samples"]):
        return 'none'
    else:
        return 'PSCBS'

def _get_genome_build(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    if any(builds in str({wildcards.genome_build}) for builds in ['grch38', 'hg38']):
        return "hg38"
    if "38" in str({wildcards.genome_build}):
        return "hg38"
    if any(builds in str({wildcards.genome_build}) for builds in ['grch37', 'hg19', 'hs37d5']):
        return "hg19"
    if "19" in str({wildcards.genome_build}):
        return "hg19"

if CFG["cnvkit_seg"] == True:
    rule _purecn_run:
        input:
            cnr = str(rules._purecn_symlink_cnvkit_cnr.output.cnr),
            seg = str(rules._purecn_symlink_cnvkit_seg.output.seg),
            vcf = str(rules._purecn_mutect2_filter_passed.output.vcf),
            normal_db = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/purecn_cnvkit_normal/mapping_bias_{capture_space}_{genome_build}.rds",
            stats = str(rules._purecn_mutect2_merge_stats.output.stats),
            blacklist = str(rules._purecn_setup_blacklist_bed.output.blacklist)
        output:
            ploidy = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.csv",
            seg = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg",
            gene_cn = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_genes.csv",
            loh = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_loh.csv",
            rds = temp(CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.rds"),
            pdf = temp(CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.pdf")
        params:
            outdir = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/",
            sample_id = "{tumour_id}",
            genome_build = _get_genome_build,
            min_offtarget = CFG["options"]["pureCN"]["min_offtarget"],
            max_cn = CFG["options"]["pureCN"]["max_cn"],
            model = CFG["options"]["pureCN"]["model"],
            alpha = CFG["options"]["pureCN"]["alpha"],
            segmentation_method = _adjust_segMethod,
            opts = CFG["options"]["pureCN"]["opts"]
        conda: CFG["conda_envs"]["purecn"]
        resources:
            **CFG["resources"]["pureCN"]
        threads: CFG["threads"]["pureCN"]
        log:
            CFG["logs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.log"
        shell:
            """
                echo $CONDA_DEFAULT_ENV ;
                PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/
                Rscript --vanilla $PURECN/PureCN.R --out {params.outdir}  \
                    --sampleid {params.sample_id} \
                    --tumor {input.cnr} \
                    --seg-file {input.seg} \
                    --stats-file {input.stats} \
                    --vcf {input.vcf} \
                    --snp-blacklist {input.blacklist} \
                    --mappingbiasfile {input.normal_db} \
                    --genome {params.genome_build} \
                    --fun-segmentation {params.segmentation_method} \
                    --max-copy-number {params.max_cn} \
                    --min-fraction-offtarget {params.min_offtarget} \
                    --alpha {params.alpha} \
                    --model {params.model} \
                    --cores {threads} \
                    {params.opts}  > {log} 2>&1 
            """

if CFG["cnvkit_seg"] == True:
    rule _purecn_cleanup_xs:
        input:
            rds = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.rds",
            pdf = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.pdf",
        output:
            rds_removed = touch(CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.done")
        group: "cnvkit_post_process"
        shell:
            """
                rm {input.rds} ;
                rm {input.pdf}
            """

if CFG["cnvkit_seg"] == True:
    rule _purecn_fix_seg:
        input:
            purecn_native = str(rules._purecn_run.output.seg),
            rds_removed = str(rules._purecn_cleanup_xs.output.rds_removed)
        output:
            purecn_converted_seg = CFG["dirs"]["convert_coordinates"] + "purecn_cnvkit/fixed_seg/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg"
        params:
            tidy_pureCN_script = CFG["software"]["tidy_pureCN_script"],
            sample_id = "{tumour_id}"
        group: "cnvkit_post_process"
        shell:
            """
                {params.tidy_pureCN_script} -i {input.purecn_native} -o {output.purecn_converted_seg} -s {params.sample_id}
            """
 

# -------------------------------------------------------------------------------------------------- #
# For pureCN de novo PSCBS seg method using its own coverage files
if CFG['options']['new_normals'] == True:
    rule _purecn_set_normal_denovo:
        input:
            normal_panel = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/normal_pon_vcf/{capture_space}_normalpanel.vcf.gz",
            normal_panel_tbi = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/normal_pon_vcf/{capture_space}_normalpanel.vcf.gz.tbi",
            cov_list =  CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/cov_list.txt"
        output:
            normal_db = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/purecn_denovo_normal/mapping_bias_{capture_space}_{genome_build}.rds",
            normals = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/purecn_denovo_normal/normalDB_{capture_space}_{genome_build}.rds"
        params:
            dirOut = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/purecn_denovo_normal/",
            genome = "{genome_build}",
            platform = "{capture_space}"
        conda: CFG["conda_envs"]["purecn"]
        resources:
            **CFG["resources"]["pureCN"]
        threads: CFG["threads"]["pureCN"]
        log:
            CFG["logs"]["pon"] + "{seq_type}--{genome_build}/purecn_denovo_normaldb/{capture_space}.log"
        shell:
            """
                echo $CONDA_DEFAULT_ENV ;
                PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
                mkdir -p {params.dirOut} ;
                Rscript --vanilla $PURECN/NormalDB.R --out-dir {params.dirOut} --normal-panel {input.normal_panel} \
                --coverage-files {input.cov_list} \
                --assay {params.platform} --genome {params.genome} --force > {log} 2>&1 || true
            """


def _adjust_segMethod_denovo(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    if any(sample in str({wildcards.tumour_id}) for sample in CFG["options"]["purecn_denovo"]["flagged_samples"]):
        return 'CBS'
    elif any(sample in str({wildcards.tumour_id}) for sample in CFG["options"]["purecn_denovo"]["failed_samples"]):
        return 'none'
    else:
        return 'PSCBS'

rule _purecn_denovo_run:
    input:
        cnr = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_coverage_loess.txt.gz",
        vcf = str(rules._purecn_mutect2_filter_passed.output.vcf),
        normal_db = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/purecn_denovo_normal/mapping_bias_{capture_space}_{genome_build}.rds",
        normals = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/purecn_denovo_normal/normalDB_{capture_space}_{genome_build}.rds",
        stats = str(rules._purecn_mutect2_merge_stats.output.stats),
        blacklist = str(rules._purecn_setup_blacklist_bed.output.blacklist),
        intervals = str(rules._purecn_setinterval.output.intervals)
    output:
        ploidy = CFG["dirs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.csv",
        seg = CFG["dirs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg",
        loh = CFG["dirs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_loh.csv",
        rds = temp(CFG["dirs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.rds"),
        pdf = temp(CFG["dirs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.pdf"),
    params:
        outdir = CFG["dirs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/",
        sample_id = "{tumour_id}",
        genome_build = _get_genome_build,
        min_offtarget = CFG["options"]["pureCN"]["min_offtarget"],
        max_cn = CFG["options"]["pureCN"]["max_cn"],
        model = CFG["options"]["pureCN"]["model"],
        alpha = CFG["options"]["pureCN"]["alpha"],
        segmentation_method = _adjust_segMethod_denovo,
        opts = CFG["options"]["pureCN"]["opts"]
    conda: CFG["conda_envs"]["purecn"]
    resources:
        **CFG["resources"]["pureCN"]
    threads: CFG["threads"]["pureCN"]
    log:
        CFG["logs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.log"
    shell:
        """
            echo $CONDA_DEFAULT_ENV ;
            PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/
            Rscript --vanilla $PURECN/PureCN.R --out {params.outdir}  \
                --sampleid {params.sample_id} \
                --tumor {input.cnr} \
                --stats-file {input.stats} \
                --vcf {input.vcf} \
                --snp-blacklist {input.blacklist} \
                --mappingbiasfile {input.normal_db} \
                --normaldb {input.normals} \
                --genome {params.genome_build} \
                --intervals {input.intervals} \
                --fun-segmentation {params.segmentation_method} \
                --max-copy-number {params.max_cn} \
                --min-fraction-offtarget {params.min_offtarget} \
                --alpha {params.alpha} \
                --model betabin \
                --cores {threads} \
                {params.opts} > {log} 2>&1 
        """

rule _purecn_denovo_cleanup_xs:
    input:
        rds = CFG["dirs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.rds",
        pdf = CFG["dirs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.pdf",
    output:
        rds_removed = touch(CFG["dirs"]["pureCN"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.done")
    group: "purecn_post_process"
    shell:
        """
            rm {input.rds} ;
            rm {input.pdf}
        """

rule _purecn_denovo_fix_seg:
    input:
        purecn_native = str(rules._purecn_denovo_run.output.seg),
        rds_removed = str(rules._purecn_denovo_cleanup_xs.output.rds_removed)
    output:
        purecn_converted_seg = CFG["dirs"]["convert_coordinates"] + "purecn_denovo/fixed_seg/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg"
    params:
        tidy_pureCN_script = CFG["software"]["tidy_pureCN_script"],
        sample_id = "{tumour_id}"
    group: "purecn_post_process"
    shell:
        """
            {params.tidy_pureCN_script} -i {input.purecn_native} -o {output.purecn_converted_seg} -s {params.sample_id}
        """
        
# -------------------------------------------------------------------------------------------------- #
# Part V  - Project to other genome builds and remove capture_space wildcard
# -------------------------------------------------------------------------------------------------- #
# Projection to other genome builds
def _purecn_get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")


# Convert the coordinates of seg file to a different genome build
if CFG["cnvkit_seg"] == True:
    rule _purecn_convert_coordinates:
        input:
            purecn_native = str(rules._purecn_fix_seg.output.purecn_converted_seg),
            purecn_chain = _purecn_get_chain
        output:
            purecn_lifted = CFG["dirs"]["convert_coordinates"] + "purecn_cnvkit/from--{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.seg"
        log:
            stderr = CFG["logs"]["convert_coordinates"] + "purecn_cnvkit/from--{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.stderr.log"
        threads: 1
        params:
            liftover_script = CFG["options"]["liftover_script_path"],
            liftover_minmatch = CFG["options"]["liftover_minMatch"]
        conda:
            CFG["conda_envs"]["liftover"]
        group: "cnvkit_post_process"
        shell:
            op.as_one_line("""
            echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
            bash {params.liftover_script}
            SEG
            {input.purecn_native}
            {output.purecn_lifted}
            {input.purecn_chain}
            YES
            {params.liftover_minmatch}
            2>> {log.stderr}
            """)

# de novo pureCN
rule _purecn_denovo_convert_coordinates:
    input:
        purecn_native = str(rules._purecn_denovo_fix_seg.output.purecn_converted_seg),
        purecn_chain = _purecn_get_chain
    output:
        purecn_lifted = CFG["dirs"]["convert_coordinates"] + "purecn_denovo/from--{seq_type}--{genome_build}/{capture_space}/purecn_cnvkit_seg/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["convert_coordinates"] + "purecn_denovo/from--{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.stderr.log"
    threads: 1
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        liftover_minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    group: "purecn_post_process"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.liftover_script}
        SEG
        {input.purecn_native}
        {output.purecn_lifted}
        {input.purecn_chain}
        YES
        {params.liftover_minmatch}
        2>> {log.stderr}
        """)


def _purecn_prepare_projection(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    this_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"].tolist()
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    prefixed_projections = CFG["options"]["prefixed_projections"]
    non_prefixed_projections = CFG["options"]["non_prefixed_projections"]

    if any(substring in this_genome_build[0] for substring in prefixed_projections):
        hg38_projection = str(rules._purecn_fix_seg.output.purecn_converted_seg).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        grch37_projection = str(rules._purecn_convert_coordinates.output.purecn_lifted).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        # handle the hg19 (prefixed) separately
        if "38" in str(this_genome_build[0]):
            grch37_projection = grch37_projection.replace("{chain}", "hg38ToHg19")
        else:
            grch37_projection = grch37_projection.replace("{chain}", "hg19ToHg38")

    elif any(substring in this_genome_build[0] for substring in non_prefixed_projections):
        grch37_projection = str(rules._purecn_fix_seg.output.purecn_converted_seg).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        hg38_projection = str(rules._purecn_convert_coordinates.output.purecn_lifted).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        # handle the grch38 (non-prefixed) separately
        if "38" in str(this_genome_build[0]):
            hg38_projection = hg38_projection.replace("{chain}", "hg38ToHg19")
        else:
            hg38_projection = hg38_projection.replace("{chain}", "hg19ToHg38")
    else:
        raise AttributeError(f"The specified genome build {this_genome_build[0]} is not specified in the config under options to indicate its chr prefixing.")
    
    return{
        "grch37_projection": grch37_projection,
        "hg38_projection": hg38_projection
    }
    
    
# Fill the missing segments of seg files with neutral regions to complete the genome coverage
if CFG["cnvkit_seg"] == True:
    rule _purecn_fill_segments:
        input:
            unpack(_purecn_prepare_projection)
        output:
            grch37_filled = temp(CFG["dirs"]["fill_regions"] + "purecn_cnvkit/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"),
            hg38_filled = temp(CFG["dirs"]["fill_regions"] + "purecn_cnvkit/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg"),
        log:
            stderr = CFG["logs"]["fill_regions"] + "{seq_type}--projection/purecn_cnvkit/{tumour_id}--{normal_id}--{pair_status}.{tool}_fill_segments.stderr.log"
        threads: 1
        params:
            path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"]
        conda:
            CFG["conda_envs"]["bedtools"]
        group: "cnvkit_post_process"
        shell:
            op.as_one_line("""
            echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
            echo "Filling grch37 projection" >> {log.stderr};
            bash {params.path}fill_segments.sh
            {params.path}src/chromArm.grch37.bed
            {input.grch37_projection}
            {params.path}src/blacklisted.grch37.bed
            {output.grch37_filled}
            {wildcards.tumour_id}
            SEG
            2>> {log.stderr};
            echo "Filling hg38 projection" >> {log.stderr};
            bash {params.path}fill_segments.sh
            {params.path}src/chromArm.hg38.bed
            {input.hg38_projection}
            {params.path}src/blacklisted.hg38.bed
            {output.hg38_filled}
            {wildcards.tumour_id}
            SEG
            2>> {log.stderr};
            """)


def _purecn_denovo_prepare_projection(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    this_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"].tolist()
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    prefixed_projections = CFG["options"]["prefixed_projections"]
    non_prefixed_projections = CFG["options"]["non_prefixed_projections"]

    if any(substring in this_genome_build[0] for substring in prefixed_projections):
        hg38_denovo_projection = str(rules._purecn_denovo_fix_seg.output.purecn_converted_seg).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        grch37_denovo_projection = str(rules._purecn_denovo_convert_coordinates.output.purecn_lifted).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        # handle the hg19 (prefixed) separately
        if "38" in str(this_genome_build[0]):
            grch37_denovo_projection = grch37_denovo_projection.replace("{chain}", "hg38ToHg19")
        else:
            grch37_denovo_projection = grch37_denovo_projection.replace("{chain}", "hg19ToHg38")

    elif any(substring in this_genome_build[0] for substring in non_prefixed_projections):
        grch37_denovo_projection = str(rules._purecn_denovo_fix_seg.output.purecn_converted_seg).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        hg38_denovo_projection = str(rules._purecn_denovo_convert_coordinates.output.purecn_lifted).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        # handle the grch38 (non-prefixed) separately
        if "38" in str(this_genome_build[0]):
            hg38_denovo_projection = hg38_denovo_projection.replace("{chain}", "hg38ToHg19")
        else:
            hg38_denovo_projection = hg38_denovo_projection.replace("{chain}", "hg19ToHg38")
    else:
        raise AttributeError(f"The specified genome build {this_genome_build[0]} is not specified in the config under options to indicate its chr prefixing.")

    return{
        "grch37_denovo_projection": grch37_denovo_projection,
        "hg38_denovo_projection": hg38_denovo_projection,
    }


rule _purecn_denovo_fill_segments:
    input:
        unpack(_purecn_denovo_prepare_projection)
    output:
        grch37_denovo_filled = temp(CFG["dirs"]["fill_regions"] + "purecn_denovo/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"),
        hg38_denovo_filled = temp(CFG["dirs"]["fill_regions"] + "purecn_denovo/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg")
    log:
        stderr = CFG["logs"]["fill_regions"] + "{seq_type}--projection/purecn_denovo/{tumour_id}--{normal_id}--{pair_status}.{tool}_fill_segments.stderr.log"
    threads: 1
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"]
    conda:
        CFG["conda_envs"]["bedtools"]
    group: "purecn_post_process"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        echo "Filling grch37 de novo projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.grch37.bed
        {input.grch37_denovo_projection}
        {params.path}src/blacklisted.grch37.bed
        {output.grch37_denovo_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        echo "Filling hg38 de novo projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.hg38.bed
        {input.hg38_denovo_projection}
        {params.path}src/blacklisted.hg38.bed
        {output.hg38_denovo_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        """)


def _purecn_determine_projection(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    if any(substring in wildcards.projection for substring in ["hg19", "grch37", "hs37d5"]):
        this_file = CFG["dirs"]["fill_regions"] + "purecn_cnvkit/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"
    elif any(substring in wildcards.projection for substring in ["hg38", "grch38"]):
        this_file = CFG["dirs"]["fill_regions"] + "purecn_cnvkit/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg"
    return (this_file)


def _purecn_denovo_determine_projection(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    if any(substring in wildcards.projection for substring in ["hg19", "grch37", "hs37d5"]):
        this_file = CFG["dirs"]["fill_regions"] + "purecn_denovo/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"
    elif any(substring in wildcards.projection for substring in ["hg38", "grch38"]):
        this_file = CFG["dirs"]["fill_regions"] + "purecn_denovo/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg"
    return (this_file)


# Normalize chr prefix of the output file
if CFG["cnvkit_seg"] == True:
    rule _purecn_normalize_projection:
        input:
            filled = _purecn_determine_projection,
            chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
        output:
            projection = CFG["dirs"]["normalize"] + "purecn_cnvkit/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
        resources:
            **CFG["resources"]["post_purecn"]
        threads: 1
        group: "cnvkit_post_process"
        wildcard_constraints: 
            projection = "|".join(CFG["output"]["requested_projections"]),
            # purecn_version = "|".join(CFG["output"]["purecn_versions"]),
            tool = "purecn"
        run:
            # read the main chromosomes file of the projection
            chromosomes = pd.read_csv(input.chrom_file, sep = "\t", names=["chromosome"], header=None)
            # handle chr prefix
            if "chr" in chromosomes["chromosome"][0]:
                seg_open = pd.read_csv(input.filled, sep = "\t")
                chrom = list(seg_open['chrom'])
                # avoid cases of chrchr1 if the prefix already there
                for i in range(len(chrom)):
                    if 'chr' not in str(chrom[i]):
                        chrom[i]='chr'+str(chrom[i])
                seg_open.loc[:, 'chrom']=chrom
                seg_open.to_csv(output.projection, sep="\t", index=False)
            else:
                # remove chr prefix
                seg_open = pd.read_csv(input.filled, sep = "\t")
                seg_open["chrom"] = seg_open["chrom"].astype(str).str.replace('chr', '')
                seg_open.to_csv(output.projection, sep="\t", index=False)

rule _purecn_denovo_normalize_projection:
    input:
        filled = _purecn_denovo_determine_projection,
        chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
    output:
        projection = CFG["dirs"]["normalize"] + "purecn_denovo/seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    resources:
        **CFG["resources"]["post_purecn"]
    threads: 1
    group: "purecn_post_process"
    wildcard_constraints: 
        projection = "|".join(CFG["output"]["requested_projections"]),
        # purecn_version = "|".join(CFG["output"]["purecn_versions"]),
        tool = "purecn"
    run:
        # read the main chromosomes file of the projection
        chromosomes = pd.read_csv(input.chrom_file, sep = "\t", names=["chromosome"], header=None)
        # handle chr prefix
        if "chr" in chromosomes["chromosome"][0]:
            seg_open = pd.read_csv(input.filled, sep = "\t")
            chrom = list(seg_open['chrom'])
            # avoid cases of chrchr1 if the prefix already there
            for i in range(len(chrom)):
                if 'chr' not in str(chrom[i]):
                    chrom[i]='chr'+str(chrom[i])
            seg_open.loc[:, 'chrom']=chrom
            seg_open.to_csv(output.projection, sep="\t", index=False)
        else:
            # remove chr prefix
            seg_open = pd.read_csv(input.filled, sep = "\t")
            seg_open["chrom"] = seg_open["chrom"].astype(str).str.replace('chr', '')
            seg_open.to_csv(output.projection, sep="\t", index=False)

# Symlinks the final output files into the module results directory (under '99-outputs/')
if CFG["cnvkit_seg"] == True:
    rule _purecn_cnvkit_output_projection:
        input:
            projection = str(rules._purecn_normalize_projection.output.projection)
        output:
            projection = CFG["output"]["cnvkit_seg"]["projection"]
        threads: 1
        group: "cnvkit_post_process"
        wildcard_constraints: 
            projection = "|".join(CFG["output"]["requested_projections"]), 
            pair_status = "|".join(set(CFG["runs"]["pair_status"].tolist())),
            purecn_version = "|".join(CFG["output"]["purecn_versions"]),
            tool = "purecn"
        run:
            op.relative_symlink(input.projection, output.projection, in_module = True)

rule _purecn_denovo_output_projection:
    input:
        projection_denovo = str(rules._purecn_denovo_normalize_projection.output.projection),
    output:
        projection_denovo = CFG["output"]["denovo_seg"]["projection"]
    threads: 1
    group: "purecn_post_process"
    wildcard_constraints: 
        projection = "|".join(CFG["output"]["requested_projections"]), 
        pair_status = "|".join(set(CFG["runs"]["pair_status"].tolist())),
        purecn_version = "|".join(CFG["output"]["purecn_versions"]),
        tool = "purecn"
    run:
        op.relative_symlink(input.projection_denovo, output.projection_denovo, in_module = True)

# Output files without a capture_space wildcard in them for consistency
def _purecn_cnvkit_drop_capture_space_wc(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    cnvkit_ploidy = str(rules._purecn_run.output.ploidy).replace("{capture_space}", this_space[0])
    cnvkit_gene_cn = str(rules._purecn_run.output.gene_cn).replace("{capture_space}", this_space[0])
    cnvkit_loh = str(rules._purecn_run.output.loh).replace("{capture_space}", this_space[0])
    return{
        "cnvkit_ploidy": cnvkit_ploidy,
        "cnvkit_gene_cn": cnvkit_gene_cn,
        "cnvkit_loh": cnvkit_loh
    }

def _purecn_denovo_drop_capture_space_wc(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    denovo_ploidy = str(rules._purecn_denovo_run.output.ploidy).replace("{capture_space}", this_space[0])
    denovo_loh = str(rules._purecn_denovo_run.output.loh).replace("{capture_space}", this_space[0])
    return{
        "denovo_ploidy": denovo_ploidy,
        "denovo_loh": denovo_loh
    }


# Symlinks the final output files into the module results directory (under '99-outputs/')
if CFG["cnvkit_seg"] == True:
    rule _purecn_cnvkit_output_files:
        input:
            unpack(_purecn_cnvkit_drop_capture_space_wc)
        output:
            cnvkit_ploidy = CFG["output"]["cnvkit_ploidy"]["info"],
            cnvkit_gene_cn = CFG["output"]["cnvkit_gene_cn"]["cnvkit_gene_cn"],
            cnvkit_loh = CFG["output"]["cnvkit_loh"]["cnvkit_loh"]
        group: "cnvkit_post_process"
        wildcard_constraints: 
            projection = "|".join(CFG["output"]["requested_projections"]), 
            pair_status = "|".join(set(CFG["runs"]["pair_status"].tolist())),
            purecn_version = "|".join(CFG["output"]["purecn_versions"])
        run:
            op.relative_symlink(input.cnvkit_ploidy, output.cnvkit_ploidy, in_module = True)
            op.relative_symlink(input.cnvkit_gene_cn, output.cnvkit_gene_cn, in_module = True)
            op.relative_symlink(input.cnvkit_loh, output.cnvkit_loh, in_module = True)
        
rule _purecn_denovo_output_files:
    input:
        unpack(_purecn_denovo_drop_capture_space_wc)
    output:
        denovo_ploidy = CFG["output"]["denovo_ploidy"]["info"],
        denovo_loh = CFG["output"]["denovo_loh"]["denovo_loh"]
    group: "purecn_post_process"
    wildcard_constraints: 
        projection = "|".join(CFG["output"]["requested_projections"]), 
        pair_status = "|".join(set(CFG["runs"]["pair_status"].tolist())),
        purecn_version = "|".join(CFG["output"]["purecn_versions"])
    run:
        op.relative_symlink(input.denovo_ploidy, output.denovo_ploidy, in_module = True)
        op.relative_symlink(input.denovo_loh, output.denovo_loh, in_module = True)


def _purecn_cnvkit_drop_capture_space_wc_seg(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    cnvkit_seg = str(rules._purecn_fix_seg.output.purecn_converted_seg).replace("{capture_space}", this_space[0])

    return{
        "cnvkit_seg": cnvkit_seg
    }

def _purecn_denovo_drop_capture_space_wc_seg(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    denovo_seg = str(rules._purecn_denovo_fix_seg.output.purecn_converted_seg).replace("{capture_space}", this_space[0])

    return{
        "denovo_seg": denovo_seg
    }
    
if CFG["cnvkit_seg"] == True:
    rule _purecn_cnvkit_output_seg:
        input:
            unpack(_purecn_cnvkit_drop_capture_space_wc_seg)
        output:
            cnvkit_seg = CFG["output"]["cnvkit_seg"]["original"]
        group: "cnvkit_post_process"
        wildcard_constraints: 
            projection = "|".join(CFG["output"]["requested_projections"]), 
            pair_status = "|".join(set(CFG["runs"]["pair_status"].tolist())),
            purecn_version = "|".join(CFG["output"]["purecn_versions"])
        run:
            op.relative_symlink(input.cnvkit_seg, output.cnvkit_seg, in_module = True)

rule _purecn_denovo_output_seg:
    input:
        unpack(_purecn_denovo_drop_capture_space_wc_seg)
    output:
        denovo_seg = CFG["output"]["denovo_seg"]["original"],
    group: "purecn_post_process"
    wildcard_constraints: 
        projection = "|".join(CFG["output"]["requested_projections"]), 
        pair_status = "|".join(set(CFG["runs"]["pair_status"].tolist())),
        purecn_version = "|".join(CFG["output"]["purecn_versions"])
    run:
        op.relative_symlink(input.denovo_seg, output.denovo_seg, in_module = True)

# Select the seg file that had the least amount of deviation (noise) - as measured by the MAD value
def _purecn_take_lowest_MAD(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    projection = CFG["output"]["requested_projections"]

    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    cnvkit_seg = CFG["dirs"]["outputs"] + "purecn_cnvkit/seg/" + wildcards.seq_type + "--projection/" + wildcards.tumour_id + "--" + wildcards.normal_id + "--" + wildcards.pair_status + "." + wildcards.tool + "." + wildcards.projection + ".seg"
    
    denovo_seg = CFG["dirs"]["outputs"] + "purecn_denovo/seg/" + wildcards.seq_type + "--projection/" + wildcards.tumour_id + "--" + wildcards.normal_id + "--" + wildcards.pair_status + "." + wildcards.tool + "." + wildcards.projection + ".seg"
    
    
    cnvkit_dir = glob.glob(CFG["dirs"]["pureCN_cnvkit"] + wildcards.seq_type + "--*/*/" + wildcards.tumour_id)
    denovo_dir = glob.glob(CFG["dirs"]["pureCN"] + wildcards.seq_type + "--*/*/" + wildcards.tumour_id)
    
    if any(sample in str({wildcards.tumour_id}) for sample in CFG["output"]["best_seg_manual"]["purecn_denovo"]):
        if os.path.isfile(str(denovo_seg)):
            best_seg = denovo_seg
            
    elif any(sample in str({wildcards.tumour_id}) for sample in CFG["output"]["best_seg_manual"]["purecn_cnvkit"]):
        if os.path.isfile(str(cnvkit_seg)):
            best_seg = cnvkit_seg
    
    if os.path.isfile(str(cnvkit_seg)) and os.path.isfile(str(denovo_seg)):
        
        cnvkit_mapd = list()
        with open(cnvkit_dir[0] + "/" + wildcards.tumour_id + ".log") as file:
            for line in file.readlines():
                if 'Mean standard deviation of log-ratios' in line:
                    cnvkit_mapd.append(str(line))
        
        # Take the 8th element - representing the MAD score
        cnvkit_mapd_value = float(cnvkit_mapd[len(cnvkit_mapd)-1].split()[8])

        denovo_mapd = list()
        with open(denovo_dir[0] + "/" + wildcards.tumour_id + ".log") as file:
            for line in file.readlines():
                if 'Mean standard deviation of log-ratios' in line:
                    denovo_mapd.append(str(line))
                    
        # Take the 8th element - representing the MAD score
        denovo_mapd_value = float(denovo_mapd[len(denovo_mapd)-1].split()[8])
        
        if (denovo_mapd_value < cnvkit_mapd_value):
            # Account for seg files that are over-segmented but have a lower MAD
            count_denovo = []
            with open(denovo_seg) as file:
                next(file)
                for line in file:
                    line = line.rstrip('\n').rstrip('\r')
                    cols = line.split('\t')
                    if (float(cols[5]) > 0.5 or float(cols[5]) <-0.5 ):
                        count_denovo.append(line)
            count_denovo = len(count_denovo)
            
            count_cnvkit = []
            with open(cnvkit_seg) as file:
                next(file)
                for line in file:
                    line = line.rstrip('\n').rstrip('\r')
                    cols = line.split('\t')
                    if (float(cols[5]) > 0.5 or float(cols[5]) <-0.5 ):
                        count_cnvkit.append(line)
            count_cnvkit = len(count_cnvkit)
            if (count_cnvkit > 0):
                if (float(count_denovo/count_cnvkit) > 2.8 and count_denovo > 400):
                    best_seg = str(cnvkit_seg)
                else:
                    best_seg = str(denovo_seg)
            else:
                best_seg = str(denovo_seg)
        else:
            best_seg = str(cnvkit_seg)
        
    elif os.path.isfile(str(cnvkit_seg)):
        best_seg = cnvkit_seg
        
    else:
        best_seg = denovo_seg
        
    return{
        "best_seg": best_seg
    }


rule _purecn_best_seg:
    input:
        unpack(_purecn_take_lowest_MAD)
    output:
        best_seg = CFG["output"]["best_seg"]
    wildcard_constraints:
        projection = "|".join(CFG["output"]["requested_projections"]), 
        pair_status = "|".join(set(CFG["runs"]["pair_status"].tolist())),
        purecn_version = "|".join(CFG["output"]["purecn_versions"])
    run:
        op.relative_symlink(input.best_seg, output.best_seg, in_module = True)

# Generates the target sentinels for each run, which generate the symlinks
if CFG["cnvkit_seg"] == True:
    rule _purecn_all:
        input:
            expand(
                [
                    str(rules._purecn_cnvkit_output_seg.output.cnvkit_seg),
                    str(rules._purecn_cnvkit_output_files.output.cnvkit_ploidy),
                    str(rules._purecn_cnvkit_output_files.output.cnvkit_gene_cn),
                    str(rules._purecn_cnvkit_output_files.output.cnvkit_loh),
                    str(rules._purecn_denovo_output_seg.output.denovo_seg),
                    str(rules._purecn_denovo_output_files.output.denovo_ploidy),
                    str(rules._purecn_denovo_output_files.output.denovo_loh)
                ],
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["runs"]["tumour_seq_type"],
                genome_build=CFG["runs"]["tumour_genome_build"],
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                pair_status=CFG["runs"]["pair_status"]
            ),
            expand(
                expand(
                [
                    str(rules._purecn_cnvkit_output_projection.output.projection),
                    str(rules._purecn_denovo_output_projection.output.projection_denovo)
                ],
                zip,  # Run expand() with zip(), not product()
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                seq_type=CFG["runs"]["tumour_seq_type"],
                pair_status=CFG["runs"]["pair_status"],
                allow_missing=True),
                tool = "purecn",
                projection=CFG["output"]["requested_projections"],
                purecn_version=CFG["output"]["purecn_versions"]
            ),
            expand(
                expand(
                [
                    str(rules._purecn_best_seg.output.best_seg)
                ],
                zip,  # Run expand() with zip(), not product()
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                seq_type=CFG["runs"]["tumour_seq_type"],
                pair_status=CFG["runs"]["pair_status"],
                allow_missing=True),
                tool = "purecn",
                projection=CFG["output"]["requested_projections"])

if CFG["cnvkit_seg"] == False:
    rule _purecn_denovo_all:
        input:
            expand(
                [
                    str(rules._purecn_denovo_output_seg.output.denovo_seg),
                    str(rules._purecn_denovo_output_files.output.denovo_ploidy),
                    str(rules._purecn_denovo_output_files.output.denovo_loh)
                ],
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["runs"]["tumour_seq_type"],
                genome_build=CFG["runs"]["tumour_genome_build"],
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                pair_status=CFG["runs"]["pair_status"]
            ),
            expand(
                expand(
                [
                    str(rules._purecn_denovo_output_projection.output.projection_denovo)
                ],
                zip,  # Run expand() with zip(), not product()
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                seq_type=CFG["runs"]["tumour_seq_type"],
                pair_status=CFG["runs"]["pair_status"],
                allow_missing=True),
                tool = "purecn",
                projection=CFG["output"]["requested_projections"],
                purecn_version=CFG["output"]["purecn_versions"]
            )
            # expand(
            #     expand(
            #     [
            #         str(rules._purecn_best_seg.output.best_seg)
            #     ],
            #     zip,  # Run expand() with zip(), not product()
            #     tumour_id=CFG["runs"]["tumour_sample_id"],
            #     normal_id=CFG["runs"]["normal_sample_id"],
            #     seq_type=CFG["runs"]["tumour_seq_type"],
            #     pair_status=CFG["runs"]["pair_status"],
            #     allow_missing=True),
            #     tool = "purecn",
            #     projection=CFG["output"]["requested_projections"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
