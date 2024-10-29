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
# `CFG` is a shortcut to `config["lcr-modules"]["whatshap"]`
CFG = op.setup_module(
    name = "whatshap",
    version = "1.0",
    subdirectories = ["inputs", "phase_vcf" , "phase_bam", "split_bam", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _whatshap_input_bam,
    _whatshap_input_vcf,
    _whatshap_input_regions,
    _whatshap_input_chrs,
    _whatshap_output_phased,
    _whatshap_output_split, 
    _whatshap_output_vcf,
    _whatshap_all

##### SAMPLES TABLES #####
# This module can use a mix of seq_types
# Split the samples table in two based on specified 
# seq_types from the config. 

VCF_SAMPLES = op.filter_samples(CFG["samples"], seq_type = CFG["options"]["vcf_seq_type"])
BAM_SAMPLES = op.filter_samples(CFG["samples"], seq_type = CFG["options"]["bam_seq_type"])

# Ensure the genome_build for each biopsy in VCF_SAMPLES matches the genome_build in BAM_SAMPLES
import pandas as pd

def update_genome_build(vcf_samples: pd.DataFrame, bam_samples: pd.DataFrame) -> pd.DataFrame:
    # Create a dictionary from bam_samples for quick lookup
    bam_genome_build = bam_samples.set_index('biopsy_id')['genome_build'].to_dict()
    
    # Update genome_build in vcf_samples
    vcf_samples['genome_build'] = vcf_samples['biopsy_id'].map(bam_genome_build)
    
    # Check for any biopsy_id in vcf_samples not found in bam_samples
    if vcf_samples['genome_build'].isnull().any():
        missing_ids = vcf_samples[vcf_samples['genome_build'].isnull()]['biopsy_id'].tolist()
        raise ValueError(f"biopsy_id(s) {missing_ids} in vcf_samples not found in bam_samples")
    
    return vcf_samples

VCF_SAMPLES = update_genome_build(VCF_SAMPLES, BAM_SAMPLES)

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _whatshap_input_bam:
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

rule _whatshap_input_vcf:
    input:
        vcf = CFG["inputs"]["sample_vcf"],
        tbi = CFG["inputs"]["sample_vcf_tbi"]
    output:
        vcf = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.vcf.gz",
        tbi = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.vcf.gz.tbi"
    run:
        op.absolute_symlink(input.vcf, output.vcf)
        op.absolute_symlink(input.tbi, output.tbi)

# Populate the whole_genome key with the path to the genome bed file if no other region is specified        
if "whole_genome" in CFG["inputs"]["regions_bed"].keys() or not CFG["inputs"]["regions_bed"]: 
    CFG["inputs"]["regions_bed"]["whole_genome"] = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.bed")

rule _whatshap_input_regions: 
    input: 
        bed = lambda w: expand(config["lcr-modules"]["whatshap"]["inputs"]["regions_bed"][w.regions_bed], genome_build = w.genome_build)
    output: 
        regions = CFG["dirs"]["inputs"] + "regions_bed/{regions_bed}.{genome_build}.regions.txt"
    params: 
        padding = CFG["options"]["region_padding"]
    run: 
        regions_list = []
        with open(input.bed[0]) as f, open(output.regions, "w") as o:
            for line in f:
                line = line.rstrip("\n").rstrip("\r")  # Handle line endings
                # Split by tab delimiter
                cols = line.split("\t")
                # Skip header lines if present
                if line.startswith("chrom") or line.startswith("Chrom"): 
                    continue
                try:
                    chrom = cols[0]
                    start = int(cols[1])
                    end = int(cols[2])
                except (IndexError, ValueError) as e:
                    raise AttributeError("Input bed file '%s' appears to be malformed" % input.bed) from e
                region = f"{chrom}:{max(0, start - params.padding)}-{end + params.padding}"
                regions_list.append(region)
            o.write(" ".join(regions_list))
    
# Split vcf phasing per-chromosome
checkpoint _whatshap_input_chrs:
    input:
        chrs = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        chrs = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.absolute_symlink(input.chrs, output.chrs)
        
rule _whatshap_split_vcf:
    input:
        vcf = str(rules._whatshap_input_vcf.output.vcf)
    output:
        vcf = temp(CFG["dirs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.vcf.gz.tbi")
    conda:
        CFG["conda_envs"]["bcftools"]
    threads: CFG["threads"]["split_vcf"]
    resources: **CFG["resources"]["split_vcf"]
    log:
        stderr = CFG["logs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.split_vcf.stderr.log"
    wildcard_constraints: 
        sample_id = "|".join(VCF_SAMPLES["sample_id"].unique())
    group: "split_and_phase"
    shell:
        op.as_one_line("""
            bcftools view -Oz -o {output.vcf} -r {wildcards.chrom} --threads {threads} {input.vcf}  2> {log.stderr} &&
            bcftools index -t {output.vcf}
        """)

def get_vcf(w):
    CFG = config["lcr-modules"]["whatshap"]
    # Get the biopsy_id of the long read bam
    this_bam = op.filter_samples(BAM_SAMPLES, sample_id = w.sample_id, genome_build = w.genome_build, seq_type = w.seq_type)
    this_bx = this_bam["biopsy_id"].tolist()
    this_vcf = op.filter_samples(VCF_SAMPLES, biopsy_id = this_bx[0], seq_type = CFG["options"]["vcf_seq_type"])
    this_vcf = expand(
        [
            rules._whatshap_split_vcf.output.vcf
        ], 
        zip, 
        sample_id = this_vcf["sample_id"], 
        genome_build = this_vcf["genome_build"], 
        seq_type = CFG["options"]["vcf_seq_type"], 
        allow_missing = True
    )
    if not this_vcf:
        raise ValueError(f"No VCF files found for the given biopsy_id {this_bx} and seq_type {CFG['options']['vcf_seq_type']}.")
    return {"vcf": this_vcf, "tbi": [f"{vcf}.tbi" for vcf in this_vcf]}

rule _whatshap_phase_vcf:
    input:
        unpack(get_vcf),
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = temp(CFG["dirs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.phased.vcf.gz"),
        index = temp(CFG["dirs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.phased.vcf.gz.tbi")
    params: 
        options = CFG["options"]["phase_vcf"]
    conda:
        CFG["conda_envs"]["whatshap"] 
    resources: 
        **CFG["resources"]["phase_vcf"] 
    threads: CFG["threads"]["phase_vcf"]       
    wildcard_constraints: 
        sample_id = "|".join(BAM_SAMPLES["sample_id"].unique())
    group: "split_and_phase"
    log:
        CFG["logs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.phase_vcf.stderr.log"    
    shell:
        op.as_one_line(""" 
            echo "Running whatshap phase on $HOSTNAME" > {log} &&
            whatshap phase 
                --ignore-read-groups 
                --reference {input.fasta} 
                --chromosome {wildcards.chrom} 
                --output {output.vcf} 
                {params.options} 
                {input.vcf} 
                {input.bam} >> {log} 2>&1 && 
            tabix -p vcf {output.vcf}
        """)

def _whatshap_get_chr_vcf(wildcards):
    CFG = config["lcr-modules"]["whatshap"]
    chrs = checkpoints._whatshap_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    vcfs = expand(
        rules._whatshap_phase_vcf.output.vcf,
        chrom = chrs, 
        allow_missing = True
    )
    return vcfs

# Have to make redundant functions for this because 
# unpacking the input function doesn't work if the 
# checkpoint rule hasn't run yet.
def _whatshap_get_chr_tbi(wildcards):
    CFG = config["lcr-modules"]["whatshap"]
    chrs = checkpoints._whatshap_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    tbis = expand(
        rules._whatshap_phase_vcf.output.index,
        chrom = chrs, 
        allow_missing = True
    )
    return tbis

rule _whatshap_merge_vcf: 
    input: 
        vcfs = _whatshap_get_chr_vcf, 
        tbis = _whatshap_get_chr_tbi
    output: 
        vcf = CFG["dirs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz",
        index = CFG["dirs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}/whatshap_merge_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["merge_vcf"]
    resources:
        **CFG["resources"]["merge_vcf"]
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    wildcard_constraints: 
        sample_id = "|".join(BAM_SAMPLES["sample_id"].unique())
    shell:
        op.as_one_line("""
        bcftools concat --threads {threads} -a -O z {input.vcfs} 2> {log.stderr}
            |
        bcftools sort -m {params.mem_mb}M -O z -o {output.vcf} 2>> {log.stderr} 
            &&
        bcftools index -t --threads {threads} {output.vcf} 2>> {log.stderr}
        """)

rule _whatshap_stats:
    input:
        vcf = str(rules._whatshap_merge_vcf.output.vcf)
    output:
        gtf = CFG["dirs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.phased.blocks.gtf", 
        blocks = CFG["dirs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.phased.block_list.tsv",
        stats = CFG["dirs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}.phased.stats.tsv"
    params: 
        options = CFG["options"]["stats"]
    conda:
        CFG["conda_envs"]["whatshap"] 
    resources: 
        **CFG["resources"]["stats"]  
    threads: CFG["threads"]["stats"]      
    log: CFG["logs"]["phase_vcf"] + "{seq_type}--{genome_build}/{sample_id}/whatshap_phase_vcf.stderr.log"    
    shell:
        op.as_one_line(""" 
            whatshap stats 
                --block-list={output.blocks} 
                --tsv={output.stats} 
                --gtf={output.gtf} 
                {input.vcf} 
                > {log} 2>&1
        """)        

rule _whatshap_phase_bam:
    input:
        vcf = str(rules._whatshap_merge_vcf.output.vcf),
        bam = str(rules._whatshap_input_bam.output.bam),
        regions = str(rules._whatshap_input_regions.output.regions),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["phase_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.phased.bam"),
        bai = temp(CFG["dirs"]["phase_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.phased.bam.bai"), 
        haptag = CFG["dirs"]["phase_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.haplotag.txt"
    conda:
        CFG["conda_envs"]["whatshap"]
    resources: 
        **CFG["resources"]["phase_bam"]
    threads:
        CFG["threads"]["phase_bam"]           
    log:
        stderr = CFG["logs"]["phase_bam"] + "{seq_type}--{genome_build}/{sample_id}/{regions_bed}.whatshap_phase_bam.stderr.log"  
    shell:
        op.as_one_line(""" 
            whatshap haplotag 
                --output-threads={threads} 
                --output={output.bam}
                --reference={input.fasta} 
                --output-haplotag-list={output.haptag} 
                --regions "$(head -n1 {input.regions})"
                --ignore-read-groups 
                {input.vcf} 
                {input.bam}  
            && samtools index -@{threads} {output.bam} 2> {log.stderr}
        """)
        
rule _whatshap_cram_phased: 
    input: 
        bam = rules._whatshap_phase_bam.output.bam,
        bai = rules._whatshap_phase_bam.output.bai,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: 
        cram = CFG["dirs"]["phase_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.phased.cram",
        crai = CFG["dirs"]["phase_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.phased.cram.crai"
    conda: 
        CFG["conda_envs"]["whatshap"]
    resources: 
        **CFG["resources"]["cram"]
    threads: CFG["threads"]["cram"]
    shell: 
        op.as_one_line("""
            samtools view -@ {threads} -h -C -T {input.fasta} -o {output.cram} {input.bam} && 
            samtools index -@ {threads} {output.cram}
        """)

rule _whatshap_split_bam:
    input:
        bam = str(rules._whatshap_phase_bam.output.bam), 
        bai = str(rules._whatshap_phase_bam.output.bai),
        haptag = str(rules._whatshap_phase_bam.output.haptag)
    output:
        h1 = temp(CFG["dirs"]["split_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.H1.bam"),
        h2 = temp(CFG["dirs"]["split_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.H2.bam"), 
        unphased = temp(CFG["dirs"]["split_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.unphased.bam"), 
        hist = CFG["dirs"]["split_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.read_length_histogram.tsv"
    params: 
        options = CFG["options"]["split_bam"]
    conda:
        CFG["conda_envs"]["whatshap"]
    resources: 
        **CFG["resources"]["split_bam"] 
    threads: CFG["threads"]["split_bam"]
    log:
        stderr = CFG["logs"]["split_bam"] + "{seq_type}--{genome_build}/{sample_id}/{regions_bed}.whatshap_split_bam.stderr.log"  
    shell:
        op.as_one_line(""" 
            whatshap split 
                --output-h1 {output.h1} 
                --output-h2 {output.h2} 
                --output-untagged {output.unphased} 
                --read-lengths-histogram {output.hist}
                {params.options}
                {input.bam} 
                {input.haptag} 
                2> {log.stderr}
        """)

rule _whatshap_cram_split: 
    input: 
        bam = CFG["dirs"]["split_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.{haplotype}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: 
        cram = CFG["dirs"]["split_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.{haplotype}.cram", 
        crai = CFG["dirs"]["split_bam"] + "{seq_type}--{genome_build}/{sample_id}.{regions_bed}.{haplotype}.cram.crai"
    conda: 
        CFG["conda_envs"]["whatshap"]
    resources: 
        **CFG["resources"]["cram"]
    threads: CFG["threads"]["cram"]
    shell: 
        op.as_one_line("""
            samtools view -@ {threads} -h -C -T {input.fasta} -o {output.cram} {input.bam} && 
            samtools index -@ {threads} {output.cram}
        """)
    

# Symlinks the final output files into the module results directory (under '99-outputs/')

rule _whatshap_output_vcf:
    input:
        vcf = str(rules._whatshap_merge_vcf.output.vcf),
        index = str(rules._whatshap_merge_vcf.output.index),
        gtf = str(rules._whatshap_stats.output.gtf), 
        blocks = str(rules._whatshap_stats.output.blocks),
        stats = str(rules._whatshap_stats.output.stats)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz",
        index = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.phased.vcf.gz.tbi",
        gtf = CFG["dirs"]["outputs"] + "stats/{seq_type}--{genome_build}/{sample_id}.phased.blocks.gtf", 
        blocks = CFG["dirs"]["outputs"] + "stats/{seq_type}--{genome_build}/{sample_id}.phased.block_list.tsv",
        stats = CFG["dirs"]["outputs"] + "stats/{seq_type}--{genome_build}/{sample_id}.phased.stats.tsv"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True),
        op.relative_symlink(input.index, output.index, in_module= True),
        op.relative_symlink(input.gtf, output.gtf, in_module= True),
        op.relative_symlink(input.blocks, output.blocks, in_module= True),
        op.relative_symlink(input.stats, output.stats, in_module= True)

rule _whatshap_output_split:
    input:
        cram = str(rules._whatshap_cram_split.output.cram),
        crai = str(rules._whatshap_cram_split.output.crai)
    output:
        cram = CFG["dirs"]["outputs"] + "cram/{seq_type}--{genome_build}/{sample_id}.{regions_bed}.{haplotype}.split.cram",
        crai = CFG["dirs"]["outputs"] + "cram/{seq_type}--{genome_build}/{sample_id}.{regions_bed}.{haplotype}.split.cram.crai"
    wildcard_constraints: 
        regions_bed = "|".join(CFG["inputs"]["regions_bed"].keys()), 
        haplotag = "|".join(["H1", "H2", "unphased"])
    run:
        op.relative_symlink(input.cram, output.cram, in_module= True)
        op.relative_symlink(input.crai, output.crai, in_module= True)

rule _whatshap_output_phased:
    input:
        cram = str(rules._whatshap_cram_phased.output.cram),
        crai = str(rules._whatshap_cram_phased.output.crai)
    output:
        cram = CFG["dirs"]["outputs"] + "cram/{seq_type}--{genome_build}/{sample_id}.{regions_bed}.phased.cram",
        crai = CFG["dirs"]["outputs"] + "cram/{seq_type}--{genome_build}/{sample_id}.{regions_bed}.phased.cram.crai"
    wildcard_constraints: 
        regions_bed = "|".join(CFG["inputs"]["regions_bed"].keys())
    run:
        op.relative_symlink(input.cram, output.cram, in_module= True)
        op.relative_symlink(input.crai, output.crai, in_module= True)

# Generates the target sentinels for each run, which generate the symlinks

# Specify this as the target rule if you only want phased vcf and associated stats
rule _whatshap_vcf_all: 
    input: 
        expand(
            rules._whatshap_output_vcf.output, 
            zip, 
            seq_type=BAM_SAMPLES["seq_type"],
            genome_build=BAM_SAMPLES["genome_build"],
            sample_id=BAM_SAMPLES["sample_id"]
        )

rule _whatshap_all:
    input:
        rules._whatshap_vcf_all.input, 
        expand(
            expand(
                rules._whatshap_output_phased.output,
                zip,  
                seq_type=BAM_SAMPLES["seq_type"],
                genome_build=BAM_SAMPLES["genome_build"],
                sample_id=BAM_SAMPLES["sample_id"], 
                allow_missing = True
            ), 
            regions_bed = CFG["inputs"]["regions_bed"].keys()
        ), 
        expand(
            expand(
                rules._whatshap_output_split.output,
                zip,  
                seq_type=BAM_SAMPLES["seq_type"],
                genome_build=BAM_SAMPLES["genome_build"],
                sample_id=BAM_SAMPLES["sample_id"], 
                allow_missing = True
            ), 
            haplotype = ["H1", "H2", "unphased"], 
            regions_bed = CFG["inputs"]["regions_bed"].keys()
        )



##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)  
