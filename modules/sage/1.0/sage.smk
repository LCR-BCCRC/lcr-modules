#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.12"
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
# `CFG` is a shortcut to `config["lcr-modules"]["sage"]`
CFG = op.setup_module(
    name = "sage",
    version = "1.0",
    subdirectories = ["inputs", "sage", "vcf", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _sage_input_bam,
    _input_references,
    _download_sage_references,
    _sage_filter_vcf,
    _sage_split_vcf,
    _sage_output_vcf,
    _sage_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _sage_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bam+ ".bai", output.bai)


# Setup shared reference files. Symlinking these files to 00-inputs to ensure index and dictionary are present
# before pipeline starts, otherwise on a fresh directory dictionary is not created and worflow exits.
rule _input_references:
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

# Download specific reference files
rule _download_sage_references:
    input:
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        hotspots = CFG["dirs"]["inputs"] + "references/{genome_build}/sage/KnownHotspots.vcf.gz",
        panel_bed = CFG["dirs"]["inputs"] + "references/{genome_build}/sage/ActionableCodingPanel.somatic.bed.gz",
        high_conf_bed = CFG["dirs"]["inputs"] + "references/{genome_build}/sage/HighConfidence.bed.gz",
        cache = directory(CFG["dirs"]["inputs"] + "references/{genome_build}/ensembl_cache/"),
        complete = touch(CFG["dirs"]["inputs"] + "references/{genome_build}/ensembl_cache/cache.complete")
    params:
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/sage",
        url_cache = "www.bcgsc.ca/downloads/morinlab/hmftools-references/ensembl_data_cache",
        build = lambda w: "hg38" if "38" in str({w.genome_build}) else "hg19",
        cache_build = lambda w: "38" if "38" in str({w.genome_build}) else "37"
    conda:
        CFG["conda_envs"]["wget"]
    shell:
        op.as_one_line("""
        wget -O {output.hotspots} {params.url}/KnownHotspots.{params.build}.vcf.gz
          &&
        wget -O {output.panel_bed} {params.url}/ActionableCodingPanel.somatic.{params.build}.bed.gz
          &&
        wget -O {output.high_conf_bed} {params.url}/HighConfidence.{params.build}.bed.gz
          &&
        wget -O {output.cache}/{params.cache_build}.zip {params.url_cache}/{params.cache_build}.zip
          &&
        unzip -d {output.cache} {output.cache}/{params.cache_build}.zip
          &&
        touch {output.complete}
        """)

def _sage_get_capspace(wildcards):
    CFG = config["lcr-modules"]["sage"]
    try:
        # Get the appropriate capture space for this sample
        this_bed = op.get_capture_space(CFG, wildcards.tumour_id, wildcards.genome_build, wildcards.seq_type, "bed.gz")
        this_bed = reference_files(this_bed)
    except NameError:
        # If we are using an older version of the reference workflow, use the same region file as the genome sample
        this_bed = rules._download_sage_references.output.panel_bed
    # If this is a genome sample, return a BED file listing all chromosomes
    if wildcards.seq_type != "capture":
        this_bed = rules._download_sage_references.output.panel_bed
    return this_bed

# Variant calling rule
rule _run_sage:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = str(rules._input_references.output.genome_fa),
        hotspots = rules._download_sage_references.output.hotspots,
        high_conf_bed = str(rules._download_sage_references.output.high_conf_bed),
        cache = str(rules._download_sage_references.output.cache),
        panel_bed = _sage_get_capspace,
        main_chromosomes = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        vcf = temp(CFG["dirs"]["sage"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}--{normal_id}--{pair_status}.vcf"),
        vcf_gz = temp(CFG["dirs"]["sage"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}--{normal_id}--{pair_status}.vcf.gz")
    log:
        stdout = CFG["logs"]["sage"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/run_sage.stdout.log",
        stderr = CFG["logs"]["sage"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/run_sage.stderr.log"
    params:
        opts = CFG["options"]["sage_run"],
        assembly = lambda w: "38" if "38" in str({w.genome_build}) else "37",
        cache_dir = lambda w: config["lcr-modules"]["sage"]["dirs"]["inputs"] + "references/" + w.genome_build + "/ensembl_cache/" + str("38" if "38" in str({w.genome_build}) else "37"),
        sage= "$(dirname $(readlink -e $(which SAGE)))/sage.jar",
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["sage"]
    threads:
        CFG["threads"]["sage_run"]
    resources:
        **CFG["resources"]["sage_run"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname)" > {log.stdout}
        &&
        SAGE_CHROMOSOMES=$(cat {input.main_chromosomes} | paste -sd ";" -)
        &&
        java -Xms1G -Xmx{params.jvmheap}m
        -cp {params.sage} com.hartwig.hmftools.sage.SageApplication
        -threads {threads}
        {params.opts}
        -specific_chr $SAGE_CHROMOSOMES
        -reference {wildcards.normal_id}
        -reference_bam {input.normal_bam}
        -tumor {wildcards.tumour_id}
        -tumor_bam {input.tumour_bam}
        -ref_genome_version {params.assembly}
        -ref_genome {input.fasta}
        -hotspots {input.hotspots}
        -panel_bed {input.panel_bed}
        -high_confidence_bed {input.high_conf_bed}
        -ensembl_data_dir {params.cache_dir}
        -out {output.vcf}
        >>  {log.stdout} 2>> {log.stderr}
        && bgzip -c {output.vcf} > {output.vcf_gz}
        """)

# Filter resulting VCF file on PASS variants
rule _sage_filter_vcf:
    input:
        vcf = str(rules._run_sage.output.vcf_gz)
    output:
        vcf_passed = CFG["dirs"]["vcf"] + "combined/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/passed_combined.vcf.gz",
        vcf_passed_tbi = CFG["dirs"]["vcf"] + "combined/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/passed_combined.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/filter_passed.stdout.log",
        stderr = CFG["logs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/filter_passed.stderr.log"
    params:
        heap_mem = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["filter"]
    resources:
        **CFG["resources"]["filter"]
    shell:
        op.as_one_line("""
        bcftools view -f ".,PASS" {input.vcf} -Ov
          |
        bcftools sort --max-mem {params.heap_mem}M -Oz -o {output.vcf_passed}
        >> {log.stdout} 2>> {log.stderr}
          &&
        tabix -p vcf {output.vcf_passed} >> {log.stdout} 2>> {log.stderr}
        """)

# Split filtered VCF file into snvs and indels
rule _sage_split_vcf:
    input:
        vcf_passed = str(rules._sage_filter_vcf.output.vcf_passed)
    output:
        indels = CFG["dirs"]["vcf"] + "indels/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/indels.vcf.gz",
        indels_tbi = CFG["dirs"]["vcf"] + "indels/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/indels.vcf.gz.tbi",
        snvs = CFG["dirs"]["vcf"] + "snvs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/snvs.vcf.gz",
        snvs_tbi = CFG["dirs"]["vcf"] + "snvs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/snvs.vcf.gz.tbi",
    log:
        stdout = CFG["logs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/split_passed.stdout.log",
        stderr = CFG["logs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/split_passed.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["filter"]
    resources:
        **CFG["resources"]["filter"]
    shell:
        op.as_one_line("""
        bcftools view -v indels {input.vcf_passed} -Oz -o {output.indels} >> {log.stdout} 2>> {log.stderr}
          &&
        tabix -p vcf {output.indels} >> {log.stdout} 2>> {log.stderr}
          &&
        bcftools view -v snps {input.vcf_passed} -Oz -o {output.snvs} >> {log.stdout} 2>> {log.stderr}
          &&
        tabix -p vcf {output.snvs} >> {log.stdout} 2>> {log.stderr}
        """)



# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sage_output_vcf:
    input:
        combined = str(rules._sage_filter_vcf.output.vcf_passed),
        indels = str(rules._sage_split_vcf.output.indels),
        snvs = str(rules._sage_split_vcf.output.snvs),
    output:
        combined = CFG["dirs"]["outputs"] + "combined/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz",
        combined_tbi = CFG["dirs"]["outputs"] + "combined/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz.tbi",
        indels = CFG["dirs"]["outputs"] + "indels/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.indels.vcf.gz",
        indels_tbi = CFG["dirs"]["outputs"] + "indels/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.indels.vcf.gz.tbi",
        snvs = CFG["dirs"]["outputs"] + "snvs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.snvs.vcf.gz",
        snvs_tbi = CFG["dirs"]["outputs"] + "snvs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.snvs.vcf.gz.tbi"
    run:
        op.relative_symlink(input.combined, output.combined,in_module=True)
        op.relative_symlink(input.combined+".tbi", output.combined+".tbi",in_module=True)
        op.relative_symlink(input.indels, output.indels,in_module=True)
        op.relative_symlink(input.indels+".tbi", output.indels+".tbi",in_module=True)
        op.relative_symlink(input.snvs, output.snvs,in_module=True)
        op.relative_symlink(input.snvs+".tbi", output.snvs+".tbi",in_module=True)


# Generates the targets each run
rule _sage_all:
    input:
        expand(
            [
            str(rules._sage_output_vcf.output.combined),
            str(rules._sage_output_vcf.output.indels),
            str(rules._sage_output_vcf.output.snvs),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
