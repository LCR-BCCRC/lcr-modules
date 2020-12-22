#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["sage"]`
CFG = op.setup_module(
    name = "sage",
    version = "1.0",
    subdirectories = ["inputs", "sage", "vcf", "outputs"]
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _sage_input_bam,
    _run_sage,
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
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bam+ ".bai", output.bai)


# Setup shared reference files
rule _input_references: 
    input: 
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        genome_fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        genome_dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output: 
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.fa", 
        genome_fai = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.fa.fai", 
        genome_dict = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.dict"
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
        high_conf_bed = CFG["dirs"]["inputs"] + "references/{genome_build}/sage/HighConfidence.bed.gz"
    params:
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/sage",
        build = lambda w: "hg38" if "38" in str({w.genome_build}) else "hg19"
    shell: 
        op.as_one_line("""
        wget -O {output.hotspots} {params.url}/KnownHotspots.{params.build}.vcf.gz
          &&
        wget -O {output.panel_bed} {params.url}/ActionableCodingPanel.somatic.{params.build}.bed.gz
          &&
        wget -O {output.high_conf_bed} {params.url}/HighConfidence.{params.build}.bed.gz
        """)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _run_sage:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = str(rules._input_references.output.genome_fa),
        hotspots = str(rules._download_sage_references.output.hotspots),
        panel_bed = str(rules._download_sage_references.output.panel_bed),
        high_conf_bed = str(rules._download_sage_references.output.high_conf_bed)
    output:
        vcf = temp(CFG["dirs"]["sage"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}--{normal_id}--{pair_status}.vcf"),
        vcf_gz = CFG["dirs"]["sage"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}--{normal_id}--{pair_status}.vcf.gz"
    log:
        stdout = CFG["logs"]["sage"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/run_sage.stdout.log",
        stderr = CFG["logs"]["sage"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/run_sage.stderr.log"
    params:
        opts = CFG["options"]["sage_run"],
        assembly = lambda w: "hg38" if "38" in str({w.genome_build}) else "hg19",
        sage= "$(dirname $(readlink -e $(which SAGE)))/sage.jar"
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
        java -Xms4G -Xmx{resources.mem_mb}
        -cp {params.sage} com.hartwig.hmftools.sage.SageApplication
        -threads {threads}
        {params.opts}
        -reference {wildcards.normal_id}
        -reference_bam {input.normal_bam}
        -tumor {wildcards.tumour_id} 
        -tumor_bam {input.tumour_bam}
        -assembly {params.assembly}
        -ref_genome {input.fasta}
        -hotspots {input.hotspots}
        -panel_bed {input.panel_bed}
        -high_confidence_bed {input.high_conf_bed}
        -out {output.vcf}
        >>  {log.stdout} 2> {log.stderr}
        && bgzip {output.vcf}
        """)


rule _sage_filter_vcf:
    input:
        vcf = str(rules._run_sage.output.vcf_gz)
    output:
        vcf_passed = CFG["dirs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/passed_output.vcf.gz",
        vcf_passed_tbi = CFG["dirs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/passed_output.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/filter_passed.stdout.log",
        stderr = CFG["logs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/filter_passed.stderr.log"
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
        bcftools sort --max-mem {resources.mem_mb}M -Oz -o {output.vcf_passed}
        > {log.stdout} 2> {log.stderr}
          &&
        tabix -p vcf {output.vcf_passed} >> {log.stdout} 2>> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sage_output_vcf:
    input:
        vcf = str(rules._sage_filter_vcf.output.vcf_passed)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz",
        vcf_tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.vcf+".tbi", output.vcf+".tbi")


# Generates the target sentinels for each run, which generate the symlinks
rule _sage_all:
    input:
        expand(
            [
                str(rules._sage_output_vcf.output.vcf),
                str(rules._sage_output_vcf.output.vcf_tbi)
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
