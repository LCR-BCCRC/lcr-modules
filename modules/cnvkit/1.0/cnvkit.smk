#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Jasper Wong
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["cnvkit"]`
CFG = op.setup_module(
    name = "cnvkit",
    version = "1.0",
    subdirectories = ["inputs", "coverage", "fix", "cns", "SNPs", "call", "plots", "breaks", "gene_metrics", "seg", "convert_coordinates", "fill_regions", "normalize", "metrics", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _cnvkit_input_bam,
    _cnvkit_input_chroms,
    _cnvkit_symlink_beds,
    _cnvkit_symlink_pon_reference,
    _cnvkit_output,
    _cnvkit_output_projection,
    _cnvkit_all


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _cnvkit_input_bam:
    input:
        bam = ancient(CFG["inputs"]["sample_bam"])
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)

# Pulls in list of chromosomes for the genome builds
checkpoint _cnvkit_input_chroms:
    input:
        txt = ancient(reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt"))
    output:
        txt = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes_withY.txt"
    run:
        op.absolute_symlink(input.txt, output.txt)

# Recreate index so the timestamp will always be later then the bam
rule _cnvkit_index_bam:
    input:
        bam = str(rules._cnvkit_input_bam.output.bam)
    output:
        bai = temp(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.bam.bai")
    log:
        log = CFG["logs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_index.log"
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["samtools"]
    shell:
        op.as_one_line("""
        samtools index -@ {threads} {input.bam} 2> {log.log} &&
        cd $(dirname {input.bam});
        if [[ -e {wildcards.tumour_id}.bam.crai ]];
        then
            ln -s {wildcards.tumour_id}.bam.crai {wildcards.tumour_id}.bam.bai;
        fi
        """)

rule _cnvkit_symlink_beds:
    input:
        target = CFG["inputs"]["target_bed"],
        antitarget = CFG["inputs"]["antitarget_bed"]
    output:
        target = CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/{capture_space}_target_sites.bed",
        antitarget = CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/{capture_space}_antitarget_sites.bed"
    run:
        op.relative_symlink(input.target, output.target)
        op.relative_symlink(input.antitarget, output.antitarget)

rule _cnvkit_symlink_pon_reference:
    input:
        pon =  CFG["inputs"]["pon_reference"]
    output:
        pon =  CFG["dirs"]["inputs"] + "pon/{seq_type}--{genome_build}/{capture_space}_normal_reference.cnn"
    run:
        op.relative_symlink(input.pon, output.pon)

# Coverage for each sample
rule _cnvkit_coverage_target:
    input:
        bam = str(rules._cnvkit_input_bam.output.bam),
        bai = str(rules._cnvkit_index_bam.output.bai),
        bed = str(rules._cnvkit_symlink_beds.output.target),
    output:
        cov = CFG["dirs"]["coverage"] + "target/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.targetcoverage.cnn"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["coverage"]
    resources:
        **CFG["resources"]["coverage"]
    log:
        log = CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_target.log"
    shell:
        op.as_one_line("""
        cnvkit.py coverage {input.bam} {input.bed} -o {output.cov} -p {threads}
         &> {log.log}
        """)

rule _cnvkit_coverage_antitarget:
    input:
        bam = str(rules._cnvkit_input_bam.output.bam),
        bai = str(rules._cnvkit_index_bam.output.bai),
        bed = str(rules._cnvkit_symlink_beds.output.antitarget),
    output:
        cov = CFG["dirs"]["coverage"] + "antitarget/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.antitargetcoverage.cnn"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["coverage"]
    resources:
        **CFG["resources"]["coverage"]
    log:
        log = CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_antitarget.log"
    shell:
        op.as_one_line("""
        cnvkit.py coverage {input.bam} {input.bed} -o {output.cov} -p {threads}
         &> {log.log}
        """)

# Fix coverage using panel_of_normals reference cnn
rule _cnvkit_fix:
    input:
        targetcov = str(rules._cnvkit_coverage_target.output.cov),
        antitargetcov = str(rules._cnvkit_coverage_antitarget.output.cov),
        pon_reference = str(rules._cnvkit_symlink_pon_reference.output.pon)
    output:
        cnr = CFG["dirs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["fix"]
    log:
        stdout = CFG["logs"]["fix"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.log"
    shell:
        op.as_one_line("""
            cnvkit.py fix {input.targetcov} {input.antitargetcov} {input.pon_reference} -o {output.cnr} &> {log.stdout}
        """)

rule _cnvkit_segment:
    input:
        cnr = str(rules._cnvkit_fix.output.cnr)
    output:
        cns = CFG["dirs"]["cns"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cns"
    params:
        method = CFG["options"]["cns"]["method"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["cns"]
    resources:
        **CFG["resources"]["cns"]
    log:
        log = CFG["logs"]["cns"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.segment.log"
    shell:
        op.as_one_line("""
            cnvkit.py segment {input.cnr} -o {output.cns} -p {threads} --drop-low-coverage -m {params.method} &> {log.log}
        """)

# need SNPs not SNVs (i.e. get germline calls using a dbSNP vcf)
rule _cnvkit_dbsnp_to_bed:
    input:
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        bed = CFG["dirs"]["inputs"] + "dbsnp/{genome_build}/dbsnp.common_all-151.bed"
    log:
        stderr = CFG["logs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/dbsnp_to_bed.log"
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["SNPs"]
    shell:
        op.as_one_line("""
        gunzip -c {input.vcf} | awk {{'printf ("%s\\t%s\\t%s\\n", $1,$2-1,$2)'}} | zgrep -v -h "^#" > {output.bed} 2> {log.stderr}
        """)

# vcf needs DP, GT, AD - bcftools -mv calls multiallelic variants (will annotate GT)
# without it, GT will not be annotated
rule _cnvkit_mpileup_per_chrom:
    input:
        bam = str(rules._cnvkit_input_bam.output.bam),
        bai = str(rules._cnvkit_index_bam.output.bai),
        fastaFile = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        bed = str(rules._cnvkit_dbsnp_to_bed.output.bed)
    output: # creates a temporary file for mpileup
        vcf = temp(CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.{chrom}.vcf.gz.tbi")
    params:
        quality = CFG["options"]["SNPs"]["quality"],
        opts = CFG["options"]["SNPs"]["opts"]
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["SNPs"]
    resources:
        **CFG["resources"]["SNPs"]
    group: "cnvkit"
    log:
        stderr = CFG["logs"]["SNPs"] + "{capture_space}/{seq_type}--{genome_build}/{tumour_id}/{chrom}.vcf.stderr.log"
    shell:
         op.as_one_line("""
        bcftools mpileup
         --threads {threads}
          -T {input.bed}
          -r {wildcards.chrom}
          -f {input.fastaFile}
          -Q {params.quality}
          {params.opts}
          -Ou {input.bam} |
         bcftools call -mv -Oz --threads {threads} -o {output.vcf} 2> {log.stderr} &&
         tabix -@ {threads} -p vcf {output.vcf} 2>> {log.stderr}
        """)

# Collect vcfs
def _cnvkit_get_chr_mpileups(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    with open(checkpoints._cnvkit_input_chroms.get(**wildcards).output.txt) as f:
        mains_chroms = f.read().rstrip("\n").split("\n")
    vcfs = expand(
        CFG["dirs"]["SNPs"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{tumour_id}}.{chrom}.vcf.gz",
        chrom = chrs
    )
    tbis = expand(
        CFG["dirs"]["SNPs"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{tumour_id}}.{chrom}.vcf.gz.tbi",
        chrom = chrs
    )
    return {
        "vcf" = vcfs,
        "tbi" = tbis
    }

rule _cnvkit_concatenate_vcf:
    input:
        unpack(_cnvkit_get_chr_mpileups)
    output:
        vcf = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.vcf.gz",
        tbi = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["SNPs"] + "{capture_space}/{seq_type}--{genome_build}/{tumour_id}.vcf.stderr.log"
    threads:
        CFG["threads"]["SNPs"]
    resources:
        **CFG["resources"]["SNPs"]
    group: "cnvkit"
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools concat {input.vcf} -Oz --threads {threads} -o {output.vcf} 2> {log.stderr} &&
         tabix -p vcf -@ {threads} {output.vcf} 2>> {log.stderr}
        """)


# ----------------------------------------------------------------------------------------------- #
# Integrating cnvkit with BAF to call absolute CN
# ----------------------------------------------------------------------------------------------- #
# Adds extra columns - One-sample t-test of bin log2 ratios versus 0.0 and ci high and ci low to be able to use filtering by ci in the next step
# Note that the t-test is not used in filtration step, but the ci is
rule _cnvkit_segmetrics_ttest:
    input:
        cnr = str(rules._cnvkit_fix.output.cnr),
        cns = str(rules._cnvkit_segment.output.cns)
    output:
        cns = temp(CFG["dirs"]["call"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.temp.cns")
    params:
        add_col = CFG["options"]["segmetrics"]["add_col"]
    log:
        log = CFG["logs"]["call"] + "segmetrics/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_segmetrics.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["call"]
    shell:
        op.as_one_line("""
        cnvkit.py segmetrics {input.cnr} -s {input.cns} -o {output.cns}
         {params.add_col} &> {log.log}
        """)


rule _cnvkit_call:
    input:
        cns = str(rules._cnvkit_segmetrics_ttest.output.cns),
        vcf = str(rules._cnvkit_concatenate_vcf.output.vcf),
        tbi = str(rules._cnvkit_concatenate_vcf.output.tbi)
    output:
        cns =  CFG["dirs"]["call"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cns"
    params:
        rescale = CFG["options"]["call"]["rescale"],
        min_depth = CFG["options"]["call"]["min_depth"],
        filter_by = CFG["options"]["call"]["filter_by"],
        male_ref = CFG["options"]["male_ref"],
        opts = CFG["options"]["call"]["opts"]
    log:
        log = CFG["logs"]["call"] + "call/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_call.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["call"]
    shell:
        op.as_one_line("""
        cnvkit.py call {input.cns} --output {output.cns} -v {input.vcf}
         --min-variant-depth {params.min_depth} -m {params.rescale} --filter
          {params.filter_by} {params.male_ref} {params.opts} &> {log.log}
        """)


# plot a scatter plot of amps and dels, also BAF
rule _cnvkit_scatter:
    input:
        cnr = str(rules._cnvkit_fix.output.cnr),
        cns = str(rules._cnvkit_call.output.cns),
        vcf = str(rules._cnvkit_concatenate_vcf.output.vcf),
        tbi = str(rules._cnvkit_concatenate_vcf.output.tbi)
    output:
        pdf = CFG["dirs"]["plots"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_scatter.pdf"
    params:
        min_depth = CFG["options"]["scatter"]["min_depth"],
        ymax = CFG["options"]["scatter"]["ymax"],
        ymin = CFG["options"]["scatter"]["ymin"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["plots"]
    log:
        log = CFG["logs"]["plots"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_scatter.log"
    shell:
        op.as_one_line("""
        cnvkit.py scatter {input.cnr} -s {input.cns} --output
         {output.pdf} -v {input.vcf} --min-variant-depth {params.min_depth} --y-max
         {params.ymax} --y-min {params.ymin} &> {log.log}
        """)


# plot chromosome diagrams highlighting these amps/dels and also key genes that are located in these CNVs
rule _cnvkit_diagram:
    input:
        cnr = str(rules._cnvkit_fix.output.cnr),
        cns = str(rules._cnvkit_call.output.cns)
    output:  # only pdf works
        pdf = CFG["dirs"]["plots"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_diagram.pdf"
    params:
        threshold = CFG["options"]["diagram"]["threshold"], # to only label genes in high level amps and dels
        male_ref = CFG["options"]["male_ref"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["plots"]
    log:
        log = CFG["logs"]["plots"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_diagram.log"
    shell:
        op.as_one_line("""
        cnvkit.py diagram {input.cnr} -s {input.cns} --output {output.pdf}
         -t {params.threshold} {params.male_ref} &> {log.log}
        """)


# find potential breakpoints across the CNVs or regions with large CN signal inconsistencies
rule _cnvkit_breaks:
    input:
        cnr = str(rules._cnvkit_fix.output.cnr),
        cns = str(rules._cnvkit_call.output.cns)
    output:
        breaks = CFG["dirs"]["breaks"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.genebreaks.txt"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["breaks"]
    log:
        log = CFG["logs"]["breaks"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_breaks.log"
    shell:
        op.as_one_line("""
        cnvkit.py breaks {input.cnr} {input.cns} > {output.breaks} 2> {log.log}
        """)


# with segments (cns) as input, the minimum probes option defines the segment's bin count
# without cns as input, the gene's weighted bin counts are used instead
rule _cnvkit_genemetrics_seg:
    input:
        cnr = str(rules._cnvkit_fix.output.cnr),
        cns = str(rules._cnvkit_call.output.cns)
    output:
        genemetrics = CFG["dirs"]["gene_metrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/segment.gene_cn.txt"
    log:
        log = CFG["logs"]["gene_metrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_genemetrics_seg.log"
    params:
        threshold = CFG["options"]["gene_metrics"]["threshold"],
        min_segments = CFG["options"]["gene_metrics"]["min_segments"],
        male_ref = CFG["options"]["male_ref"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["gene_metrics"]
    shell:
        op.as_one_line("""
        cnvkit.py genemetrics {input.cnr} -s {input.cns} --threshold
         {params.threshold} --min-probes {params.min_segments}
         {params.male_ref} > {output.genemetrics} 2> {log.log}
        """)

rule _cnvkit_genemetrics_gene:
    input:
        cnr = str(rules._cnvkit_fix.output.cnr)
    output:
        genemetrics = CFG["dirs"]["gene_metrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/bin.gene_cn.txt"
    log:
        log = CFG["logs"]["gene_metrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_genemetrics_gene.log"
    params:
        threshold = CFG["options"]["gene_metrics"]["threshold"],
        min_segments = CFG["options"]["gene_metrics"]["min_segments"], # to remove false positives that cover a small number of bins
        male_ref = CFG["options"]["male_ref"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["gene_metrics"]
    shell:
        op.as_one_line("""
        cnvkit.py genemetrics {input.cnr} --threshold
         {params.threshold} --min-probes {params.min_segments}
         {params.male_ref} > {output.genemetrics} 2> {log.log}
        """)


# can take the intersection of the two methods to filter for a list of genes that confidently have CN change
rule _cnvkit_trusted_genes_cna:
    input:
        genemetrics_seg = str(rules._cnvkit_genemetrics_seg.output.genemetrics),
        genemetrics_gene = str(rules._cnvkit_genemetrics_gene.output.genemetrics)
    output:
        trusted_genes = CFG["dirs"]["gene_metrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/trusted_genes.txt"
    threads: 1 # does not use parallelization, but still needs to be submitted
    shell:
        op.as_one_line("""
        comm -12 <(tail -n+2 {input.genemetrics_seg} | cut -f1 | sort ) <(tail -n+2 {input.genemetrics_gene} | cut -f1 | sort ) > {output.trusted_genes}
        """)


rule _cnvkit_infer_sex:
    input:
        targetcov = str(rules._cnvkit_coverage_target.output.cov),
        antitargetcov = str(rules._cnvkit_coverage_antitarget.output.cov),
        cnr = str(rules._cnvkit_fix.output.cnr),
        cns = str(rules._cnvkit_segment.output.cns),
        call = str(rules._cnvkit_call.output.cns)
    output:
        sex = CFG["dirs"]["gene_metrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/inferred_sex.txt"
    log:
        log = CFG["logs"]["gene_metrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_infer_sex.log"
    params:
        male_ref = CFG["options"]["male_ref"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["gene_metrics"]
    shell:
        op.as_one_line("""
        cnvkit.py sex {input.targetcov} {input.antitargetcov} {input.cnr} {input.cns} {input.call} {params.male_ref} > {output.sex} 2> {log.log}
        """)


rule _cnvkit_metrics:
    input:
        cnr = str(rules._cnvkit_fix.output.cnr),
        cns = str(rules._cnvkit_call.output.cns)
    output:
        metrics = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.metrics.txt"
    log:
        log = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_metrics.log"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads: 1 # does not use parallelization, but still needs to be submitted
    resources:
        **CFG["resources"]["gene_metrics"]
    shell:
        op.as_one_line("""
        cnvkit.py metrics {input.cnr} -s {input.cns} > {output.metrics} 2> {log.log}
        """)

rule _cnvkit_cnv2igv:
    input:
        cns = str(rules._cnvkit_call.output.cns),
        cnv2igv = ancient(CFG["inputs"]["cnv2igv"])
    output:
        seg = CFG["dirs"]["cnv2igv"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg"
    log:
        stderr = CFG["logs"]["cnv2igv"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_cnv2igv.stderr.log"
    conda:
        CFG["conda_envs"]["cnv2igv"]
    threads: 1
    shell:
        op.as_one_line("""
        python {input.cnv2igv} --mode cnvkit {params.opts} --sample {wildcards.tumour_id} {input.cns} > {output.seg} 2>> {log.stderr}
        """)


def _cnvkit_get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")


# Convert the coordinates of seg file to a different genome build
rule _cnvkit_convert_coordinates:
    input:
        cnvkit_native = str(rules._cnvkit_cnv2igv.output.seg),
        cnvkit_chain = _cnvkit_get_chain
    output:
        cnvkit_lifted = CFG["dirs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{capture_space}/{tumour_id}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}/{tumour_id}.lifted_{chain}.stderr.log"
    threads: 1
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        liftover_minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    shell:
        op.as_one_line("""
        bash {params.liftover_script}
        SEG
        {input.cnvkit_native}
        {output.cnvkit_lifted}
        {input.cnvkit_chain}
        YES
        {params.liftover_minmatch}
        2>> {log.stderr}
        """)


def _cnvkit_prepare_projection(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    tbl = CFG["runs"]
    this_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"].tolist()
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    if "38" in this_genome_build[0]:
        hg38_projection = str(rules._cnvkit_cnv2igv.output.seg).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        grch37_projection = str(rules._cnvkit_convert_coordinates.output.cnvkit_lifted).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        grch37_projection = grch37_projection.replace("{chain}", "hg38ToHg19")
    else:
        grch37_projection = str(rules._cnvkit_cnv2igv.output.seg).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        hg38_projection = str(rules._cnvkit_convert_coordinates.output.cnvkit_lifted).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        hg38_projection = hg38_projection.replace("{chain}", "hg19ToHg38")
    return{
        "grch37_projection": grch37_projection,
        "hg38_projection": hg38_projection
    }


# Fill the missing segments of seg files with neutral regions to complete the genome coverage
rule _cnvkit_fill_segments:
    input:
        unpack(_cnvkit_prepare_projection)
    output:
        grch37_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}.{tool}.grch37.seg"),
        hg38_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}.{tool}.hg38.seg")
    log:
        stderr = CFG["logs"]["fill_regions"] + "{seq_type}--projection/{tumour_id}.{tool}_fill_segments.stderr.log"
    threads: 1
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"]
    conda:
        CFG["conda_envs"]["bedtools"]
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

# Normalize chr prefix of the output file
rule _cnvkit_normalize_projection:
    input:
        filled = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{masked}/{tumour_id}.{tool}.{projection}.seg",
        chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
    output:
        projection = CFG["dirs"]["normalize"] + "seg/{seq_type}--projection/{tumour_id}.{tool}.{projection}.seg"
    resources:
        **CFG["resources"]["post_cnvkit"]
    threads: 1
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
            seg_open.to_csv(output.projection, sep="\t", index=False, na_rep='NA')
        else:
            # remove chr prefix
            seg_open = pd.read_csv(input.filled, sep = "\t")
            seg_open["chrom"] = seg_open["chrom"].astype(str).str.replace('chr', '')
            seg_open.to_csv(output.projection, sep="\t", index=False, na_rep='NA')

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _cnvkit_output_projection:
    input:
        projection = str(rules._cnvkit_normalize_projection.output.projection)
    output:
        projection = CFG["dirs"]["outputs"] + "seg/{seq_type}--projection/{masked}/{tumour_id}.{tool}.{projection}.seg"
    run:
        op.relative_symlink(input.projection, output.projection, in_module = True)


# The rest of the outputs will keep the capture_space wildcard
rule _cnvkit_output:
    input:
        cns = str(rules._cnvkit_call.output.cns),
        scatter = str(rules._cnvkit_scatter.output.pdf),
        diagram = str(rules._cnvkit_diagram.output.pdf),
        breaks = str(rules._cnvkit_breaks.output.breaks),
        genemetrics_seg = str(rules._cnvkit_genemetrics_seg.output.genemetrics),
        genemetrics_gene = str(rules._cnvkit_genemetrics_gene.output.genemetrics),
        trusted_genes = str(rules._cnvkit_trusted_genes_cna.output.trusted_genes),
        sex = str(rules._cnvkit_infer_sex.output.sex),
        metrics = str(rules._cnvkit_metrics.output.metrics),
        seg = str(rules._cnvkit_cnv2igv.output.seg)
    output:
        cns =  CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.cns",
        scatter = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_scatter.pdf",
        diagram = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_diagram.pdf",
        breaks = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_genebreaks.txt",
        genemetrics_seg = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.segment.gene_cn.txt",
        genemetrics_gene = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.bin.gene_cn.txt",
        sex = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_inferred_sex.txt",
        metrics = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.metrics.txt"
        seg = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.seg"
    run:
        op.relative_symlink(input.cns, output.cns, in_module = True)
        op.relative_symlink(input.scatter, output.scatter, in_module = True)
        op.relative_symlink(input.diagram, output.diagram, in_module = True)
        op.relative_symlink(input.breaks, output.breaks, in_module = True)
        op.relative_symlink(input.genemetrics_seg, output.genemetrics_seg, in_module = True)
        op.relative_symlink(input.genemetrics_gene, output.genemetrics_gene, in_module = True)
        op.relative_symlink(input.sex, output.sex, in_module = True)
        op.relative_symlink(input.metrics, output.metrics, in_module = True)
        op.relative_symlink(input.seg, output.seg, in_module = True)



# Generates the target sentinels for each run, which generate the symlinks
rule _cnvkit_all:
    input:
        expand(
            [
                str(rules._cnvkit_output.output.cns),
                str(rules._cnvkit_output.output.scatter),
                str(rules._cnvkit_output.output.diagram),
                str(rules._cnvkit_output.output.breaks),
                str(rules._cnvkit_output.output.genemetrics_seg),
                str(rules._cnvkit_output.output.genemetrics_gene),
                str(rules._cnvkit_output.output.sex),
                str(rules._cnvkit_output.output.metrics),
                str(rules._cnvkit_output.output.seg)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            capture_space=CFG["runs"]["tumour_capture_space"]
        ),
        expand(
            expand(
            [
                str(rules._cnvkit_output_projection.output.projection)
            ],
            zip,  # Run expand() with zip(), not product()
            tumour_id=CFG["runs"]["tumour_sample_id"],
            seq_type=CFG["runs"]["tumour_seq_type"],
            allow_missing=True),
            tool = "cnvkit",
            projection=CFG["requested_projections"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
