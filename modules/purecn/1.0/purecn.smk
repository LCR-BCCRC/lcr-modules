#!/usr/bin/env snakemake


##### ATTRIBUTION #####
# Translated from the CCSRI_1500 exomes PureCN pipeline
# (/projects/dscott_prj/CCSRI_1500/exomes/snakefile_purecn.smk, J. Wong).
#
# Tumour-only purity / ploidy / integer copy number / LOH with PureCN. Consumes:
#   - the cnvkit module (.cnr + BAF call.cns), and
#   - the panel_of_normals module (intervals, mapping_bias, denovo normalDB,
#     Mutect2 PON), all keyed by {seq_type}--{genome_build}/{capture_space}.
# Runs PureCN in two modes per tumour (cnvkit mode + denovo mode).


##### SETUP #####

import oncopipe as op

CFG = op.setup_module(
    name = "purecn",
    version = "1.0",
    subdirectories = ["inputs", "coverage", "mutect2", "purecn_cnvkit",
                      "purecn_denovo", "outputs"],
)

GENOME_MAP = CFG["options"]["genome_builds_map"]


def _purecn_chroms(genome_build):
    prefix = "chr" if "hg" in str(genome_build) else ""
    return [prefix + str(c) for c in list(range(1, 23)) + ["X", "Y"]]


localrules:
    _purecn_input_bam,
    _purecn_input_cnvkit,
    _purecn_setup_blacklist,
    _purecn_export_seg,
    _purecn_output,
    _purecn_all


##### 00 INPUTS #####

rule _purecn_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"],
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam.crai",
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)

rule _purecn_input_cnvkit:
    input:
        cnr = CFG["inputs"]["cnvkit_cnr"],
        cns = CFG["inputs"]["cnvkit_cns"],
    output:
        cnr = CFG["dirs"]["inputs"] + "cnvkit/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr",
        cns = CFG["dirs"]["inputs"] + "cnvkit/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.call.cns",
    run:
        op.absolute_symlink(input.cnr, output.cnr)
        op.absolute_symlink(input.cns, output.cns)

# SNP blacklist (repeatmasker, reduced to 3 columns).
rule _purecn_setup_blacklist:
    input:
        blacklist = lambda w: reference_files("genomes/" + w.genome_build + "/repeatmasker/repeatmasker.{genome_build}.bed".format(genome_build = w.genome_build)),
    output:
        blacklist = CFG["dirs"]["inputs"] + "references/{genome_build}/repeatmasker.3col.bed"
    shell:
        """awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3}}' {input.blacklist} > {output.blacklist}"""


##### 01 TUMOUR COVERAGE (PureCN Coverage.R, for denovo mode) #####

rule _purecn_coverage:
    input:
        bam = str(rules._purecn_input_bam.output.bam),
        intervals = CFG["inputs"]["pon_intervals"],
    output:
        coverage = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}_coverage_loess.txt.gz"
    params:
        outdir = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}",
        opts = CFG["options"]["coverage"]["opts"],
    conda:
        CFG["conda_envs"]["purecn"]
    threads:
        CFG["threads"]["coverage"]
    resources:
        **CFG["resources"]["coverage"]
    log:
        CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.log"
    shell:
        op.as_one_line("""
        PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
        Rscript --vanilla $PURECN/Coverage.R --out-dir {params.outdir}
            --bam {input.bam} --intervals {input.intervals} --force {params.opts} > {log} 2>&1 ;
        mv {params.outdir}/{wildcards.sample_id}*_loess.txt.gz {output.coverage}
        """)


##### 02 MUTECT2 ON TUMOUR (BAF VCF, scored against the PON) #####

rule _purecn_mutect2_chrom:
    input:
        bam = str(rules._purecn_input_bam.output.bam),
        fasta = lambda w: reference_files("genomes/{genome_build}/genome_fasta/genome.fa".format(genome_build = w.genome_build)),
        gnomad = lambda w: reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz".format(genome_build = w.genome_build)),
        pon = CFG["inputs"]["pon_mutect2_pon"],
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.vcf.gz.tbi"),
        stats = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.vcf.gz.stats"),
        f1r2 = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.f1r2.tar.gz"),
    params:
        mem_mb = lambda w, resources: int(resources.mem_mb * 0.8),
    conda:
        CFG["conda_envs"]["gatk"]
    resources:
        **CFG["resources"]["mutect2"]
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{chrom}.log"
    shell:
        op.as_one_line("""
        gatk Mutect2 --java-options "-Xmx{params.mem_mb}m" --genotype-germline-sites true
            --genotype-pon-sites true --interval-padding 50 --germline-resource {input.gnomad}
            -R {input.fasta} -L {wildcards.chrom} -pon {input.pon} -I {input.bam}
            -O {output.vcf} --f1r2-tar-gz {output.f1r2} > {log} 2>&1
        """)

def _purecn_tumour_chr(wildcards, which):
    out = getattr(rules._purecn_mutect2_chrom.output, which)
    return expand(str(out), chrom = _purecn_chroms(wildcards.genome_build), **wildcards)

rule _purecn_mutect2_concat:
    input:
        vcf = lambda w: _purecn_tumour_chr(w, "vcf"),
    output:
        vcf = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.vcf.gz",
        tbi = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.vcf.gz.tbi",
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    shell:
        "bcftools concat {input.vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}"

rule _purecn_mutect2_merge_stats:
    input:
        stats = lambda w: _purecn_tumour_chr(w, "stats"),
    output:
        stats = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.vcf.gz.stats",
    conda:
        CFG["conda_envs"]["gatk"]
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/merge_stats.log"
    shell:
        op.as_one_line("""
        gatk MergeMutectStats $(for s in {input.stats}; do echo -n "-stats $s "; done)
            -O {output.stats} > {log} 2>&1
        """)

rule _purecn_mutect2_orient:
    input:
        f1r2 = lambda w: _purecn_tumour_chr(w, "f1r2"),
    output:
        model = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/read-orientation-model.tar.gz",
    params:
        mem_mb = lambda w, resources: int(resources.mem_mb * 0.8),
    conda:
        CFG["conda_envs"]["gatk"]
    resources:
        **CFG["resources"]["mutect2"]
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/orient.log"
    shell:
        op.as_one_line("""
        gatk LearnReadOrientationModel --java-options "-Xmx{params.mem_mb}m"
            $(for f in {input.f1r2}; do echo -n "-I $f "; done) -O {output.model} > {log} 2>&1
        """)

rule _purecn_mutect2_pileup:
    input:
        bam = str(rules._purecn_input_bam.output.bam),
        snps = lambda w: reference_files("genomes/{genome_build}/gatk/mutect2_small_exac.{genome_build}.vcf.gz".format(genome_build = w.genome_build)),
        fasta = lambda w: reference_files("genomes/{genome_build}/genome_fasta/genome.fa".format(genome_build = w.genome_build)),
    output:
        pileup = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/pileupSummary.table",
    params:
        mem_mb = lambda w, resources: int(resources.mem_mb * 0.8),
    conda:
        CFG["conda_envs"]["gatk"]
    resources:
        **CFG["resources"]["mutect2"]
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/pileup.log"
    shell:
        op.as_one_line("""
        gatk GetPileupSummaries --java-options "-Xmx{params.mem_mb}m" -I {input.bam}
            -R {input.fasta} -V {input.snps} -L {input.snps} -O {output.pileup} > {log} 2>&1
        """)

rule _purecn_mutect2_contam:
    input:
        pileup = str(rules._purecn_mutect2_pileup.output.pileup),
    output:
        segments = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/segments.table",
        contamination = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/contamination.table",
    params:
        mem_mb = lambda w, resources: int(resources.mem_mb * 0.8),
    conda:
        CFG["conda_envs"]["gatk"]
    resources:
        **CFG["resources"]["mutect2"]
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/contam.log"
    shell:
        op.as_one_line("""
        gatk CalculateContamination --java-options "-Xmx{params.mem_mb}m" -I {input.pileup}
            -tumor-segmentation {output.segments} -O {output.contamination} > {log} 2>&1
        """)

rule _purecn_mutect2_filter:
    input:
        vcf = str(rules._purecn_mutect2_concat.output.vcf),
        stats = str(rules._purecn_mutect2_merge_stats.output.stats),
        segments = str(rules._purecn_mutect2_contam.output.segments),
        contamination = str(rules._purecn_mutect2_contam.output.contamination),
        model = str(rules._purecn_mutect2_orient.output.model),
        fasta = lambda w: reference_files("genomes/{genome_build}/genome_fasta/genome.fa".format(genome_build = w.genome_build)),
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/unfilt.vcf.gz"),
    params:
        mem_mb = lambda w, resources: int(resources.mem_mb * 0.8),
    conda:
        CFG["conda_envs"]["gatk"]
    resources:
        **CFG["resources"]["mutect2"]
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/filter.log"
    shell:
        op.as_one_line("""
        gatk FilterMutectCalls --java-options "-Xmx{params.mem_mb}m" -V {input.vcf}
            -R {input.fasta} --tumor-segmentation {input.segments}
            --contamination-table {input.contamination} --ob-priors {input.model}
            -O {output.vcf} > {log} 2>&1
        """)

rule _purecn_mutect2_pass:
    input:
        vcf = str(rules._purecn_mutect2_filter.output.vcf),
    output:
        vcf = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}_passed.vcf.gz",
        tbi = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}_passed.vcf.gz.tbi",
    params:
        filter_for = CFG["options"]["mutect2"]["filter_for"],
        filter_out = CFG["options"]["mutect2"]["filter_out"],
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools view {params.filter_for} -e "{params.filter_out}" -Oz -o {output.vcf} {input.vcf} &&
        tabix -p vcf {output.vcf}
        """)


##### 03 PURECN #####

rule _purecn_export_seg:
    input:
        cns = str(rules._purecn_input_cnvkit.output.cns),
    output:
        seg = CFG["dirs"]["inputs"] + "cnvkit/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg",
    conda:
        CFG["conda_envs"]["cnvkit"]
    shell:
        op.as_one_line("""
        cnvkit.py export seg {input.cns} -o {output.seg}.tmp --enumerate-chroms &&
        sed 's/.call//g' {output.seg}.tmp > {output.seg} && rm {output.seg}.tmp
        """)

rule _purecn_filter_cnr:
    input:
        cnr = str(rules._purecn_input_cnvkit.output.cnr),
    output:
        cnr = temp(CFG["dirs"]["inputs"] + "cnvkit/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.filt.cnr"),
    params:
        min_log2 = CFG["options"]["min_cnr_log2"],
    shell:
        """awk 'NR==1 || ($5 > {params.min_log2})' {input.cnr} > {output.cnr}"""

# cnvkit mode: CNVkit .cnr + .seg + tumour VCF + PON mapping bias.
rule _purecn_run_cnvkit:
    input:
        cnr = str(rules._purecn_filter_cnr.output.cnr),
        seg = str(rules._purecn_export_seg.output.seg),
        vcf = str(rules._purecn_mutect2_pass.output.vcf).replace("{sample_id}", "{tumour_id}"),
        stats = str(rules._purecn_mutect2_merge_stats.output.stats).replace("{sample_id}", "{tumour_id}"),
        mapping_bias = CFG["inputs"]["pon_mapping_bias"],
        blacklist = str(rules._purecn_setup_blacklist.output.blacklist),
    output:
        ploidy = CFG["dirs"]["purecn_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.csv",
        seg = CFG["dirs"]["purecn_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg",
        loh = CFG["dirs"]["purecn_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_loh.csv",
    params:
        outdir = CFG["dirs"]["purecn_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}",
        genome = lambda w: GENOME_MAP[w.genome_build],
        seg_method = CFG["options"]["purecn"]["fun_segmentation"],
        max_cn = CFG["options"]["purecn"]["max_copy_number"],
        min_offtarget = CFG["options"]["purecn"]["min_fraction_offtarget"],
        alpha = CFG["options"]["purecn"]["alpha"],
        model = CFG["options"]["purecn"]["model"],
        opts = CFG["options"]["purecn"]["opts"],
    conda:
        CFG["conda_envs"]["purecn"]
    threads:
        CFG["threads"]["purecn"]
    resources:
        **CFG["resources"]["purecn"]
    log:
        CFG["logs"]["purecn_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.log"
    shell:
        op.as_one_line("""
        PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
        Rscript --vanilla $PURECN/PureCN.R --out {params.outdir} --sampleid {wildcards.tumour_id}
            --tumor {input.cnr} --seg-file {input.seg} --stats-file {input.stats}
            --vcf {input.vcf} --snp-blacklist {input.blacklist}
            --mappingbiasfile {input.mapping_bias} --genome {params.genome}
            --fun-segmentation {params.seg_method} --max-copy-number {params.max_cn}
            --min-fraction-offtarget {params.min_offtarget} --alpha {params.alpha}
            --model {params.model} --cores {threads} {params.opts} > {log} 2>&1
        """)

# denovo mode: PureCN Coverage.R output + PON denovo NormalDB + intervals.
rule _purecn_run_denovo:
    input:
        cnr = str(rules._purecn_coverage.output.coverage).replace("{sample_id}", "{tumour_id}"),
        vcf = str(rules._purecn_mutect2_pass.output.vcf).replace("{sample_id}", "{tumour_id}"),
        stats = str(rules._purecn_mutect2_merge_stats.output.stats).replace("{sample_id}", "{tumour_id}"),
        mapping_bias = CFG["inputs"]["pon_denovo_bias"],
        normal_db = CFG["inputs"]["pon_normal_db"],
        intervals = CFG["inputs"]["pon_intervals"],
        blacklist = str(rules._purecn_setup_blacklist.output.blacklist),
    output:
        ploidy = CFG["dirs"]["purecn_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.csv",
        seg = CFG["dirs"]["purecn_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg",
        loh = CFG["dirs"]["purecn_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_loh.csv",
    params:
        outdir = CFG["dirs"]["purecn_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}",
        genome = lambda w: GENOME_MAP[w.genome_build],
        seg_method = CFG["options"]["purecn"]["fun_segmentation"],
        max_cn = CFG["options"]["purecn"]["max_copy_number"],
        min_offtarget = CFG["options"]["purecn"]["min_fraction_offtarget"],
        alpha = CFG["options"]["purecn"]["alpha"],
        model = CFG["options"]["purecn"]["model"],
        opts = CFG["options"]["purecn"]["opts"],
    conda:
        CFG["conda_envs"]["purecn"]
    threads:
        CFG["threads"]["purecn"]
    resources:
        **CFG["resources"]["purecn"]
    log:
        CFG["logs"]["purecn_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.log"
    shell:
        op.as_one_line("""
        PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
        Rscript --vanilla $PURECN/PureCN.R --out {params.outdir} --sampleid {wildcards.tumour_id}
            --tumor {input.cnr} --stats-file {input.stats} --vcf {input.vcf}
            --snp-blacklist {input.blacklist} --mappingbiasfile {input.mapping_bias}
            --normaldb {input.normal_db} --intervals {input.intervals} --genome {params.genome}
            --fun-segmentation {params.seg_method} --max-copy-number {params.max_cn}
            --min-fraction-offtarget {params.min_offtarget} --alpha {params.alpha}
            --model {params.model} --cores {threads} {params.opts} > {log} 2>&1
        """)


##### 99 OUTPUTS #####

rule _purecn_output:
    input:
        seg = str(rules._purecn_run_cnvkit.output.seg),
        ploidy = str(rules._purecn_run_cnvkit.output.ploidy),
        loh = str(rules._purecn_run_cnvkit.output.loh),
        seg_d = str(rules._purecn_run_denovo.output.seg),
        ploidy_d = str(rules._purecn_run_denovo.output.ploidy),
        loh_d = str(rules._purecn_run_denovo.output.loh),
    output:
        seg = CFG["dirs"]["outputs"] + "purecn_cnvkit/seg/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_dnacopy.seg",
        ploidy = CFG["dirs"]["outputs"] + "purecn_cnvkit/ploidy/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.csv",
        loh = CFG["dirs"]["outputs"] + "purecn_cnvkit/loh/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_loh.csv",
        seg_d = CFG["dirs"]["outputs"] + "purecn_denovo/seg/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_dnacopy.seg",
        ploidy_d = CFG["dirs"]["outputs"] + "purecn_denovo/ploidy/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.csv",
        loh_d = CFG["dirs"]["outputs"] + "purecn_denovo/loh/{seq_type}--{genome_build}/{capture_space}/{tumour_id}_loh.csv",
    run:
        op.relative_symlink(input.seg, output.seg, in_module = True)
        op.relative_symlink(input.ploidy, output.ploidy, in_module = True)
        op.relative_symlink(input.loh, output.loh, in_module = True)
        op.relative_symlink(input.seg_d, output.seg_d, in_module = True)
        op.relative_symlink(input.ploidy_d, output.ploidy_d, in_module = True)
        op.relative_symlink(input.loh_d, output.loh_d, in_module = True)


rule _purecn_all:
    input:
        expand(
            [
                str(rules._purecn_output.output.seg),
                str(rules._purecn_output.output.ploidy),
                str(rules._purecn_output.output.seg_d),
                str(rules._purecn_output.output.ploidy_d),
            ],
            zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            capture_space = CFG["runs"]["tumour_capture_space"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
        )


##### CLEANUP #####
op.cleanup_module(CFG)
