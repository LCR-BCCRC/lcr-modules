#!/usr/bin/env snakemake


##### ATTRIBUTION #####
# Translated from the CCSRI_1500 exomes PureCN pipeline
# (/projects/dscott_prj/CCSRI_1500/exomes/snakefile_purecn.smk, J. Wong).
#
# Builds the PureCN-side panel of normals, ONCE per
# {seq_type}--{genome_build}/{capture_space}, from every tissue_status==normal in
# the samples table (so an external/unpaired set such as the normal exomes in
# gambl_metadata_normals.tsv can be injected by the wrapper). Emits:
#   - baits intervals (IntervalFile.R)
#   - mapping_bias.rds            (NormalDB.R --normal-panel)
#   - denovo normalDB.rds         (NormalDB.R --normal-panel --coverage-files)
#   - mutect2_pon.vcf.gz          (GenomicsDBImport -> CreateSomaticPanelOfNormals)
# The cnvkit coverage reference (normal_reference.cnn) is built by the cnvkit
# module from the same normals (cnvkit reads its PON members from the samples
# table too), so both tools draw from one normal set.


##### SETUP #####

import oncopipe as op

CFG = op.setup_module(
    name = "panel_of_normals",
    version = "1.0",
    subdirectories = ["inputs", "references", "coverage", "mutect2", "merged",
                      "normaldb", "pon", "outputs"],
)

GENOME_MAP = CFG["options"]["genome_builds_map"]


def _pon_chroms(genome_build):
    prefix = "chr" if "hg" in str(genome_build) else ""
    return [prefix + str(c) for c in list(range(1, 23)) + ["X", "Y"]]


# distinct PON groups (seq_type, genome_build, capture_space) from the normals
def _pon_groups():
    s = CFG["samples"]
    n = s[s["tissue_status"] == "normal"]
    return list(dict.fromkeys(zip(n["seq_type"], n["genome_build"], n["capture_space"])))


def _pon_get_normals(wildcards, template):
    s = CFG["samples"]
    n = s[(s["tissue_status"] == "normal") &
          (s["seq_type"] == wildcards.seq_type) &
          (s["genome_build"] == wildcards.genome_build) &
          (s["capture_space"] == wildcards.capture_space)]
    return expand(template, sample_id = list(dict.fromkeys(n["sample_id"])),
                  seq_type = wildcards.seq_type, genome_build = wildcards.genome_build,
                  capture_space = wildcards.capture_space)


localrules:
    _pon_input_bam,
    _pon_get_mappability,
    _pon_coverage_list,
    _pon_samples_map,
    _pon_output,
    _pon_all


##### 00 INPUTS / REFERENCES #####

rule _pon_input_bam:
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


def _pon_mappability_src(wildcards):
    return CFG["options"]["mappability_bw"].get(wildcards.genome_build, "")

rule _pon_get_mappability:
    output:
        bw = CFG["dirs"]["references"] + "mappability/{genome_build}.bw"
    params:
        src = _pon_mappability_src
    conda:
        CFG["conda_envs"]["wget"]
    shell:
        op.as_one_line("""
        if [[ "{params.src}" == http* ]]; then wget -qO {output.bw} "{params.src}";
        else ln -srf "{params.src}" {output.bw}; fi
        """)

rule _pon_setinterval:
    input:
        bed = lambda w: reference_files("genomes/" + w.genome_build + "/capture_space/" + w.capture_space + ".bed"),
        fasta = lambda w: reference_files("genomes/" + w.genome_build + "/genome_fasta/genome.fa"),
        bw = lambda w: (str(rules._pon_get_mappability.output.bw).replace("{genome_build}", w.genome_build)
                        if CFG["options"]["mappability_bw"].get(w.genome_build, "") else []),
    output:
        intervals = CFG["dirs"]["references"] + "{seq_type}--{genome_build}/{capture_space}/baits_intervals.txt"
    params:
        genome = lambda w: GENOME_MAP[w.genome_build],
        mapflag = lambda w, input: ("--mappability " + str(input.bw)) if input.bw else "",
        opts = CFG["options"]["setinterval"]["opts"],
    conda:
        CFG["conda_envs"]["purecn"]
    threads:
        CFG["threads"]["setinterval"]
    resources:
        **CFG["resources"]["setinterval"]
    log:
        CFG["logs"]["references"] + "{seq_type}--{genome_build}/{capture_space}/setinterval.log"
    shell:
        op.as_one_line("""
        PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
        Rscript --vanilla $PURECN/IntervalFile.R --in-file {input.bed} --fasta {input.fasta}
            --out-file {output.intervals} --genome {params.genome} {params.mapflag}
            --force {params.opts} > {log} 2>&1
        """)


##### 01 PER-NORMAL COVERAGE (PureCN Coverage.R) #####

rule _pon_coverage:
    input:
        bam = str(rules._pon_input_bam.output.bam),
        intervals = str(rules._pon_setinterval.output.intervals),
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


##### 02 MUTECT2 GERMLINE ON NORMALS #####

rule _pon_mutect2_chrom:
    input:
        bam = str(rules._pon_input_bam.output.bam),
        fasta = lambda w: reference_files("genomes/{genome_build}/genome_fasta/genome.fa".format(genome_build = w.genome_build)),
        gnomad = lambda w: reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz".format(genome_build = w.genome_build)),
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.vcf.gz.tbi"),
        stats = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.{chrom}.vcf.gz.stats"),
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
            --genotype-pon-sites true --interval-padding 50 --max-mnp-distance 0
            --germline-resource {input.gnomad} -R {input.fasta} -L {wildcards.chrom}
            -I {input.bam} -O {output.vcf} > {log} 2>&1
        """)

def _pon_chr_vcfs(wildcards):
    return expand(str(rules._pon_mutect2_chrom.output.vcf),
                  chrom = _pon_chroms(wildcards.genome_build), **wildcards)

rule _pon_mutect2_concat:
    input:
        vcf = _pon_chr_vcfs,
    output:
        vcf = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.vcf.gz",
        tbi = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}/{sample_id}.vcf.gz.tbi",
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    shell:
        "bcftools concat {input.vcf} -Oz -o {output.vcf} && tabix -p vcf {output.vcf}"

rule _pon_mutect2_pass:
    input:
        vcf = str(rules._pon_mutect2_concat.output.vcf),
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


##### 03 MERGED NORMAL PANEL + NORMAL DATABASES #####

rule _pon_merge_normal_vcfs:
    input:
        vcfs = lambda w: _pon_get_normals(w, str(rules._pon_mutect2_pass.output.vcf)),
        tbis = lambda w: _pon_get_normals(w, str(rules._pon_mutect2_pass.output.tbi)),
    output:
        vcf = CFG["dirs"]["merged"] + "{seq_type}--{genome_build}/{capture_space}/normal_panel.vcf.gz",
        tbi = CFG["dirs"]["merged"] + "{seq_type}--{genome_build}/{capture_space}/normal_panel.vcf.gz.tbi",
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        "bcftools merge {input.vcfs} -Oz -o {output.vcf} --force-samples && tabix -p vcf {output.vcf}"

rule _pon_coverage_list:
    input:
        cov = lambda w: _pon_get_normals(w, str(rules._pon_coverage.output.coverage)),
    output:
        lst = CFG["dirs"]["normaldb"] + "{seq_type}--{genome_build}/{capture_space}/cov_list.txt",
    shell:
        "ls {input.cov} > {output.lst}"

rule _pon_normaldb_mappingbias:
    input:
        normal_panel = str(rules._pon_merge_normal_vcfs.output.vcf),
        tbi = str(rules._pon_merge_normal_vcfs.output.tbi),
    output:
        mapping_bias = CFG["dirs"]["normaldb"] + "{seq_type}--{genome_build}/{capture_space}/mapping_bias.rds",
    params:
        outdir = CFG["dirs"]["normaldb"] + "{seq_type}--{genome_build}/{capture_space}",
        genome = lambda w: GENOME_MAP[w.genome_build],
        assay = "{capture_space}",
    conda:
        CFG["conda_envs"]["purecn"]
    threads:
        CFG["threads"]["normaldb"]
    resources:
        **CFG["resources"]["normaldb"]
    log:
        CFG["logs"]["normaldb"] + "{seq_type}--{genome_build}/{capture_space}/mapping_bias.log"
    shell:
        op.as_one_line("""
        PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
        Rscript --vanilla $PURECN/NormalDB.R --out-dir {params.outdir}
            --normal-panel {input.normal_panel} --assay {params.assay}
            --genome {params.genome} --force > {log} 2>&1 ;
        mv {params.outdir}/mapping_bias_*.rds {output.mapping_bias}
        """)

rule _pon_normaldb_denovo:
    input:
        normal_panel = str(rules._pon_merge_normal_vcfs.output.vcf),
        tbi = str(rules._pon_merge_normal_vcfs.output.tbi),
        cov_list = str(rules._pon_coverage_list.output.lst),
    output:
        mapping_bias = CFG["dirs"]["normaldb"] + "{seq_type}--{genome_build}/{capture_space}/denovo/mapping_bias.rds",
        normal_db = CFG["dirs"]["normaldb"] + "{seq_type}--{genome_build}/{capture_space}/denovo/normalDB.rds",
    params:
        outdir = CFG["dirs"]["normaldb"] + "{seq_type}--{genome_build}/{capture_space}/denovo",
        genome = lambda w: GENOME_MAP[w.genome_build],
        assay = "{capture_space}",
    conda:
        CFG["conda_envs"]["purecn"]
    threads:
        CFG["threads"]["normaldb"]
    resources:
        **CFG["resources"]["normaldb"]
    log:
        CFG["logs"]["normaldb"] + "{seq_type}--{genome_build}/{capture_space}/denovo.log"
    shell:
        op.as_one_line("""
        PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
        Rscript --vanilla $PURECN/NormalDB.R --out-dir {params.outdir}
            --normal-panel {input.normal_panel} --coverage-files {input.cov_list}
            --assay {params.assay} --genome {params.genome} --force > {log} 2>&1 ;
        mv {params.outdir}/mapping_bias_*.rds {output.mapping_bias} ;
        mv {params.outdir}/normalDB_*.rds {output.normal_db}
        """)


##### 04 MUTECT2 PON VCF #####

rule _pon_samples_map:
    input:
        vcfs = lambda w: _pon_get_normals(w, str(rules._pon_mutect2_concat.output.vcf)),
    output:
        map = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/samples_map.txt",
    shell:
        op.as_one_line("""
        rm -f {output.map};
        for v in {input.vcfs}; do n=$(basename $v .vcf.gz);
            printf "%s\\t%s\\n" "$n" "$v" >> {output.map}; done
        """)

rule _pon_genomicsdb:
    input:
        vcfs = lambda w: _pon_get_normals(w, str(rules._pon_mutect2_concat.output.vcf)),
        tbis = lambda w: _pon_get_normals(w, str(rules._pon_mutect2_concat.output.tbi)),
        map = str(rules._pon_samples_map.output.map),
        bed = lambda w: reference_files("genomes/" + w.genome_build + "/capture_space/" + w.capture_space + ".bed"),
    output:
        done = touch(CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/genomicsdb.done"),
    params:
        mem_mb = lambda w, resources: int(resources.mem_mb * 0.8),
        db = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/genomicsdb",
    conda:
        CFG["conda_envs"]["gatk"]
    resources:
        **CFG["resources"]["mutect2"]
    log:
        CFG["logs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/genomicsdb.log"
    shell:
        op.as_one_line("""
        gatk GenomicsDBImport --java-options "-Xmx{params.mem_mb}m" -L {input.bed}
            --sample-name-map {input.map} --genomicsdb-workspace-path {params.db}
            --lenient --merge-input-intervals TRUE
            --overwrite-existing-genomicsdb-workspace TRUE > {log} 2>&1
        """)

rule _pon_create_pon:
    input:
        done = str(rules._pon_genomicsdb.output.done),
        fasta = lambda w: reference_files("genomes/{genome_build}/genome_fasta/genome.fa".format(genome_build = w.genome_build)),
    output:
        pon = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/mutect2_pon.vcf.gz",
    params:
        mem_mb = lambda w, resources: int(resources.mem_mb * 0.8),
        db = "gendb://" + CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/genomicsdb",
    conda:
        CFG["conda_envs"]["gatk"]
    resources:
        **CFG["resources"]["mutect2"]
    log:
        CFG["logs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/create_pon.log"
    shell:
        op.as_one_line("""
        gatk CreateSomaticPanelOfNormals --java-options "-Xmx{params.mem_mb}m"
            --reference {input.fasta} --variant {params.db} -O {output.pon} > {log} 2>&1
        """)


##### 99 OUTPUTS #####

rule _pon_output:
    input:
        intervals = str(rules._pon_setinterval.output.intervals),
        mapping_bias = str(rules._pon_normaldb_mappingbias.output.mapping_bias),
        denovo_bias = str(rules._pon_normaldb_denovo.output.mapping_bias),
        normal_db = str(rules._pon_normaldb_denovo.output.normal_db),
        mutect2_pon = str(rules._pon_create_pon.output.pon),
    output:
        intervals = CFG["dirs"]["outputs"] + "intervals/{seq_type}--{genome_build}/{capture_space}/baits_intervals.txt",
        mapping_bias = CFG["dirs"]["outputs"] + "mapping_bias/{seq_type}--{genome_build}/{capture_space}/mapping_bias.rds",
        denovo_bias = CFG["dirs"]["outputs"] + "denovo/{seq_type}--{genome_build}/{capture_space}/mapping_bias.rds",
        normal_db = CFG["dirs"]["outputs"] + "denovo/{seq_type}--{genome_build}/{capture_space}/normalDB.rds",
        mutect2_pon = CFG["dirs"]["outputs"] + "mutect2_pon/{seq_type}--{genome_build}/{capture_space}/mutect2_pon.vcf.gz",
    run:
        op.relative_symlink(input.intervals, output.intervals, in_module = True)
        op.relative_symlink(input.mapping_bias, output.mapping_bias, in_module = True)
        op.relative_symlink(input.denovo_bias, output.denovo_bias, in_module = True)
        op.relative_symlink(input.normal_db, output.normal_db, in_module = True)
        op.relative_symlink(input.mutect2_pon, output.mutect2_pon, in_module = True)


rule _pon_all:
    input:
        [expand(
            [
                str(rules._pon_output.output.intervals),
                str(rules._pon_output.output.mapping_bias),
                str(rules._pon_output.output.normal_db),
                str(rules._pon_output.output.mutect2_pon),
            ],
            seq_type = st, genome_build = gb, capture_space = cs)
         for (st, gb, cs) in _pon_groups()]


##### CLEANUP #####
op.cleanup_module(CFG)
