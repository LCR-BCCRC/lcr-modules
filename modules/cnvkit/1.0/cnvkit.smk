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
    subdirectories = ["inputs", "cnvkit", "cn_calls", "SNPs", "BAF", "scatter", "seg", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _cnvkit_input_bam,
    _cnvkit_symlink_cns,
    _cnvkit_output,
    _cnvkit_all,


##### RULES #####
if CFG["options"]["hard_masked"] == True:
    CFG["runs"]["masked"] = "_masked"
else:
    CFG["runs"]["masked"] = ""


# cnvkit reference files
# gene annotation files
rule _get_refFlat:
    output:
        refFlat = CFG["dirs"]["inputs"] + "{genome_build}/refFlat.txt"
    params:
        url = "http://hgdownload.soe.ucsc.edu/goldenPath/",
        build = lambda w: "hg38" if "38" in str({w.genome_build}) else "hg19",
        gzip = CFG["dirs"]["inputs"] + "ref/refFlat.txt.gz",
        prefix = lambda w: "" if "grch" in str({w.genome_build}) else "chr",
        tmp = CFG["dirs"]["inputs"] + "ref/refFlat.txt.tmp"
    shell:
        """
            wget -O {params.gzip} {params.url}{params.build}/database/refFlat.txt.gz &&
            gzip -d {params.gzip} &&
            sed 's/chr/{params.prefix}/g' {output.refFlat} > {params.tmp} &&
            mv {params.tmp} {output.refFlat}
        """
    

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _cnvkit_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam.crai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
        op.relative_symlink(input.bai, output.crai)


def get_normals(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    capture_space = CFG["runs"][CFG["runs"]["normal_capture_space"].isin([wildcards.capture_space])]
    normals = expand(CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}.bam", 
                    zip,
                    capture_space = capture_space["normal_capture_space"],
                    seq_type = capture_space["normal_seq_type"],
                    genome_build = capture_space["normal_genome_build"],
                    normal_id = capture_space["normal_sample_id"])
    normals = list(dict.fromkeys(normals))
    return normals

def _cnvkit_get_capspace(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    if str(wildcards.capture_space) in CFG["options"]["target_bed"]:
        bed = op.switch_on_wildcard("capture_space", CFG["options"]["target_bed"])
    else:
        bed = op.get_capture_space(CFG, wildcards.normal_id, wildcards.genome_build, wildcards.seq_type, "bed")
    return bed
    

rule _run_cnvkit_batch:
    input:
        samples = expand(CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.bam", 
                        zip,
                        tumour_id = CFG["runs"]["tumour_sample_id"],
                        capture_space = CFG["runs"]["tumour_capture_space"],
                        genome_build = CFG["runs"]["tumour_genome_build"],
                        seq_type = CFG["runs"]["tumour_seq_type"]),
        controls = get_normals,
        refFlat = str(rules._get_refFlat.output.refFlat),
        fasta =  reference_files("genomes/{genome_build}{masked}/genome_fasta/genome.fa"),
    output: 
        panel_of_norms = CFG["dirs"]["cnvkit"] + "{seq_type}--{genome_build}{masked}/{capture_space}/normal_reference.cnn",
        cnr = expand(
                CFG["dirs"]["cnvkit"] + "{{seq_type}}--{{genome_build}}{{masked}}/{{capture_space}}/{tumour_id}.cnr", zip,
                tumour_id = CFG["runs"]["tumour_sample_id"]),
        cns = expand(
                CFG["dirs"]["cnvkit"] + "{{seq_type}}--{{genome_build}}{{masked}}/{{capture_space}}/{tumour_id}.cns", zip,
                tumour_id = CFG["runs"]["tumour_sample_id"]),
        call_cns = expand(
                CFG["dirs"]["cnvkit"] + "{{seq_type}}--{{genome_build}}{{masked}}/{{capture_space}}/{tumour_id}.call.cns", zip,
                tumour_id = CFG["runs"]["tumour_sample_id"]),
    params:
        outdir = CFG["dirs"]["cnvkit"] + "{seq_type}--{genome_build}{masked}/{capture_space}/",
        # targets = op.switch_on_wildcard("capture_space", CFG["options"]["target_bed"]),
        targets = _cnvkit_get_capspace,
        opts = CFG["options"]["batch"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["batch"]
    resources:
        **CFG["resources"]["batch"]
    wildcard_constraints:
        masked = ".{0}|_masked",
        genome_build = ".+(?<!masked)"
    log: 
        stdout = CFG["logs"]["cnvkit"] + "{seq_type}--{genome_build}{masked}/{capture_space}_cnr.log"
    shell:
        op.as_one_line("""cnvkit.py batch {input.samples} --normal {input.controls} --fasta {input.fasta} --annotate {input.refFlat} --processes {threads} --output-reference {output.panel_of_norms} --targets {params.targets} --output-dir {params.outdir} {params.opts} > {log} 2&>1 """)


rule _cnvkit_symlink_cns:
    input:
        call_cns = CFG["dirs"]["cnvkit"] + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.call.cns"
    output:
        call_cns = CFG["dirs"]["cn_calls"] + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.call.cns"
    wildcard_constraints:
        masked = ".{0}|_masked",
        genome_build = ".+(?<!masked)"
    run:
        op.relative_symlink(input.call_cns, output.call_cns, in_module = True)


# need SNPs not SNVs (i.e. get germline calls using a dbSNP vcf)
rule _cnvkit_dbsnp_to_bed:
    input:
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        bed = CFG["dirs"]["SNPs"] + "{genome_build}/dbsnp.common_all-151.bed"
    resources: **CFG["resources"]["SNPs"]
    shell:
        op.as_one_line(""" gunzip -c {input.vcf} | awk {{'printf ("%s\\t%s\\t%s\\n", $1,$2-1,$2)'}} | zgrep -v -h "^#" > {output.bed} """)

#### set-up mpileups for BAF calling ####
def _cnvkit_get_chr_mpileups(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    mpileups = expand(
        CFG["dirs"]["SNPs"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}.{chrom}.vcf.gz", 
        chrom = chrs
    )
    return(mpileups)
    
def _cnvkit_get_chr_mpileups_tbi(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    chrs = reference_files("genomes/" + wildcards.genome_build + "/genome_fasta/main_chromosomes_withY.txt")
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    mpileups = expand(
        CFG["dirs"]["SNPs"] + "{{seq_type}}--{{genome_build}}/{{capture_space}}/{{sample_id}}.{chrom}.vcf.gz.tbi", 
        chrom = chrs
    )
    return(mpileups)

# vcf needs DP, GT, AD - bcftools -mv calls multiallelic variants (will annotate GT)
# without it, GT will not be annotated
rule _cnvkit_mpileup_per_chrom:
    input:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam",
        fastaFile = reference_files("genomes/{genome_build}{masked}/genome_fasta/genome.fa"),
        bed = str(rules._cnvkit_dbsnp_to_bed.output.bed)
    output: # creates a temporary file for mpileup
        vcf = temp(CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}{masked}/{capture_space}/{sample_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}{masked}/{capture_space}/{sample_id}.{chrom}.vcf.gz.tbi")
    params:
        quality = CFG["options"]["SNPs"]["quality"],
        opts = CFG["options"]["SNPs"]["opts"]
    conda:
        CFG["conda_envs"]["bcftools"]
    resources: 
        **CFG["resources"]["SNPs"]
    group: "cnvkit"
    wildcard_constraints:
        masked = ".{0}|_masked",
        genome_build = ".+(?<!masked)"
    log:
        stderr = CFG["logs"]["SNPs"] + "{capture_space}/{seq_type}--{genome_build}{masked}/{sample_id}/{chrom}.vcf.stderr.log",
    shell:
        """
            bcftools mpileup -T {input.bed} -r {wildcards.chrom} -f {input.fastaFile} -Q {params.quality} {params.opts} -Ou {input.bam} | bcftools call -mv -Oz -o {output.vcf} 2> {log.stderr} && 
            tabix -p vcf {output.vcf}
        """
        

rule _cnvkit_concatenate_vcf:
    input: 
        vcf = _cnvkit_get_chr_mpileups,
        tbi = _cnvkit_get_chr_mpileups_tbi,
    output: 
        vcf = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}{masked}/{capture_space}/{sample_id}.vcf.gz",
        tbi = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}{masked}/{capture_space}/{sample_id}.vcf.gz.tbi"
    resources: 
        **CFG["resources"]["SNPs"]
    group: "cnvkit"
    conda:
        CFG["conda_envs"]["bcftools"]
    wildcard_constraints:
        masked = ".{0}|_masked",
        genome_build = ".+(?<!masked)"
    shell: 
        """
            bcftools concat {input.vcf} -Oz -o {output.vcf} && 
            tabix -p vcf {output.vcf}
        """


# ----------------------------------------------------------------------------------------------- #
# Integrating cnvkit with BAF
# ----------------------------------------------------------------------------------------------- #
rule _run_cnvkit_call_vcf:
    input:        
        cns = CFG["dirs"]["cn_calls"] + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.call.cns",
        vcf = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.vcf.gz",
        tbi = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.vcf.gz.tbi"
    output: 
        cns =  CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.call.cns"
    params:
        rescale = CFG["options"]["BAF"]["rescale"],
        min_depth = CFG["options"]["BAF"]["min_depth"],
        opts = CFG["options"]["BAF"]["opts"]
    log: 
        CFG["logs"]["BAF"] + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}_call.log"
    group: "cnvkit"
    conda: 
        CFG["conda_envs"]["cnvkit"]
    wildcard_constraints:
        masked = ".{0}|_masked",
        genome_build = ".+(?<!masked)"
    shell:
        """
            cnvkit.py call {input.cns} --output {output.cns} -v {input.vcf} --min-variant-depth {params.min_depth} -m {params.rescale} {params.opts} &> {log}
        """


rule _run_cnvkit_scatter:
    input:
        cnr = CFG["dirs"]["cnvkit"] + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.cnr",
        cns =  CFG["dirs"]["BAF"] + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.call.cns", 
        vcf = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.vcf.gz",
        tbi = CFG["dirs"]["SNPs"]  + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.vcf.gz.tbi"
    output: 
        png = CFG["dirs"]["scatter"] + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}_scatter.png"
    params:
        min_depth = CFG["options"]["scatter"]["min_depth"],
        ymax = CFG["options"]["scatter"]["ymax"],
        ymin = CFG["options"]["scatter"]["ymin"]
    conda: 
        CFG["conda_envs"]["cnvkit"]
    group: "cnvkit"
    wildcard_constraints:
        masked = ".{0}|_masked",
        genome_build = ".+(?<!masked)"
    shell:
        "cnvkit.py scatter {input.cnr} -s {input.cns} --output {output.png} -v {input.vcf} --min-variant-depth {params.min_depth} --y-max {params.ymax} --y-min {params.ymin}"


rule _cnvkit_to_seg:
    input:
        cns = str(rules._run_cnvkit_call_vcf.output.cns),
    output:
        seg = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.seg"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["seg"]
    resources:
        **CFG["resources"]["seg"]
    wildcard_constraints:
        masked = ".{0}|_masked",
        genome_build = ".+(?<!masked)"
    shell:
        """cnvkit.py export seg {input.cns} -o {output.seg} """


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _cnvkit_output:
    input:
        call_cns = str(rules._run_cnvkit_call_vcf.output.cns),
        scatter = str(rules._run_cnvkit_scatter.output.png),
        seg = str(rules._cnvkit_to_seg.output.seg)
    output:
        call_cns = CFG["dirs"]["outputs"] + "capture_space/BAF_cns/{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.call.cns",
        scatter = CFG["dirs"]["outputs"] + "capture_space/plots/{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}_scatter.png",
        seg = CFG["dirs"]["outputs"] + "capture_space/seg/{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.seg"
    wildcard_constraints:
        masked = ".{0}|_masked",
        genome_build = ".+(?<!masked)"
    run:
        op.relative_symlink(input.call_cns, output.call_cns, in_module = True)
        op.relative_symlink(input.scatter, output.scatter, in_module = True)
        op.relative_symlink(input.seg, output.seg, in_module = True)


# this is to collapse the wildcard {capture_space} to streamline downstream analyses
rule _cnvkit_output_no_capture_space:
    input:
        call_cns = str(rules._run_cnvkit_call_vcf.output.cns),
        scatter = str(rules._run_cnvkit_scatter.output.png),
        seg = str(rules._cnvkit_to_seg.output.seg)
    output:
        call_cns = touch(CFG["dirs"]["outputs"] + "no_capture_space/done/BAF_cns/{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.call.cns"),
        scatter = touch(CFG["dirs"]["outputs"] + "no_capture_space/done/plots/{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}_scatter.png"),
        seg = touch(CFG["dirs"]["outputs"] + "no_capture_space/done/seg/{seq_type}--{genome_build}{masked}/{capture_space}/{tumour_id}.seg")
    params:
        call_cns = CFG["dirs"]["outputs"] + "no_capture_space/BAF_cns/{seq_type}--{genome_build}{masked}/{tumour_id}.call.cns",
        scatter = CFG["dirs"]["outputs"] + "no_capture_space/plots/{seq_type}--{genome_build}{masked}/{tumour_id}_scatter.png",
        seg = CFG["dirs"]["outputs"] + "no_capture_space/seg/{seq_type}--{genome_build}{masked}/{tumour_id}.seg"
    wildcard_constraints:
        masked = ".{0}|_masked",
        genome_build = ".+(?<!masked)"
    run:
        op.relative_symlink(input.call_cns, params.call_cns, in_module = True)
        op.relative_symlink(input.scatter, params.scatter, in_module = True)
        op.relative_symlink(input.seg, params.seg, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _cnvkit_all:
    input:
        expand(
            [
                str(rules._cnvkit_output_no_capture_space.output.call_cns),
                str(rules._cnvkit_output_no_capture_space.output.scatter),
                str(rules._cnvkit_output_no_capture_space.output.seg),
                str(rules._cnvkit_output.output.call_cns),
                str(rules._cnvkit_output.output.scatter),
                str(rules._cnvkit_output.output.seg),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"],
            capture_space=CFG["runs"]["tumour_capture_space"],
            masked=CFG["runs"]["masked"],
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
