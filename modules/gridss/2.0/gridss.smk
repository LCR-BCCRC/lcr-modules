#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

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
# `CFG` is a shortcut to `config["lcr-modules"]["gridss"]`
CFG = op.setup_module(
    name = "gridss",
    version = "2.0",
    subdirectories = ["inputs", "preprocess", "gridss", "repeatmasker", "gripss", "outputs"],
)

VERSION_MAP_GRIDSS = CFG["options"]["version_map"]

possible_genome_builds = ", ".join(list(VERSION_MAP_GRIDSS.keys()))
for genome_build in CFG["runs"]["tumour_genome_build"]:
    assert genome_build in possible_genome_builds, (
        f"Samples table includes genome builds not yet compatible with this module. "
        f"This module is currently only compatible with {possible_genome_builds}. "
    )

sample_ids = list(CFG['samples']['sample_id'])
unmatched_normal_ids = list(config["lcr-modules"]["_shared"]["unmatched_normal_ids"].values())
all_other_ids = list(set(sample_ids) - set(unmatched_normal_ids))

# Define rules to be run locally when using a compute cluster
localrules:
    _gridss_input_bam,
    _gridss_input_references,
    _gridss_setup_references,
    _gridss_get_pon,
    _gridss_symlink_preprocessed_normal,
    _gridss_filter_gripss,
    _gridss_gripss_to_bedpe,
    _gridss_output_somatic_vcf,
    _gridss_all



##### RULES #####

# Symlink genome fasta with bwa and .fai indices to the same directory
rule _gridss_input_references:
    input:
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        genome_bwa_prefix = reference_files("genomes/{genome_build}/bwa_index/bwa-0.7.17/genome.fa"),
    output:
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.fa",
    shell:
        op.as_one_line("""
        ln -sf {input.genome_fa} {output.genome_fa} &&
        ln -sf {input.genome_fa}.fai {output.genome_fa}.fai &&
        ln -sf {input.genome_bwa_prefix}.* `dirname {output.genome_fa}`
        """)

# Download the panel of normals
rule _gridss_get_pon:
    output:
        pon_breakpoint = CFG["dirs"]["inputs"] + "references/{genome_build}/pon/gridss_pon_breakpoint.bedpe",
        pon_breakend = CFG["dirs"]["inputs"] + "references/{genome_build}/pon/gridss_pon_single_breakend.bed",
        known_pairs = CFG["dirs"]["inputs"] + "references/{genome_build}/pon/KnownFusionPairs.bedpe"
    params:
        alt_build = lambda w: VERSION_MAP_GRIDSS[w.genome_build],
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/gridss/pon"
    shell:
        op.as_one_line("""
        wget -O {output.pon_breakpoint} {params.url}/gridss_pon_breakpoint.{params.alt_build}.bedpe;
        wget -O {output.pon_breakend} {params.url}/gridss_pon_single_breakend.{params.alt_build}.bed;
        wget -O {output.known_pairs} {params.url}/KnownFusionPairs.{params.alt_build}.bedpe
        """)


# Generage genome.fa.img file
rule _gridss_setup_references:
    input:
        fasta = str(rules._gridss_input_references.output.genome_fa),
    output:
        genome_img = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.fa.img"
    params:
        steps = "setupreference"
    conda:
        CFG["conda_envs"]["gridss"]
    resources:
        mem_mb = 4000
    threads: 8
    shell:
        op.as_one_line("""
        gridss
        --reference {input.fasta}
        --threads {threads}
        --jvmheap 3G
        --steps {params.steps}
        --workingdir `dirname {output.genome_img}`
        """)


# Symlink the input files into the module results directory (under '00-inputs/')
rule _gridss_input_bam:
    input:
        sample_bam = ancient(CFG["inputs"]["sample_bam"]),
        sample_bai = ancient(CFG["inputs"]["sample_bai"])
    output:
        sample_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        sample_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.sample_bam, output.sample_bam)
        op.absolute_symlink(input.sample_bai, output.sample_bai)

# Preprocess unmatched normal bams
rule _gridss_preprocess_unmatched_normal:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = str(rules._gridss_input_references.output.genome_fa),
        fasta_img = str(rules._gridss_setup_references.output.genome_img)
    output:
        workdir = directory(CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}.bam.gridss.working")
    log: CFG["logs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/preprocess.log"
    params:
        opts = CFG["options"]["gridss"],
        steps = "preprocess",
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gridss"]
    threads:
        CFG["threads"]["gridss"]
    resources:
        **CFG["resources"]["gridss"]
    priority: 1
    wildcard_constraints:
        sample_id="|".join(unmatched_normal_ids)
    shell:
        op.as_one_line("""
        gridss
        --reference {input.fasta}
        --workingdir $(dirname {output.workdir})
        --threads {threads}
        --jvmheap {params.mem_mb}m
        --steps {params.steps}
        {params.opts}
        {input.bam}
        2>&1 | tee -a {log}
        """)

# Symlink preprocessed sv.bam directories

rule _gridss_symlink_preprocessed_normal:
    input:
        workdir = str(rules._gridss_preprocess_unmatched_normal.output.workdir)
    output:
        workdir = temp(directory(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/{sample_id}.bam.gridss.working"))
    priority: 0
    wildcard_constraints:
        sample_id = "|".join(unmatched_normal_ids)
    run:
        op.absolute_symlink(input.workdir, output.workdir)

# Preprocess all other bams as part of the group job
rule _gridss_preprocess:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = str(rules._gridss_input_references.output.genome_fa),
        fasta_img = str(rules._gridss_setup_references.output.genome_img)
    output:
        workdir = temp(directory(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/{sample_id}.bam.gridss.working"))
    log: CFG["logs"]["preprocess"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/{sample_id}/preprocess.log"
    params:
        opts = CFG["options"]["gridss"],
        steps = "preprocess",
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gridss"]
    threads:
        CFG["threads"]["preprocess"]
    resources:
        **CFG["resources"]["preprocess"]
    group: "enormous_bam"
    wildcard_constraints:
        sample_id = "|".join(all_other_ids)
    shell:
        op.as_one_line("""
        gridss
        --reference {input.fasta}
        --workingdir $(dirname {output.workdir})
        --threads {threads}
        --jvmheap {params.mem_mb}m
        --steps {params.steps}
        {params.opts}
        {input.bam}
        2>&1 | tee -a {log}
        """)

def get_input_per_patient(wildcards):
    CFG = config['lcr-modules']['gridss']
    PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id, tumour_seq_type = wildcards.seq_type)
    if wildcards.pair_status in ["matched", "unmatched"]:
        SAMPLES = PATIENT['normal_sample_id'].unique().tolist() + PATIENT['tumour_sample_id'].tolist()
        bams = expand(
            [
                str(rules._gridss_input_bam.output.sample_bam)
            ],
            zip,
            sample_id = SAMPLES,
            allow_missing = True
        )
        preproc = expand(
            [
                str(rules._gridss_preprocess.output.workdir)
            ],
            zip,
            sample_id = SAMPLES,
            allow_missing = True
        )
    elif wildcards.pair_status == "no_normal":
        bams = expand(
            [
                str(rules._gridss_input_bam.output.sample_bam)
            ],
            zip,
            sample_id = PATIENT["tumour_sample_id"],
            allow_missing = True
        )
        preproc = expand(
            [
                str(rules._gridss_preprocess.output.workdir)
            ],
            zip,
            sample_id = PATIENT["tumour_sample_id"],
            allow_missing = True
        )
    return {'bams': ancient(bams), 'preproc': preproc}

def get_input_sample_ids(wildcards):
    CFG = config['lcr-modules']['gridss']
    PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id, tumour_seq_type = wildcards.seq_type)
    if wildcards.pair_status in ["matched", "unmatched"]:
        ids = ",".join([",".join(PATIENT['normal_sample_id'].unique().tolist()), ",".join(PATIENT['tumour_sample_id'].tolist())])
    elif wildcards.pair_status == "no_normal":
        ids = ",".join(PATIENT['tumour_sample_id'])
    return ids

# Run GRIDSS in paired mode
rule _gridss_run:
    input:
        unpack(get_input_per_patient),
        fasta = str(rules._gridss_input_references.output.genome_fa),
        fasta_img = str(rules._gridss_setup_references.output.genome_img),
        blacklist = reference_files("genomes/{genome_build}/encode/encode-blacklist.{genome_build}.bed")
    output:
        vcf = temp(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/gridss_raw.vcf.gz"),
        assembly = temp(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/assembly.bam"),
        assembly_dir = temp(directory(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/assembly.bam.gridss.working")),
        vcf_dir = temp(directory(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/gridss_raw.vcf.gz.gridss.working"))
    log: CFG["logs"]["gridss"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/gridss.log"
    params:
        ids = lambda wildcards: get_input_sample_ids(wildcards),
        opts = CFG["options"]["gridss"],
        steps = "assemble,call",
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gridss"]
    threads:
        CFG["threads"]["gridss"]
    resources:
        **CFG["resources"]["gridss"]
    group: "enormous_bam"
    shell:
        op.as_one_line("""
        gridss
        --reference {input.fasta}
        --output {output.vcf}
        --workingdir `dirname {output.vcf}`
        --assembly {output.assembly}
        --blacklist {input.blacklist}
        --threads {threads}
        --jvmheap {params.mem_mb}m
        --labels "{params.ids}"
        --steps {params.steps}
        {params.opts}
        {input.bams}
        2>&1 | tee -a {log}
        """)

# Annotate GRIDSS VCF with Repeatmasker
rule _gridss_annotate_repeatmasker:
    input:
        vcf = str(rules._gridss_run.output.vcf)
    output:
        vcf = temp(CFG["dirs"]["repeatmasker"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/gridss_repeatmasker.vcf.gz"),
        tbi = temp(CFG["dirs"]["repeatmasker"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/gridss_repeatmasker.vcf.gz.tbi")
    log: CFG["logs"]["repeatmasker"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/gridss_repeatmasker.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gridss"]
    threads:
        CFG["threads"]["repeatmasker"]
    resources:
        **CFG["resources"]["repeatmasker"]
    shell:
        op.as_one_line("""
        gridss_annotate_vcf_repeatmasker
        -o {output.vcf}
        -t {threads}
        -w $(dirname {output.vcf})
        {input.vcf}
        > {log} 2>&1
        """)

def get_split_ids(wildcards):
    CFG = config['lcr-modules']['gridss']
    if wildcards.normal_id == "None":
        return wildcards.tumour_id
    else:
        return wildcards.normal_id + "," + wildcards.tumour_id

rule _gridss_split_vcf:
    input:
        vcf = str(rules._gridss_run.output.vcf)
    output:
        vcf = temp(CFG['dirs']['repeatmasker'] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/{tumour_id}--{normal_id}--{pair_status}.gridss_split.vcf.gz"),
        tbi = temp(CFG['dirs']['repeatmasker'] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/{tumour_id}--{normal_id}--{pair_status}.gridss_split.vcf.gz.tbi")
    log: CFG["logs"]["repeatmasker"] + "{seq_type}--{genome_build}/{patient_id}--{pair_status}/{tumour_id}--{normal_id}--{pair_status}.gridss_split_vcf.log"
    params:
        ids = lambda wildcards: get_split_ids(wildcards),
    conda:
        CFG["conda_envs"]["bcftools"]
    threads: CFG['threads']['split']
    resources:
        **CFG['resources']['split']
    shell:
        op.as_one_line("""
        bcftools view -s {params.ids} -Oz -o {output.vcf} {input.vcf} 2> {log} &&
        tabix -p vcf {output.vcf}
        """)

def get_split_vcf(wildcards):
    CFG = config['lcr-modules']['gridss']
    TUMOUR = op.filter_samples(CFG['runs'], tumour_sample_id = wildcards.tumour_id, tumour_seq_type = wildcards.seq_type)
    vcf = expand(
        str(rules._gridss_split_vcf.output.vcf),
        patient_id = TUMOUR['tumour_patient_id'],
        allow_missing = True
    )
    return {'vcf': vcf}

def get_gripss_sample_id_cli(wildcards):
    CFG = config['lcr-modules']['gridss']
    TUMOUR = op.filter_samples(CFG["runs"], tumour_sample_id = wildcards.tumour_id, tumour_seq_type = wildcards.seq_type)
    if wildcards.pair_status in ["matched", "unmatched"]:
        return "-tumor " + str("".join(TUMOUR['tumour_sample_id'])) + " -reference " + str("".join(TUMOUR['normal_sample_id']))
    elif wildcards.pair_status == "no_normal":
        return "-tumor " + str("".join(TUMOUR['tumour_sample_id']))

# Perform somatic filtering against the panel of normals
rule _gridss_run_gripss:
    input:
        unpack(get_split_vcf),
        fasta = str(rules._gridss_input_references.output.genome_fa),
        pon_breakend = str(rules._gridss_get_pon.output.pon_breakend),
        pon_breakpoint = str(rules._gridss_get_pon.output.pon_breakpoint),
        known_pairs = str(rules._gridss_get_pon.output.known_pairs)
    output:
        vcf = CFG["dirs"]["gripss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz",
        tbi = CFG["dirs"]["gripss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz.tbi"
    log: log = CFG["logs"]["gripss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gripss.log"
    resources:
        **CFG["resources"]["gripss"]
    params:
        cli = lambda wildcards: get_gripss_sample_id_cli(wildcards),
        opts = CFG["options"]["gripss"],
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gripss"]
    threads:
        CFG["threads"]["gripss"]
    shell:
        op.as_one_line("""
        gripss -Xms4G -Xmx{params.mem_mb}m
        -ref_genome {input.fasta}
        -breakend_pon {input.pon_breakend}
        -breakpoint_pon {input.pon_breakpoint}
        -breakpoint_hotspot {input.known_pairs}
        -input_vcf {input.vcf}
        -output_vcf {output.vcf}
        {params.cli}
        {params.opts}
        2>&1 | tee -a {log}
        """)

rule _gridss_filter_gripss:
    input:
        vcf = str(rules._gridss_run_gripss.output.vcf)
    output:
        vcf = CFG["dirs"]["gripss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.vcf.gz",
        tbi = CFG["dirs"]["gripss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        zcat {input.vcf} |
            awk '$7 == "PASS" || $1 ~ /^#/ ' |
            bcftools view -Oz -o {output.vcf} &&
        tabix -p vcf {output.vcf}
        """)

rule _gridss_gripss_to_bedpe:
    input:
        vcf = str(rules._gridss_filter_gripss.output.vcf)
    output:
        bedpe = CFG["dirs"]["gripss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.bedpe"
    conda:
        CFG["conda_envs"]["svtools"]
    shell:
        op.as_one_line("""
        if zcat {input.vcf} |  awk '$1 ~ /^#/ || $5 ~ /:/' | tail -1  | grep -q "^#CHROM";
        then
            touch {output.bedpe};
        else
            zcat {input.vcf} |
                awk '$1 ~ /^#/ || $5 ~ /:/' |
                svtools vcftobedpe | grep -v "##" > {output.bedpe};
        fi
        """)


# Symlink the final output files into the module results directory (under '99-outputs/')
rule _gridss_output_somatic_vcf:
    input:
        filtered = str(rules._gridss_filter_gripss.output.vcf),
        filtered_tbi = str(rules._gridss_filter_gripss.output.tbi),
        somatic = str(rules._gridss_run_gripss.output.vcf),
        somatic_tbi = str(rules._gridss_run_gripss.output.tbi),
        bedpe = str(rules._gridss_gripss_to_bedpe.output.bedpe)
    output:
        somatic = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_somatic.vcf.gz",
        somatic_tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_somatic.vcf.gz.tbi",
        filtered = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_somatic_filtered.vcf.gz",
        filtered_tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_somatic_filtered.vcf.gz.tbi",
        bedpe = CFG["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_somatic_filtered.bedpe"
    run:
        op.relative_symlink(input.somatic, output.somatic, in_module=True)
        op.relative_symlink(input.somatic_tbi, output.somatic_tbi, in_module=True)
        op.relative_symlink(input.filtered, output.filtered, in_module=True)
        op.relative_symlink(input.filtered_tbi, output.filtered_tbi, in_module=True)
        op.relative_symlink(input.bedpe, output.bedpe, in_module=True)




# Generates the target sentinels for each run, which generate the symlinks
rule _gridss_all:
    input:
        expand(
            [
                str(rules._gridss_output_somatic_vcf.output.filtered),
                str(rules._gridss_output_somatic_vcf.output.filtered_tbi),
                str(rules._gridss_output_somatic_vcf.output.somatic),
                str(rules._gridss_output_somatic_vcf.output.somatic_tbi),
                str(rules._gridss_output_somatic_vcf.output.bedpe)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
