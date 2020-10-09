#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op


# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["gridss"]`
CFG = op.setup_module(
    name = "gridss",
    version = "1.0",
    subdirectories = ["inputs", "preprocess", "gridss", "viral_annotation", "gripss", "outputs"],
)

possible_genome_builds = ["grch37", "hs37d5", "hg38"]
for genome_build in CFG["runs"]["tumour_genome_build"]:
    assert genome_build in possible_genome_builds, (
        "Samples table includes genome builds not yet compatible with this module. "
        "This module is currently only compatible with {possible_genome_builds}. "
    )

sample_ids = list(CFG['samples']['sample_id'])
unmatched_normal_ids = list(config["lcr-modules"]["_shared"]["unmatched_normal_ids"].values())
all_other_ids = list(set(sample_ids) - set(unmatched_normal_ids))

# Define rules to be run locally when using a compute cluster
localrules:
    _gridss_input_bam,
    _gridss_input_references,
    _gridss_setup_references,
    _gridss_input_viral_ref,
    _gridss_setup_viral_ref,
    _gridss_symlink_preprocessed_normal, 
    _gridss_unpaired_filter,
    _gridss_filter_gripss,
    _gridss_unpaired_to_bedpe, 
    _gridss_gripss_to_bedpe, 
    _gridss_output_viral_vcf, 
    _gridss_output_somatic_vcf,
    _gridss_dispatch,
    _gridss_all

wildcard_constraints: 
    genome_build = "|".join(possible_genome_builds),
    pair_status = "matched|unmatched|no_normal"


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

# Symlink genome fasta with bwa and .fai indices to the same directory
rule _gridss_input_viral_ref: 
    input: 
        viral_fa = CFG["references"]["viral_fa"],
        viral_bwa_prefix = CFG["references"]["viral_bwa_prefix"]
    output: 
        viral_fa = CFG["dirs"]["inputs"] + "references/human_virus/human_virus.fa"
    shell: 
        op.as_one_line("""
        ln -sf {input.viral_fa} {output.viral_fa} &&
        ln -sf {input.viral_fa}.fai {output.viral_fa}.fai &&
        ln -sf {input.viral_bwa_prefix}.* `dirname {output.viral_fa}`
        """)

# Generage human_virus.fa.img file
rule _gridss_setup_viral_ref: 
    input: 
        fasta = str(rules._gridss_input_viral_ref.output)
    output: 
        viral_img = CFG["dirs"]["inputs"] + "references/human_virus/human_virus.fa.img"
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
        --workingdir `dirname {output.viral_img}` 
        """)


# Symlink the input files into the module results directory (under '00-inputs/')
rule _gridss_input_bam:
    input:
        sample_bam = CFG["inputs"]["sample_bam"], 
        sample_bai = CFG["inputs"]["sample_bai"] 
    output:
        sample_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam", 
        sample_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai" 
    run:
        op.relative_symlink(input.sample_bam, output.sample_bam)
        op.relative_symlink(input.sample_bai, output.sample_bai)

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
        # 'touch {output.workdir} && sleep 15'
        op.as_one_line("""
        function samtools() {{
            if [[ "$1" == "sort" ]]; then
                shift 1;
                command samtools sort -m {params.mem_mb}M "$@";
            else
                command samtools "$@";
            fi;
        }} && 
        export -f samtools &&
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
        workdir = temp(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{sample_id}.bam.gridss.working")
        # workdir = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{sample_id}.bam.gridss.working"
    priority: 0
    wildcard_constraints: 
        sample_id = "|".join(unmatched_normal_ids)
    run: 
        op.relative_symlink(input.workdir, output.workdir)

# Preprocess all other bams as part of the group job
rule _gridss_preprocess:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = str(rules._gridss_input_references.output.genome_fa),
        fasta_img = str(rules._gridss_setup_references.output.genome_img)
    output:
        workdir = temp(directory(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{sample_id}.bam.gridss.working"))
    log: CFG["logs"]["preprocess"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{sample_id}/preprocess.log"
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
    group: "enormous_bam"
    wildcard_constraints: 
        sample_id = "|".join(all_other_ids)
    shell:
        # 'touch {output.workdir} && sleep 15'
        op.as_one_line("""
        function samtools() {{
            if [[ "$1" == "sort" ]]; then
                shift 1;
                command samtools sort -m {params.mem_mb}M "$@";
            else
                command samtools "$@";
            fi
        }} &&
        export -f samtools &&
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




# Run GRIDSS in paired mode
rule _gridss_paired:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        tumour_preproc = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam.gridss.working", 
        normal_preproc = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{normal_id}.bam.gridss.working",
        fasta = str(rules._gridss_input_references.output.genome_fa),
        fasta_img = str(rules._gridss_setup_references.output.genome_img), 
        blacklist = reference_files("genomes/{genome_build}/encode/encode-blacklist.{genome_build}.bed"), 
        repeatmasker = reference_files("genomes/{genome_build}/repeatmasker/repeatmasker.{genome_build}.bed")
    output:
        vcf = temp(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_raw.vcf.gz"),
        assembly = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/assembly.bam", 
        assembly_dir = temp(directory(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/assembly.bam.gridss.working")), 
        vcf_dir = temp(directory(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_raw.vcf.gz.gridss.working"))
    log: CFG["logs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss.log"
    params:
        opts = CFG["options"]["gridss"], 
        steps = "assemble,call", 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8) 
    conda:
        CFG["conda_envs"]["gridss"]
    threads:
        CFG["threads"]["gridss"]
    resources:
        **CFG["resources"]["gridss"]
    wildcard_constraints: 
        pair_status = "matched|unmatched"
    group: "enormous_bam"
    shell:
        # 'touch {output} && sleep 15'
        op.as_one_line("""
        function samtools() {{
            if [[ "$1" == "sort" ]]; then
                shift 1;
                command samtools sort -m {params.mem_mb}M "$@";
            else
                command samtools "$@";
            fi
        }} &&
        export -f samtools &&
        gridss
        --reference {input.fasta}
        --output {output.vcf}
        --workingdir `dirname {output.vcf}`
        --assembly {output.assembly}
        --blacklist {input.blacklist}
        --repeatmaskerbed {input.repeatmasker}
        --threads {threads}
        --jvmheap {params.mem_mb}m
        --labels "{wildcards.normal_id},{wildcards.tumour_id}"
        --steps {params.steps}
        {params.opts}
        {input.normal_bam} 
        {input.tumour_bam} 
        2>&1 | tee -a {log}
        """)
   
# Run GRIDSS in unpaired mode
rule _gridss_unpaired:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        tumour_preproc = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam.gridss.working",
        fasta = str(rules._gridss_input_references.output.genome_fa),
        fasta_img = str(rules._gridss_setup_references.output.genome_img), 
        blacklist = reference_files("genomes/{genome_build}/encode/encode-blacklist.{genome_build}.bed"), 
        repeatmasker = reference_files("genomes/{genome_build}/repeatmasker/repeatmasker.{genome_build}.bed")
    output:
        vcf = temp(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--None--no_normal/gridss_raw.vcf.gz"),
        assembly = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--None--no_normal/assembly.bam", 
        assembly_dir = temp(directory(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--None--no_normal/assembly.bam.gridss.working")), 
        vcf_dir = temp(directory(CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--None--no_normal/gridss_raw.vcf.gz.gridss.working"))
    log: CFG["logs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--None--no_normal/gridss.log"
    params:
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
    wildcard_constraints: 
        normal_id = "None",
        pair_status = "no_normal"
    shell:
        # 'touch {output} && sleep 15'
        op.as_one_line("""
        function samtools() {{
            if [[ "$1" == "sort" ]]; then
                shift 1;
                command samtools sort -m {params.mem_mb}M "$@";
            else
                command samtools "$@";
            fi
        }} &&
        export -f samtools && 
        gridss
        --reference {input.fasta}
        --output {output.vcf}
        --workingdir `dirname {output.vcf}`
        --assembly {output.assembly}
        --blacklist {input.blacklist}
        --repeatmaskerbed {input.repeatmasker}
        --threads {threads}
        --jvmheap {params.mem_mb}m
        --labels "{wildcards.tumour_id}"
        --steps {params.steps}
        {params.opts}
        {input.tumour_bam} 
        2>&1 | tee -a {log}
        """)

# Perform viral annotation of the output VCFs
rule _gridss_viral_annotation: 
    input: 
        vcf = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_raw.vcf.gz", 
        viral_ref = str(rules._gridss_input_viral_ref.output), 
        viral_img = str(rules._gridss_setup_viral_ref.output)
    output: 
        vcf = temp(CFG["dirs"]["viral_annotation"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_viral_annotation.vcf.gz")
    log: CFG["logs"]["viral_annotation"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_viral_annotation.log"
    resources: 
        **CFG["resources"]["viral_annotation"]
    params: 
        gridss_jar = "$(readlink -e $(which gridss)).jar", 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: 
        CFG["conda_envs"]["gridss"]
    threads: 
        CFG["threads"]["viral_annotation"]
    shell: 
        # 'touch {output.vcf} && sleep 15'
        op.as_one_line("""
        java -Xmx{params.mem_mb}m 
                -cp {params.gridss_jar} gridss.AnnotateInsertedSequence 
				REFERENCE_SEQUENCE={input.viral_ref} 
				INPUT={input.vcf} 
				OUTPUT={output.vcf} 
				WORKER_THREADS={threads} 
                2>&1 | tee -a {log}
        """)

# Filter unpaired VCFs and output to bedpe. 
rule _gridss_unpaired_filter: 
    input: 
        vcf = CFG["dirs"]["viral_annotation"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_viral_annotation.vcf.gz"
    output: 
        vcf = CFG["dirs"]["viral_annotation"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_viral_annotation_filtered.vcf.gz", 
        tbi = CFG["dirs"]["viral_annotation"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_viral_annotation_filtered.vcf.gz.tbi"
    conda: 
        CFG["conda_envs"]["bcftools"]
    wildcard_constraints: 
        normal_id = "None", 
        pair_status = "no_normal"
    shell: 
        # 'touch {output} && sleep 15'
        op.as_one_line("""
        zcat {input.vcf} | 
            awk '($5 ~ /:/ && $7 == "PASS") || $1 ~ /^#/' | 
            bcftools view -Oz -o {output.vcf} && 
            tabix -p vcf {output.vcf}
        """)

rule _gridss_unpaired_to_bedpe: 
    input: 
        vcf = str(rules._gridss_unpaired_filter.output.vcf)
    output: 
        bedpe = CFG["dirs"]["viral_annotation"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_viral_annotation_filtered.bedpe"
    conda:
        CFG["conda_envs"]["svtools"]
    shell: 
        # 'touch {output.bedpe} && sleep 15'
        op.as_one_line("""
        svtools vcftobedpe -i {input.vcf} -o {output.bedpe}
        """) 



# Perform somatic filtering against the panel of normals    
rule _gridss_run_gripss: 
    input: 
        vcf = str(rules._gridss_viral_annotation.output.vcf), 
        fasta = str(rules._gridss_input_references.output.genome_fa),
    output: 
        vcf = CFG["dirs"]["gripss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz", 
        tbi = CFG["dirs"]["gripss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz.tbi"
    log: log = CFG["logs"]["gripss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gripss.log"
    resources: 
        **CFG["resources"]["gripss"]
    params:
        pon_dir = CFG["references"]["pon_dir"],
        alt_build = lambda w: {
            "grch37": "hg19", 
            "hs37d5": "hg19", 
            "hg38": "hg38"}[w.genome_build], 
        opts = CFG["options"]["gripss"], 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: 
        CFG["conda_envs"]["gripss"]
    threads: 
        CFG["threads"]["gripss"]
    shell: 
        # 'touch {output} && sleep 15'
        op.as_one_line(""" 
        gripss -Xms4G -Xmx{params.mem_mb}m 
        -ref_genome {input.fasta} 
        -breakend_pon {params.pon_dir}/gridss_pon_single_breakend.{params.alt_build}.bed 
        -breakpoint_pon {params.pon_dir}/gridss_pon_breakpoint.{params.alt_build}.bedpe 
        -breakpoint_hotspot {params.pon_dir}/KnownFusionPairs.{params.alt_build}.bedpe 
        -input_vcf {input.vcf} 
        -output_vcf {output.vcf} 
        -tumor {wildcards.tumour_id} 
        -reference {wildcards.normal_id}  
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
        # 'touch {output} && sleep 15'
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
        # 'touch {output} && sleep 15'
        op.as_one_line("""
        zcat {input.vcf} | 
            awk '$1 ~ /^#/ || $5 ~ /:/' | 
            svtools vcftobedpe -o {output.bedpe}
        """)
     

# Symlink the final output files into the module results directory (under '99-outputs/')
rule _gridss_output_viral_vcf:
    input:
        vcf = str(rules._gridss_unpaired_filter.output.vcf), 
        tbi = str(rules._gridss_unpaired_filter.output.tbi), 
        bedpe = str(rules._gridss_unpaired_to_bedpe.output.bedpe)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_viral_annotation_filtered.vcf.gz", 
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_viral_annotation_filtered.vcf.gz.tbi", 
        bedpe = CFG["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_viral_annotation_filtered.bedpe"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.tbi, output.tbi)
        op.relative_symlink(input.bedpe, output.bedpe)

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
        op.relative_symlink(input.somatic, output.somatic)
        op.relative_symlink(input.somatic_tbi, output.somatic_tbi)
        op.relative_symlink(input.filtered, output.filtered)
        op.relative_symlink(input.filtered_tbi, output.filtered_tbi)
        op.relative_symlink(input.bedpe, output.bedpe)

def _gridss_predict_output(wildcards): 
    """Request symlinks for all VCF files.
    
    This function requests symlinks depending on whether samples have 
    been run in paired mode with gripss somatic filtering or in 
    unpaired mode without gripss somatic filtering. 
    """
    
    CFG = config["lcr-modules"]["gridss"]

    if wildcards.pair_status == "matched" or wildcards.pair_status == "unmatched": 
        vcf_names = ["gridss_somatic", "gridss_somatic_filtered"]
        output_vcf = expand(str(CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz"), 
            vcf_name = vcf_names, 
            **wildcards)

    else: 
        vcf_names = ["gridss_viral_annotation_filtered"]
        output_vcf = expand(str(CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz"), 
            vcf_name = vcf_names, 
            **wildcards)

    return output_vcf 

# Dispatch rule to return all symlinked files

rule _gridss_dispatch: 
    input: 
        _gridss_predict_output
    output: 
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _gridss_all:
    input:
        expand(
            [
                str(rules._gridss_dispatch.output.dispatched),
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
