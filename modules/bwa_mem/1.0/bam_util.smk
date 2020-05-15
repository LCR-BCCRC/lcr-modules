##### RULES #####

rule _bam_util_sort:
    input:
        bam = CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.bam"
    output:
        bam = temp(CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.sorted.bam")
    priority: 10
    log:
        stderr = CFG["logs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.BAMsort.stderr.log"
    params:
        opts = CFG["options"]["sort"],
        dir = CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}"
    conda:
        CFG["conda_envs"].get("samtools") or "envs/samtools-1.9.yaml"
    threads:
        CFG["threads"].get("sort") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("sort") or 5000
    shell:
        md.as_one_line("""
        samtools sort {params.opts}
        -T {params.dir} -o {output.bam} {input.bam} 2> {log.stderr}
        """)


rule _bam_util_filter:
    input:
        bam = CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.bam"
    output:
        bam = temp(CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.filtered.bam")
    priority: 10
    log:
        stderr = CFG["logs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.BAMfilter.stderr.log"
    params:
        opts = CFG["options"]["filter"]
    conda:
        CFG["conda_envs"].get("samtools") or "envs/samtools-1.9.yaml"
    threads:
        CFG["threads"].get("filter") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("filter ") or 5000
    shell:
        md.as_one_line("""
        samtools view {params.opts} {input.bam}
        > {output.bam} 2> {log.stderr}
        """)

rule _bam_util_markdup:
    input:
        bam = CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.bam"
    output:
        bam = temp(CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.markdup.bam"),
        metrics = CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.dup_metrics.bam"
    priority: 10
    log:
        stdout = CFG["logs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.BAMmarkdup.stdout.log",
        stderr = CFG["logs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.BAMmarkdup.stderr.log"
    params:
        opts = CFG["options"]["markdup"]
    conda:
        CFG["conda_envs"].get("bwa") or "envs/bwa-0.7.17.yaml"
    threads:
        CFG["threads"].get("markdup") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("markdup") or 8000
    shell:
        md.as_one_line("""
        bammarkduplicates2 {params.opts}
        I={input.bam} O={output.bam} M={output.metrics}
        > {log.stdout} 2> {log.stderr}
        """)


rule _bam_util_rmdup:
    input:
        bam = CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.bam"
    output:
        bam = temp(CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.rmdup.bam")
    priority: 10
    log:
         stderr = CFG["logs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_name}.BAMrmdup.stderr.log"
    params:
        opts = CFG["options"]["rmdup"]
    conda:
        CFG["conda_envs"].get("samtools") or "envs/samtools-1.9.yaml"
    threads:
        CFG["threads"].get("rmdup") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("rmdup") or 5000
    shell:
        md.as_one_line("""
        samtools rmdup {params.opts}
        {input.bam} {output.bam} 2> {log.stderr}
        """)

'''
rule _bam_util_output:
    input:
        bam = CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_id}." + CFG["suffix"] + ".bam",
    output:
        bam = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bam",
    shell:
        "mv {input} {output}"
#   run:
#       md.symlink(input.bam, output.bam)
'''

rule _bam_util_index:
    input:
        bam = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    output:
        bai = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bam.bai"
    log:
        stderr = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.BAMindex.stderr.log"
    conda:
        CFG["conda_envs"].get("samtools") or "envs/samtools-1.9.yaml"
    threads:
        CFG["threads"].get("filter") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("filter ") or 5000
    shell:
        "samtools index {input.bam} {output.bai} 2> {log.stderr}"

'''
rule _bam_util_all:
    input:
        bai = expand(rules._bam_util_index.output.bai, zip,
                    seq_type = list(CFG["samples"]['seq_type']),
                    genome_build = list(CFG["samples"]['genome_build']),
                    sample_id = list(CFG["samples"]['sample_id']))

##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
'''