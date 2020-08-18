#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Helena Winata
# Module Author:    Helena Winata


##### SETUP #####
# Import standard modules
import glob

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`.
# `CFG` is a shortcut to `config["lcr-modules"]["cellranger"]`
CFG = op.setup_module(
    name = "cellranger", 
    version = "1.0",
    subdirectories = ["inputs", "samplesheet", "mkfastq", "count", "vdj", "outputs"],
)

localrules: 
    _cellranger_input,
    _cellranger_create_samplesheet,
    _cellranger_all


##### RULES #####
rule _cellranger_input:
    input:
        lib = CFG["inputs"]["sample_dir"]
    output:
        lib = directory(CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{chip_id}")
    run:
        f = glob.glob(input.lib + "*" + wildcards.chip_id + "*")
        print(f)
        op.symlink(f[0], output.lib)


rule _cellranger_create_samplesheet:
    output:
        ss = CFG["dirs"]["samplesheet"] + "{chip_id}_samplesheet.csv"
    run:
        df = op.filter_samples(CFG["samples"], chip_id = wildcards.chip_id)
        ss = df[["lane", "sample_id", "index"]]
        ss.columns = ["lane", "sample", "index"]
        ss.to_csv(output.ss, sep = ",", index = False)



def _get_completion_files(raw_dir = CFG["inputs"]["sample_dir"], suffix = ["RTAComplete*", "RunInfo*", "RunParameters*"]):
    def _get_custom_files(wildcards):
        runs = glob.glob(raw_dir + f"/*{wildcards.chip_id}*")[0]
        file = []
        for f in suffix:
            file.append(glob.glob(runs + "/" + f)[0])
        return file
    return _get_custom_files


rule _cellranger_mkfastq:
    input:
        run_dir = str(rules._cellranger_input.output.lib),
        ss = str(rules._cellranger_create_samplesheet.output.ss),
        check = _get_completion_files()
    output:
        stamp = CFG["dirs"]["outputs"] + "stamps/{seq_type}--{genome_build}/{chip_id}_mkfastq.stamp"
    log:
        stdout = CFG["logs"]["mkfastq"] + "{seq_type}--{genome_build}/{chip_id}/mkfastq.stdout.log",
        stderr = CFG["logs"]["mkfastq"] + "{seq_type}--{genome_build}/{chip_id}/mkfastq.stderr.log"
    params:
        cr = CFG["software"],
        out_dir = CFG["dirs"]["mkfastq"] + "{seq_type}--{genome_build}/chip_{chip_id}",
        opts = CFG["options"]["mkfastq"]
    threads:
        CFG["threads"]["mkfastq"]
    resources: 
        mem_mb = CFG["mem_mb"]["mkfastq"]
    shell:
        op.as_one_line("""
        {params.cr} mkfastq
        {params.opts}
        --run={input.run_dir} 
        --samplesheet={input.ss}
        --output-dir={params.out_dir}
        --localcores={threads}
        --localmem=$(({resources.mem_mb}/1000))
        > {log.stdout} 2> {log.stderr}
        && touch {output.stamp}
        """)


rule _cellranger_count:
    input:
        stamp = str(rules._cellranger_mkfastq.output.stamp)
    output:
        stamp = CFG["dirs"]["outputs"] + "stamps/{seq_type}--{genome_build}/{chip_id}--{sample_id}_count.stamp"
    log:
        stdout = CFG["logs"]["count"] + "{seq_type}--{genome_build}/{chip_id}--{sample_id}_count.stdout.log",
        stderr = CFG["logs"]["count"] + "{seq_type}--{genome_build}/{chip_id}--{sample_id}_count.stderr.log"
    params:
        cr = CFG["software"],
        fastq_dir = CFG["dirs"]["mkfastq"] + "{seq_type}--{genome_build}/chip_{chip_id}/",
        opts = CFG["options"]["count"],
        ref = CFG["reference"]["transcriptome"]
    threads:
        CFG["threads"]["count"]
    resources: 
        mem_mb = CFG["mem_mb"]["count"]
    shell:
        op.as_one_line("""
        {params.cr} count
        {params.opts}
        --id={wildcards.sample_id}
        --sample={wildcards.sample_id}
        --fastqs={params.fastq_dir} 
        --transcriptome={params.ref}
        --localcores={threads}
        --localmem=$(({resources.mem_mb}/1000))
        > {log.stdout} 2> {log.stderr}
        && touch {output.stamp}
        """)


rule _cellranger_vdj:
    input:
        stamp = str(rules._cellranger_mkfastq.output.stamp)
    output:
        stamp = CFG["dirs"]["outputs"] + "stamps/{seq_type}--{genome_build}/{chip_id}--{sample_id}_vdj.stamp"
    log:
        stdout = CFG["logs"]["vdj"] + "{seq_type}--{genome_build}/{chip_id}--{sample_id}_vdj.stdout.log",
        stderr = CFG["logs"]["vdj"] + "{seq_type}--{genome_build}/{chip_id}--{sample_id}_vdj.stderr.log"
    params:
        cr = CFG["software"],
        fastq_dir = CFG["dirs"]["mkfastq"] + "{seq_type}--{genome_build}/chip_{chip_id}",
        opts = CFG["options"]["vdj"],
        ref = CFG["reference"]["vdj"]
    threads:
        CFG["threads"]["vdj"]
    resources: 
        mem_gb = CFG["mem_mb"]["vdj"]
    shell:
        op.as_one_line("""
        {params.cr} vdj
        {params.opts}
        --id={wildcards.sample_id}
        --sample={wildcards.sample_id}
        --fastqs={params.fastq_dir} 
        --reference={params.ref}
        --localcores={threads}
        --localmem=$(({resources.mem_gb}/1000))
        > {log.stdout} 2> {log.stderr}
        && touch {output.stamp}
        """)


# Generates the target sentinels for each run, which generate the symlinks
rule _cellranger_all:
    input:
        expand(expand("{{dir}}stamps/{seq_type}--{genome_build}/{chip_id}--{sample_id}_{analysis}.stamp",
               zip,  
               seq_type=CFG["samples"]["seq_type"],
               genome_build=CFG["samples"]["genome_build"],
               chip_id = CFG["samples"]["chip_id"],
               sample_id=CFG["samples"]["sample_id"],
               analysis = CFG["samples"]["analysis"]),
               dir = CFG["dirs"]["outputs"])

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk (including the samples and runs)
op.cleanup_module(CFG)
