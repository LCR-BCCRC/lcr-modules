#!/usr/bin/env snakemake


##### SETUP #####
import modutils as md
import glob

# Make sure the `CFG` variable doesn"t exist yet
assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

# Setup module and store module-specific configuration in `CFG`.
CFG = md.setup_module(
    config = config, 
    name = "cellranger", 
    version = "1.0",
    subdirs = ["inputs", "samplesheet", "mkfastq", "count", "vdj", "outputs"],
    req_references = ["transcriptome"]
)

localrules: 
    _cellranger_input,
    _cellranger_all


##### RULES #####
rule _cellranger_input:
    input:
        lib = CFG["inputs"]["sample_dir"]
    output:
        lib = directory(CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{chip_id}")
    run:
        f = glob.glob(input.lib + "*" + wildcards.chip_id + "*")
        md.symlink(f[0], output.lib)

"""
rule _cellranger_create_samplesheet:
    output:
        ss = CFG["dirs"]["samplesheet"] + "{chip_id}_samplesheet.csv"
    run:
        df = md.filter_samples(CFG["samples"], chip_id == wildcards.chip_id)
        ss = df[["lane", "sample_id", "index"]]
        ss.columns = ["lane", "sample", "index"]
        ss.to_csv({output.ss}, sep = ",", index = False)

"""
def _create_samplesheet(mconfig = CFG):
    samples = mconfig["samples"]
    ss_dir = mconfig["dirs"]["samplesheet"]
    def _custom_samplesheet(wildcards):
        df = samples[samples["chip_id"] == wildcards.chip_id]
        ss = df[["lane", "sample_id", "index"]]
        ss.columns = ["lane", "sample", "index"]
        ss_file = ( ss_dir + f"{wildcards.chip_id}_samplesheet.csv")
        ss.to_csv(ss_file, sep = ",", index = False)
        return ss_file
    return _custom_samplesheet


def _get_completion_files(raw_dir = CFG["inputs"]["sample_dir"], suffix = ["RTAComplete*", "RunInfo*", "RunParameters*"]):
    def _get_custom_files(wildcards):
        #path = raw_dir + f"{wildcards.seq_type}--{wildcards.genome_build}/*{wildcards.chip_id}*"
        runs = glob.glob(raw_dir + f"*{wildcards.chip_id}*")[0]
        file = []
        for f in suffix:
            file.append(glob.glob(runs + "/" + f)[0])
        return file
    return _get_custom_files


rule _cellranger_mkfastq:
    input:
        run_dir = rules._cellranger_input.output.lib,
        ss = _create_samplesheet(), #rules._cellranger_create_samplesheet.output.ss,
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
#    conda:
#        CFG["conda_envs"].get("cellranger") or "envs/cellranger-3.0.2.yaml"
    threads:
        CFG["threads"].get("mkfastq") or 8
    resources: 
        mem_mb = CFG["mem_mb"].get("mkfastq") or 50000
    shell:
        md.as_one_line("""
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
        stamp = rules._cellranger_mkfastq.output.stamp
    output:
        stamp = CFG["dirs"]["outputs"] + "stamps/{seq_type}--{genome_build}/{chip_id}_{sample_id}_count.stamp"
    log:
        stdout = CFG["logs"]["count"] + "{seq_type}--{genome_build}/{chip_id}_{sample_id}_count.stdout.log",
        stderr = CFG["logs"]["count"] + "{seq_type}--{genome_build}/{chip_id}_{sample_id}_count.stderr.log"
    params:
        cr = CFG["software"],
        fastq_dir = CFG["dirs"]["mkfastq"] + "{seq_type}--{genome_build}/chip_{chip_id}/",
        opts = CFG["options"]["count"],
        ref = md.get_reference(CFG, "transcriptome")
#    conda:
#        CFG["conda_envs"].get("cellranger") or "envs/cellranger-3.0.2.yaml"
    threads:
        CFG["threads"].get("count") or 8
    resources: 
        mem_mb = CFG["mem_mb"].get("count") or 50000
    shell:
        md.as_one_line("""
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
        stamp = rules._cellranger_mkfastq.output.stamp,
        fastq_dir = CFG["dirs"]["mkfastq"] + "{seq_type}--{genome_build}/{chip_id}"
    output:
        stamp = CFG["dirs"]["outputs"] + "stamps/{seq_type}--{genome_build}/{chip_id}_{sample_id}_vdj.stamp"
    log:
        stdout = CFG["logs"]["vdj"] + "{seq_type}--{genome_build}/{chip_id}_{sample_id}_vdj.stdout.log",
        stderr = CFG["logs"]["vdj"] + "{seq_type}--{genome_build}/{chip_id}_{sample_id}_vdj.stderr.log"
    params:
        cr = CFG["software"],
        opts = CFG["options"]["vdj"],
        ref = md.get_reference(CFG, "vdj")
#    conda:
#        CFG["conda_envs"].get("cellranger") or "envs/cellranger-3.0.2.yaml"
    threads:
        CFG["threads"].get("vdj") or 8
    resources: 
        mem_gb = CFG["mem_mb"].get("vdj") or 50000
    shell:
        md.as_one_line("""
        {params.cr} vdj
        {params.opts}
        --id={wildcards.sample_id}
        --sample={wildcards.sample_id}
        --fastqs={input.fastq_dir} 
        --reference={params.ref}
        --localcores={threads}
        --localmem=$(({resources.mem_mb}/1000))
        > {log.stdout} 2> {log.stderr}
        && touch {output.stamp}
        """)


# Generates the target sentinels for each run, which generate the symlinks
# TODO: Update to ask for the output of every `_cellranger_output_*` rule
rule _cellranger_all:
    input:
        expand(expand("{{dir}}stamps/{seq_type}--{genome_build}/{chip_id}_{sample_id}_{analysis}.stamp",
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
md.cleanup_module(CFG)

# Delete any local variables to avoid interfering with other code
del CFG
