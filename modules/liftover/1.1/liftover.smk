#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op


# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["liftover"]`
CFG = op.setup_module(
    name = "liftover",
    version = "1.1",
    subdirectories = ["inputs", "seg2bed", "liftover", "bed2seg", "outputs"])

# Define rules to be run locally when using a compute cluster
localrules:
    _liftover_input_seg,
    _liftover_seg_2_bed,
    _lofreq_all,
    _run_liftover,
    _liftover_sort,
    _liftover_bed_2_seg,
    _liftover_output_seg,
    _liftover_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _liftover_input_seg:
    input:
        seg = CFG["inputs"]["sample_seg"]
    output:
        seg = CFG["dirs"]["inputs"] + "{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.seg"
    run:
        op.relative_symlink(input.seg, output.seg)


# Convert initial seg file into bed format
rule _liftover_seg_2_bed:
    input:
        seg = str(rules._liftover_input_seg.output.seg)
    output:
        bed = CFG["dirs"]["seg2bed"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.bed",
        header = temp(CFG["dirs"]["seg2bed"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.bed.header")
    log:
        stderr = CFG["logs"]["seg2bed"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.stderr.log"
    params:
        opts = CFG["options"]["seg2bed2seg"],
        chr_colNum = CFG["options"]["chr_colNum"],
        start_colNum = CFG["options"]["start_colNum"],
        end_colNum = CFG["options"]["end_colNum"],
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.seg} 
        --output {output.bed} 
        --chromColnum {params.chr_colNum} 
        --startColnum {params.start_colNum} 
        --endColnum {params.end_colNum}
        {params.opts} 
        2> {log.stderr}
        """)


def get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")


# Convert the bed file in hg38 coordinates into hg19 coordinates
rule _run_liftover:
    input:
        native = rules._liftover_seg_2_bed.output.bed,
        chains = get_chain
    output:
        lifted = temp(CFG["dirs"]["liftover"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.lifted_{chain}.bed"),
        unmapped = CFG["dirs"]["liftover"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.lifted_{chain}.unmapped.bed"
    log:
        stderr = CFG["logs"]["liftover"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.lifted_{chain}.stderr.log"
    params:
        opts = CFG["options"]["liftover"],
        mismatch = CFG["options"]["min_mismatch"]
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        liftOver -minMatch={params.mismatch}
        {input.native} {input.chains} 
        {output.lifted} {output.unmapped}
        2> {log.stderr}
        """)

# Sort liftover output
rule _liftover_sort:
    input:
        lifted = expand(CFG["dirs"]["liftover"] + "from--{{genome_build}}/{{tumour_sample_id}}--{{normal_sample_id}}.{{tool}}.lifted_{chain}.bed",
                  chain="hg38ToHg19" if "38" in str({genome_build}) else "hg19ToHg38")
    output:
        lifted_sorted = CFG["dirs"]["liftover"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.lifted_{chain}.sorted.bed"
    log:
        stderr = CFG["logs"]["liftover"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.lifted_{chain}.sorted.stderr.log"
    params:
        opts = CFG["options"]["liftover"],
        mismatch = CFG["options"]["min_mismatch"]       
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        bedtools sort -i {input.lifted} > {output.lifted_sorted}
        2> {log.stderr}
        """)


# Convert the bed file in hg19 coordinates into seg format
rule _liftover_bed_2_seg:
    input:
        lifted_sorted = rules._liftover_sort.output.lifted_sorted,
        headers = rules._liftover_seg_2_bed.output.header
    output:
        seg_lifted = CFG["dirs"]["bed2seg"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["bed2seg"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.lifted_{chain}.stderr.log"
    params:
        opts = CFG["options"]["seg2bed2seg"]      
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.lifted_sorted}
        --column-header {input.headers}
        --output {output.seg_lifted} 
        {params.opts} 
        2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _liftover_output_seg:
    input:
        seg = rules._liftover_bed_2_seg.output.seg_lifted
    output:
        seg = CFG["dirs"]["outputs"] + "from--{genome_build}/{tumour_sample_id}--{normal_sample_id}.{tool}.lifted_{chain}.seg"
    run:
        op.relative_symlink(input.seg, output.seg)


# Generates the target sentinels for each run, which generate the symlinks
rule _liftover_all:
    input:
        expand(
            [
                rules._liftover_output_seg.output.seg
            ],
            zip,  # Run expand() with zip(), not product()
            tumour_sample_id=CFG["runs"]["tumour_sample_id"],
            normal_sample_id=CFG["runs"]["normal_sample_id"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            #repeat the tool name N times in expand so each pair in run is used
            tool=[CFG["tool"]] * len(CFG["runs"]["tumour_sample_id"]),
            chain="hg38ToHg19" if "38" in str({genome_build}) else "hg19ToHg38"
            )
            
            

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
