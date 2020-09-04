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
    version = "1.0",
    subdirectories = ["inputs", "seghg38tobedhg38", "bedhg38tobedhg19", "bedhg19toseghg19", "outputs"])

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _liftover_input_seg:
    input:
        seg = CFG["inputs"]["sample_seg"]
    output:
        seg = CFG["dirs"]["inputs"] + "seg/{tumour_sample_id}--{normal_sample_id}.{tool}.igv.seg"
    run:
        op.relative_symlink(input.seg, output.seg)


# Convert initial seg file into bed format
rule _hg38seg_2_hg38bed:
    input:
        seg_hg38 = rules._liftover_input_seg.output
    output:
        bed_hg38 = CFG["dirs"]["seghg38tobedhg38"] + "bed/{tumour_sample_id}--{normal_sample_id}.{tool}.igv.bed",
        header = CFG["dirs"]["seghg38tobedhg38"] + "bed/{tumour_sample_id}--{normal_sample_id}.{tool}.igv.bed.header"
    log:
        stderr = CFG["logs"]["seghg38tobedhg38"] + "{tumour_sample_id}--{normal_sample_id}.{tool}.stderr.log"
    params:
        opts = CFG["options"]["seghg38tobedhg38"],
        chr_colNum = CFG["options"]["chr_colNum"],
        start_colNum = CFG["options"]["start_colNum"],
        end_colNum = CFG["options"]["end_colNum"],
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.seg_hg38} 
        --output {output.bed_hg38} 
        --chromColnum {params.chr_colNum} 
        --startColnum {params.start_colNum} 
        --endColnum {params.end_colNum}
        {params.opts} 
        2> {log.stderr}
        """)


# Convert the bed file in hg38 coordinates into hg19 coordinates
rule _hg38bed_2_hg19bed:
    input:
        seg_hg38 = rules._hg38seg_2_hg38bed.output.bed_hg38
    output:
        bed_hg19 = CFG["dirs"]["bedhg38tobedhg19"] + "bed/{tumour_sample_id}--{normal_sample_id}.{tool}.hg19.igv.bed",
        unmapped = CFG["dirs"]["bedhg38tobedhg19"] + "unmapped/{tumour_sample_id}--{normal_sample_id}.{tool}.unmapped.hg19.igv.bed"
    log:
        stderr = CFG["logs"]["bedhg38tobedhg19"] + "{tumour_sample_id}--{normal_sample_id}.{tool}.stderr.log"
    params:
        opts = CFG["options"]["bedhg38tobedhg19"],
        mismatch = CFG["options"]["min_mismatch"]       
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        liftOver -minMatch={params.mismatch}
        {input} {params.opts} 
        {output.bed_hg19} {output.unmapped}
        2> {log.stderr}
        """)

# Convert the bed file in hg19 coordinates into seg format
rule _hg19bed_2_hg19seg:
    input:
        bed_hg19 = rules._hg38bed_2_hg19bed.output.bed_hg19,
        headers = rules._hg38seg_2_hg38bed.output.header
    output:
        seg_hg19 = CFG["dirs"]["bedhg19toseghg19"] + "seg/{tumour_sample_id}--{normal_sample_id}.{tool}.hg19.igv.seg"
    log:
        stderr = CFG["logs"]["bedhg19toseghg19"] + "{tumour_sample_id}--{normal_sample_id}.{tool}.stderr.log"
    params:
        opts = CFG["options"]["bedhg19toseghg19"]      
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.bed_hg19}
        --column-header {input.headers}
        --output {output.seg_hg19} 
        {params.opts} 
        2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _liftover_output_seg:
    input:
        seg = rules._hg19bed_2_hg19seg.output.seg_hg19
    output:
        seg = CFG["dirs"]["outputs"] + "seg/{tumour_sample_id}--{normal_sample_id}.{tool}.hg19.igv.seg"
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
            #repeat the tool name N times in expand so each pair in run is used
            tool=[CFG["tool"]] * len(CFG["runs"]["tumour_sample_id"])
            )
            
            

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
