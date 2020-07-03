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

# Define rules to be run locally when using a compute cluster
localrules:
    _liftover_input_seg,
    _hg38seg_2_hg38bed,
    _hg38bed_2_hg19bed,
    _hg19bed_2_hg19seg,
    _liftover_output_seg,
    _liftover_all,


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _liftover_input_seg:
    input:
        seg = CFG["inputs"]["sample_seg"]
    output:
        seg = CFG["dirs"]["inputs"] + "seg/{sample_id}_subclones.hg38.igv.seg"
    run:
        op.relative_symlink(input.seg, output.seg)


# Convert initial seg file into bed format
rule _hg38seg_2_hg38bed:
    input:
        seg_hg38 = CFG["dirs"]["inputs"] + "seg/{sample_id}_subclones.hg38.igv.seg"
    output:
        bed_hg38 = CFG["dirs"]["seghg38tobedhg38"] + "bed/{sample_id}_subclones.hg38.igv.bed",
        header = CFG["dirs"]["seghg38tobedhg38"] + "bed/{sample_id}_subclones.hg38.igv.bed.header"
    log:
        stderr = CFG["logs"]["seghg38tobedhg38"] + "{sample_id}_seghg38tobedhg38.stderr.log"
    params:
        opts = CFG["options"]["seghg38tobedhg38"],
        chr_colNum = 2,
        start_colNum = 3,
        end_colNum = 4,
    conda:
        CFG["conda_envs"]["liftover-366"]
    resources:
        mem_mb = CFG["mem_mb"]["seghg38tobedhg38"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.seg_hg38} 
        --output {output.bed_hg38} 
        --chromColnum {params.chr_colNum} 
        --startColnum {params.start_colNum} 
        --endColnum {params.end_colNum}
        {params.opts} 
        --threads {threads}
        2> {log.stderr}
        """)


# Convert the bed file in hg38 coordinates into hg19 coordinates
rule _hg38bed_2_hg19bed:
    input:
        seg_hg38 = rules._hg38seg_2_hg38bed.output.bed_hg38
    output:
        bed_hg19 = CFG["dirs"]["bedhg38tobedhg19"] + "bed/{sample_id}_subclones.hg19.igv.bed",
        unmapped = CFG["dirs"]["bedhg38tobedhg19"] + "unmapped/{sample_id}_subclones.unmapped.hg19.igv.bed"
    log:
        stderr = CFG["logs"]["bedhg38tobedhg19"] + "{sample_id}_bedhg38tobedhg19.stderr.log"
    params:
        opts = CFG["options"]["bedhg38tobedhg19"]
        #chain = "reference/hg38ToHg19.over.chain"
    conda:
        CFG["conda_envs"]["liftover-366"]
    resources:
        mem_mb = CFG["mem_mb"]["bedhg38tobedhg19"]
    shell:
        op.as_one_line("""
        liftOver {input} {params.opts} 
        {output.bed_hg19} {output.unmapped}
        2> {log.stderr}
        """)

# Convert the bed file in hg19 coordinates into seg format
rule _hg19bed_2_hg19seg:
    input:
        bed_hg19 = rules._hg38bed_2_hg19bed.output.bed_hg19,
        headers = rules._hg38seg_2_hg38bed.output.header
    output:
        seg_hg19 = CFG["dirs"]["bedhg19toseghg19"] + "seg/{sample_id}_subclones.hg19.igv.seg"
    log:
        stderr = CFG["logs"]["bedhg19toseghg19"] + "{sample_id}_bedhg19toseghg19.stderr.log"
    params:
        opts = CFG["options"]["bedhg19toseghg19"]      
    conda:
        CFG["conda_envs"]["liftover-366"]
    resources:
        mem_mb = CFG["mem_mb"]["seghg38tobedhg38"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.bed_hg19}
        --column-header {input.headers}
        --output {output.seg_hg19} 
        {params.opts} 
        --threads {threads}
        2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _liftover_output_seg:
    input:
        seg = rules._hg19bed_2_hg19seg.output.seg_hg19
    output:
        seg = CFG["dirs"]["outputs"] + "seg/{sample_id}_subclones.hg19.igv.seg"
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
            sample_id=CFG["runs"]["tumour_sample_id"])
            
            

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
