#!/usr/bin/env snakemake


# ---------------------------------------------------------------------------- #
##### ATTRIBUTION #####
# ---------------------------------------------------------------------------- #

# Original snakemake author: Jasper Wong
# Module author: Jasper Wong
# Additional contributors: N/A


# ---------------------------------------------------------------------------- #
##### SETUP #####
# ---------------------------------------------------------------------------- #

### Modules ###

import pandas as pd
import numpy as np
import oncopipe as op
import glob
import os

### Directories ###

resultsDir = "results/"
logDir = "log/"


# Setup module and store module-specific configuration in `CFG`.
CFG = op.setup_module(
    name = "ichorcna", 
    version = "1.0",
    subdirectories = ["inputs", "readDepth", "seg", "outputs"]
)

localrules:
    _ichorcna_input_bam,
    _ichorcna_output,
    _ichorcna_all

# ---------------------------------------------------------------------------- #
##### RULES #####
# ---------------------------------------------------------------------------- #


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _ichorcna_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai" # specific to readCounter
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
     

rule _ichorcna_read_counter:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    output:
        CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/{sample_id}.bin{binSize}.wig"
    params:
        readCounter = CFG["options"]["readcounter"]["readCounterScript"],
        binSize = CFG["options"]["readcounter"]["binSize"],
        qual = CFG["options"]["readcounter"]["qual"],
        chrs = op.switch_on_wildcard("genome_build", CFG["options"]["readcounter"]["chrs"])
    conda: CFG["conda_envs"]["ichorcna"]
    threads: CFG["threads"]["readcounter"]
    resources:
        **CFG["resources"]["readcounter"]
    log:
        CFG["logs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/{sample_id}.bin{binSize}.log"
    shell:
        "{params.readCounter} {input.bam} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"


rule _run_ichorcna:
    input:
        tum = CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}.bin{binSize}.wig",
        # norm = CFG["dirs"]["readDepth"] + "{seq_type}--{genome_build}/{binSize}/{normal_id}.bin{binSize}.wig"
    output:
        corrDepth = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.correctedDepth.txt",
        param = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.params.txt",
        cna = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cna.seg",
        segTxt = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.seg.txt",
        seg = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.seg",
        plot = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}/{tumour_id}_genomeWide.pdf",
        #rdata = "results/ichorCNA/{sample_id}/{sample_id}.RData"
    params:
        outDir = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}/",
        rscript = CFG["options"]["run"]["ichorCNA_rscript"],
        name = "{tumour_id}",
        ploidy = CFG["options"]["run"]["ichorCNA_ploidy"],
        normal = CFG["options"]["run"]["ichorCNA_normal"],
        gcwig = op.switch_on_wildcard("binSize", CFG["options"]["run"]["ichorCNA_gcWig"]),
        mapwig = op.switch_on_wildcard("binSize", CFG["options"]["run"]["ichorCNA_mapWig"]),
        normalpanel = op.switch_on_wildcard("binSize", CFG["options"]["run"]["ichorCNA_normalPanel"]),
        estimateNormal = CFG["options"]["run"]["ichorCNA_estimateNormal"],
        estimatePloidy = CFG["options"]["run"]["ichorCNA_estimatePloidy"],
        estimateClonality = CFG["options"]["run"]["ichorCNA_estimateClonality"],
        scStates = CFG["options"]["run"]["ichorCNA_scStates"],
        maxCN = CFG["options"]["run"]["ichorCNA_maxCN"],
        includeHOMD = CFG["options"]["run"]["ichorCNA_includeHOMD"],
        chrs = op.switch_on_wildcard("genome_build", CFG["options"]["run"]["ichorCNA_chrs"]),
        chrTrain = op.switch_on_wildcard("genome_build", CFG["options"]["run"]["ichorCNA_chrTrain"]),
        genomeBuild = "{genome_build}",
        genomeStyle = op.switch_on_wildcard("genome_build", CFG["options"]["run"]["ichorCNA_genomeStyle"]),
        centromere = op.switch_on_wildcard("genome_build", CFG["options"]["run"]["ichorCNA_centromere"]),
        fracReadsChrYMale = CFG["options"]["run"]["ichorCNA_fracReadsInChrYForMale"],
        minMapScore = CFG["options"]["run"]["ichorCNA_minMapScore"],
        maxFracGenomeSubclone = CFG["options"]["run"]["ichorCNA_maxFracGenomeSubclone"],
        maxFracCNASubclone = CFG["options"]["run"]["ichorCNA_maxFracCNASubclone"],
        exons = CFG["options"]["run"]["ichorCNA_exons"],
        txnE = CFG["options"]["run"]["ichorCNA_txnE"],
        txnStrength = CFG["options"]["run"]["ichorCNA_txnStrength"],
        plotFileType = CFG["options"]["run"]["ichorCNA_plotFileType"],
        plotYlim = CFG["options"]["run"]["ichorCNA_plotYlim"],
        libdir = CFG["options"]["run"]["ichorCNA_libdir"]
    conda: CFG["conda_envs"]["ichorcna"]
    threads: CFG["threads"]["run"]
    resources:
        **CFG["resources"]["run"]
    log:
        stdout = CFG["logs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}.stdout.log",
        stderr = CFG["logs"]["seg"] + "{seq_type}--{genome_build}/{binSize}/{tumour_id}--{normal_id}--{pair_status}.stderr.log"
    shell:
        "Rscript {params.rscript} --id {params.name} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log.stdout} 2> {log.stderr}"

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _ichorcna_output:
    input:
        corrDepth = str(rules._run_ichorcna.output.corrDepth),
        param = str(rules._run_ichorcna.output.param),
        cna = str(rules._run_ichorcna.output.cna),
        segTxt = str(rules._run_ichorcna.output.segTxt),
        seg = str(rules._run_ichorcna.output.seg),
        plot = str(rules._run_ichorcna.output.plot)
    output:
        corrDepth = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/corrDepth/{binSize}/{tumour_id}--{normal_id}--{pair_status}.corrDepth.txt",
        param = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/param/{binSize}/{tumour_id}--{normal_id}--{pair_status}.param.txt",
        cna = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/binCNA/{binSize}/{tumour_id}--{normal_id}--{pair_status}.cna.seg",
        segTxt = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/seg_txt/{binSize}/{tumour_id}--{normal_id}--{pair_status}.seg.txt",
        seg = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/seg/{binSize}/{tumour_id}--{normal_id}--{pair_status}.seg",
        plot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/plot/{binSize}/{tumour_id}--{normal_id}--{pair_status}_genomeWide.pdf"
    run:
        op.relative_symlink(input.corrDepth, output.corrDepth)
        op.relative_symlink(input.param, output.param)
        op.relative_symlink(input.cna, output.cna)
        op.relative_symlink(input.segTxt, output.segTxt)
        op.relative_symlink(input.seg, output.seg)
        op.relative_symlink(input.plot, output.plot)

# Generates the target sentinels for each run, which generate the symlinks
rule _ichorcna_all:
    input:
        expand(
            [
                str(rules._ichorcna_output.output.corrDepth),
                str(rules._ichorcna_output.output.param),
                str(rules._ichorcna_output.output.cna),
                str(rules._ichorcna_output.output.segTxt),
                str(rules._ichorcna_output.output.seg),
                str(rules._ichorcna_output.output.plot)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            pair_status=CFG["runs"]["pair_status"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            binSize=str(CFG["options"]["readcounter"]["binSize"]))



##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)