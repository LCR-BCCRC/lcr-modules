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
    subdirectories = ["inputs", "readDepth", "outputs"]
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
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
     

rule targets_genomic_cov_ichorCNA:
    input:
        expand(resultsDir + "ichorcna/{sample_id}/{sample_id}.cna.seg", zip,
        sample_id = config["sample_id"]),
        expand(resultsDir + "readDepth/{sample_id}.bin{binSize}.wig", zip,
        sample_id=config["sample_id"], binSize=str(config["binSize"]))
        


rule ichorcna_read_counter:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai"
    output:
        CFG["dirs"]["readDepth"] + "readDepth/{sample_id}.bin{binSize}.wig"		
    params:
        readCounter = CFG["options"]["readcounter"]["readCounterScript"],
        binSize = CFG["options"]["readcounter"]["binSize"],
        qual = CFG["options"]["readcounter"]["qual"],
        chrs = CFG["options"]["readcounter"]["chrs"]
    conda: CFG["conda_envs"]["ichorcna"]
    threads: CFG["threads"]["readcounter"]
    resources:
        mem = CFG["mem_mb"]["readcounter"]
    log:
        CFG["logs"]["readcounter"] + "{sample_id}.bin{binSize}.0.05.log"
    shell:
        "{params.readCounter} {input.bam} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"


rule _run_ichorcna:
    input:
        tum=resultsDir + "readDepth/{sample_id}.bin" + str(config["binSize"]) + ".wig",
        #norm=lambda wildcards: "results/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
    output:
        #corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
        #param="results/ichorCNA/{tumor}/{tumor}.params.txt",
        cna=resultsDir + "ichorCNA/{sample_id}/{sample_id}.cna.seg",
        #segTxt="results/ichorCNA/{sample_id}/{sample_id}.seg.txt",
        #seg="results/ichorCNA/{sample_id}/{sample_id}.seg",
        #rdata="results/ichorCNA/{sample_id}/{sample_id}.RData"
    params:
        # runRscript="/projects/dscott_prj/CCSRI_1500/lowpassWGS/envs/ichorcna/bin/Rscript",
        outDir=resultsDir + "ichorCNA/{sample_id}/",
        rscript=CFG["options"]["run"]["ichorCNA_rscript"],
        id="{sample_id}",
        ploidy=CFG["options"]["run"]["ichorCNA_ploidy"],
        normal=CFG["options"]["run"]["ichorCNA_normal"],
        gcwig=CFG["options"]["run"]["ichorCNA_gcWig"],
        mapwig=CFG["options"]["run"]["ichorCNA_mapWig"],
        normalpanel=CFG["options"]["run"]["ichorCNA_normalPanel"],
        estimateNormal=CFG["options"]["run"]["ichorCNA_estimateNormal"],
        estimatePloidy=CFG["options"]["run"]["ichorCNA_estimatePloidy"],
        estimateClonality=CFG["options"]["run"]["ichorCNA_estimateClonality"],
        scStates=CFG["options"]["run"]["ichorCNA_scStates"],
        maxCN=CFG["options"]["run"]["ichorCNA_maxCN"],
        includeHOMD=CFG["options"]["run"]["ichorCNA_includeHOMD"],
        chrs=CFG["options"]["run"]["ichorCNA_chrs"],
        chrTrain=CFG["options"]["run"]["ichorCNA_chrTrain"],
        genomeBuild=CFG["options"]["run"]["ichorCNA_genomeBuild"],
        genomeStyle=CFG["options"]["run"]["ichorCNA_genomeStyle"],
        centromere=CFG["options"]["run"]["ichorCNA_centromere"],
        fracReadsChrYMale=CFG["options"]["run"]["ichorCNA_fracReadsInChrYForMale"],
        minMapScore=CFG["options"]["run"]["ichorCNA_minMapScore"],
        maxFracGenomeSubclone=CFG["options"]["run"]["ichorCNA_maxFracGenomeSubclone"],
        maxFracCNASubclone=CFG["options"]["run"]["ichorCNA_maxFracCNASubclone"],
        exons=CFG["options"]["run"]["ichorCNA_exons"],
        txnE=CFG["options"]["run"]["ichorCNA_txnE"],
        txnStrength=CFG["options"]["run"]["ichorCNA_txnStrength"],
        plotFileType=CFG["options"]["run"]["ichorCNA_plotFileType"],
        plotYlim=CFG["options"]["run"]["ichorCNA_plotYlim"],
        libdir=CFG["options"]["run"]["ichorCNA_libdir"]
    conda: CFG["conda_envs"]["ichorcna"]
    threads: CFG["threads"]["run"]
    resources:
        mem = CFG["mem_mb"]["run"]
    log:
        logDir + "ichorCNA/{sample_id}.log"	
    shell:
        "Rscript {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"
