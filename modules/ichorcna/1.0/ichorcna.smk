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
    subdirectories = ["inputs", "outputs"]
)

# ---------------------------------------------------------------------------- #
##### RULES #####
# ---------------------------------------------------------------------------- #


rule targets_genomic_cov_ichorCNA:
	input:
		expand(resultsDir + "ichorCNA/{group}--{frac05}/{group}--{frac05}.cna.seg", zip,
		group = config["group"], frac05= bamstats['frac0.05']),
		expand(resultsDir + "readDepth/{group}--{frac05}.bin{binSize}.wig", zip,
		group=config["group"], binSize=str(config["binSize"]), frac05=bamstats['frac0.05'])
		
rule genomicCov:
	input: 
		lambda wildcards: config["group"][wildcards.group]
	output:
		bam = SCRATCHdir + "{group}--{frac05}.bam"
	params:
		samtools="/projects/dscott_prj/CCSRI_1500/lowpassWGS/envs/freec/bin/samtools",
		coverage= "{frac05}"
	conda:
		CFG["conda_envs"]["samtools"]
	log:
		logDir + "ichorCNA/{group}--{frac05}.log"
	shell:
		"{params.samtools} view -bs {params.coverage} {input} > {output} 2> {log} "

rule genomicCov_index:
	input:
		bam = SCRATCHdir + "{group}--{frac05}.bam"
	output:
		bai = SCRATCHdir + "{group}--{frac05}.bam.bai"
	params:
		samtools="/projects/dscott_prj/CCSRI_1500/lowpassWGS/envs/freec/bin/samtools"
	shell:
		"{params.samtools} index {input.bam} "

rule genomicCov_read_counter:
	input:
		bam = SCRATCHdir + "{group}--{frac05}.bam",
		bai = SCRATCHdir + "{group}--{frac05}.bam.bai"
	output:
		resultsDir + "readDepth/{group}--{frac05}.bin{binSize}.wig"		
	params:
		readCounter=config["readCounterScript"],
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"]
	resources:
		mem=4
	log:
		logDir + "ichorCNA/readDepth/{group}--{frac05}.bin{binSize}.0.05.log"
	shell:
		"{params.readCounter} {input.bam} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule genomicCov_ichorCNA:
	input:
		tum=resultsDir + "readDepth/{group}--{frac05}.bin" + str(config["binSize"]) + ".wig",
		#norm=lambda wildcards: "results/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
	output:
		#corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
		#param="results/ichorCNA/{tumor}/{tumor}.params.txt",
		cna=resultsDir + "ichorCNA/{group}--{frac05}/{group}--{frac05}.cna.seg",
		#segTxt="results/ichorCNA/{group}/{group}.seg.txt",
		#seg="results/ichorCNA/{group}/{group}.seg",
		#rdata="results/ichorCNA/{group}/{group}.RData"
	params:
		runRscript="/projects/dscott_prj/CCSRI_1500/lowpassWGS/envs/ichorcna/bin/Rscript",
		outDir=resultsDir + "ichorCNA/{group}--{frac05}/",
		rscript=config["ichorCNA_rscript"],
		id="{group}--{frac05}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
		gcwig=config["ichorCNA_gcWig"],
		mapwig=config["ichorCNA_mapWig"],
		normalpanel=config["ichorCNA_normalPanel"],
		estimateNormal=config["ichorCNA_estimateNormal"],
		estimatePloidy=config["ichorCNA_estimatePloidy"],
		estimateClonality=config["ichorCNA_estimateClonality"],
		scStates=config["ichorCNA_scStates"],
		maxCN=config["ichorCNA_maxCN"],
		includeHOMD=config["ichorCNA_includeHOMD"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		genomeBuild=config["ichorCNA_genomeBuild"],
		genomeStyle=config["ichorCNA_genomeStyle"],
		centromere=config["ichorCNA_centromere"],
		fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
		minMapScore=config["ichorCNA_minMapScore"],
		maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
		maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
		exons=config["ichorCNA_exons"],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"],
		libdir=config["ichorCNA_libdir"]
	resources:
		mem=4
	log:
		logDir + "ichorCNA/{group}--{frac05}.log"	
	shell:
		"{params.runRscript} {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"
