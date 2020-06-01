from cellranger import *

##### cellranger.py ###
import pandas as pd
import glob

# make sample sheet for each chip id
def cellranger_setup(samples, samplesheetDIR):
    '''
    generating samplesheet and sample dictionary for cellranger
    '''
    print(samples)
    chipID = list(samples.chip_id.unique())
    sample_dict = {}
    for i in chipID:
        sub_samples = samples[samples['chip_id']==i]
        samplesheet = sub_samples[['lane', 'sample_id', 'index']]
        samplesheet.columns = ['lane', 'sample', 'index']
	samplesheet.to_csv(f'{samplesheetDIR}/{i}_samplesheet.csv', sep = ',', index = False)
#        sample_dict[i] = list(sub_samples['sample_id'])
#        return(sample_dict)

def get_run_files(pattern, dir, check):
    fileDIR = dir + '/*' + pattern + '*'
    run = glob.glob(fileDIR)[0]
    file = [run]
    for suffix in check:
        file.append(glob.glob(run + '/' + suffix)[0])
    return(file)

def get_10x_raw(wildcards):
    file = get_run_files(wildcards.runID, dir = rawDIR, check = ['RTAComplete*', 'RunInfo*', 'RunParameters*'])
    # = ['RTAComplete*', 'RunInfo*', 'RunParameters*'])
    return(file)
################################

configfile: modulePATH + 'config/cellranger_config.yaml'
# ------------------------------------------------------------------------------
# Snakemake Rule
# ------------------------------------------------------------------------------
rule cellranger_all:
    input:
        expand('stamps/{runID}_{sampleID}_{type}.stamp', zip, runID = list(SAMPLES['chip_id']), sampleID = list(SAMPLES['sample_id']), type = list(SAMPLES['analysis']))

rule cellranger_mkfastq:
    input:
        get_10x_raw
    output:
        'stamps/{runID}_mkfastq.stamp'
    params:
        ss = ssDIR + '/{runID}_samplesheet.csv',
        baseArgs = config['cellranger']['args']['mkfastq']
    shell:
        '''
	    cellranger mkfastq {params.baseArgs} --id={wildcards.runID} --run={input[0]} --csv={params.ss}
	    touch {output}
	    '''

rule cellranger_count:
    input:
        'stamps/{runID}_mkfastq.stamp'
    output:
        'stamps/{runID}_{sampleID}_count.stamp'
    params:
         fastq = '{runID}/outs/fastq_path',
         baseArgs = config['cellranger']['args']['count'],
         ref = config['cellranger']['ref']['trx']
    shell:
        '''
	    cellranger count {params.baseArgs} --id={wildcards.sampleID} --sample={wildcards.sampleID} --fastqs={params.fastq} --transcriptome={params.ref}
	    touch {output}
	    '''

rule cellranger_vdj:
    input:
        'stamps/{runID}_mkfastq.stamp'
    output:
        'stamps/{runID}_{sampleID}_vdj.stamp'
    params:
        fastq = '{runID}/outs/fastq_path',
        baseArgs = config['cellranger']['args']['vdj'],
        ref = config['cellranger']['ref']['vdj']
    shell:
        '''
        cellranger vdj {params.baseArgs} --id={wildcards.sampleID} --fastqs={params.fastq} --reference={params.ref}
        touch {output}
        '''
