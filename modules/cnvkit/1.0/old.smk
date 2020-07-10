# ------------------------------------------------------------------------------
# Author: Helena Winata
# Date created: 29 Jan 2020
# ------------------------------------------------------------------------------
# Estimate copy number variant from exome data using CNVkit
#
# input:
# output:
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Python Packages and config file
# ------------------------------------------------------------------------------
from time import localtime, strftime

configfile: modulePATH + 'config/cnvKit_config.yaml'
# ------------------------------------------------------------------------------
# Snakemake Rule
# ------------------------------------------------------------------------------
rule cnvkit_all:
    input:
        config['cnv']['dir'] + '/merged_heatmap.pdf',
        expand('{cnvDIR}/{sample}_{plot}.pdf', cnvDIR = config['cnv']['dir'], sample = SAMPLE_PAIRS.values(), plot = ['scatter', 'diagram'])

rule cnvkit_batch:
    input:
        bamT = expand('{bamDIR}/{sampleT}.bam', bamDIR = (aligner + '/bam'), sampleT = SAMPLE_PAIRS.values()),
        bamN = expand('{bamDIR}/{sampleN}.bam', bamDIR = (aligner + '/bam'), sampleN = SAMPLE_PAIRS.keys())
    output:
        outCnr = expand('{cnvDIR}/{sampleT}.cnr', cnvDIR = config['cnv']['dir'], sampleT = SAMPLE_PAIRS.values()),
        outRef = config['cnv']['dir'] + '/normal_reference.cnn'
    params:
        access = config['ref'][genomeVer]['access'],
        targets = config['ref'][genomeVer]['targets'],
        fasta = config['ref'][genomeVer]['ref_fa']
    log: 'logs/' + config['cnv']['dir'] + '/all_cnvkit_run.' + strftime("%Y-%m-%d.%H-%M-%S", localtime()) + '.log'
    conda:
        config['cnv']['conda']
    shell:
        'cnvkit.py batch {input.bamT} --normal {input.bamN} --targets {params.targets} --fasta {params.fasta} --access {params.access} --processes 8 --output-reference {output.outRef} --scatter --diagram --output-dir {config[cnv][dir]} &> {log}'

# From baits and tumor/normal BAMs
cnvkit.py batch *Tumor.bam --normal *Normal.bam \
    --targets my_baits.bed --annotate refFlat.txt \
    --fasta hg19.fasta --access data/access-5kb-mappable.hg19.bed \
    --output-reference my_reference.cnn --output-dir results/ \
    --diagram --scatter

# Reusing a reference for additional samples
cnvkit.py batch *Tumor.bam -r Reference.cnn -d results/

# Reusing targets and antitargets to build a new reference, but no analysis
cnvkit.py batch -n *Normal.bam --output-reference new_reference.cnn \
    -t my_targets.bed -a my_antitargets.bed \
    -f hg19.fasta -g data/access-5kb-mappable.hg19.bed

rule cnvkit_segment:
    input:
        config['cnv']['dir'] + '/{sample}.cnr'
    output:
        config['cnv']['dir'] + '/{sample}.cns'
    log:
        'logs/' + config['cnv']['dir'] + '/{sample}_segment.' + strftime("%Y-%m-%d.%H-%M-%S", localtime()) + '.log'
    conda:
        config['cnv']['conda']
    shell:
        'cnvkit.py segment --output {output} {input} &> {log}'

rule cnvkit_call:
    input:
        config['cnv']['dir'] + '/{sample}.cns'
    output:
        config['cnv']['dir'] + '/{sample}.call'
    log:
        'logs/' + config['cnv']['dir'] + '/{sample}_call.' + strftime("%Y-%m-%d.%H-%M-%S", localtime()) + '.log'
    conda:
        config['cnv']['conda']
    shell:
        'cnvkit.py call --output {output} {input} &> {log}'

rule cnvkit_scatter:
    input:
        inCnr = config['cnv']['dir'] + '/{sample}.cnr',
        inCns = config['cnv']['dir'] + '/{sample}.cns'
    output:
        config['cnv']['dir'] + '/{sample}_scatter.pdf'
    params:
        chr = config['ref'][genomeVer]['chr_file']
    log:
        'logs/' + config['cnv']['dir'] + '/{sample}_scatter.' + strftime("%Y-%m-%d.%H-%M-%S", localtime()) + '.log'
    conda:
        config['cnv']['conda']
    shell:
        'cnvkit.py scatter {input.inCnr} --range-list {params.chr} --segment {input.inCns} --output {output} &> {log}'

rule cnvkit_diagram:
    input:
        inCnr = config['cnv']['dir'] + '/{sample}.cnr',
        inCns = config['cnv']['dir'] + '/{sample}.cns'
    output:
        config['cnv']['dir'] + '/{sample}_diagram.pdf'
    log:
        'logs/' + config['cnv']['dir'] + '/{sample}_diagram.' + strftime("%Y-%m-%d.%H-%M-%S", localtime()) + '.log'
    conda:
        config['cnv']['conda']
    shell:
        'cnvkit.py diagram {input.inCnr} --segment {input.inCns} --output {output} &> {log}'

rule cnvkit_heatmap:
    input:
        expand('{cnvDIR}/{sample}.call', cnvDIR = config['cnv']['dir'], sample = SAMPLE_PAIRS.values())
    output:
        config['cnv']['dir'] + '/merged_heatmap.pdf'
    log:
        'logs/' + config['cnv']['dir'] + '/merged_heatmap.' + strftime("%Y-%m-%d.%H-%M-%S", localtime()) + '.log'
    conda:
        config['cnv']['conda']
    shell:
        'cnvkit.py heatmap --output {output} {input} &> {log}'
