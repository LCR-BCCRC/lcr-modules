import os

PANEL_NAME = config["lcr-modules"]["artifact_alert"]["panel_name"]
OUTDIR = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "artifact_alert", "1.0", PANEL_NAME)
LOGDIR = os.path.join(OUTDIR, "logs")
MUT_SAMPLES = config["lcr-modules"]["_shared"]["samples"]

# rules 
rule generate_pileup:
    input:
        bam= config["lcr-modules"]["artifact_alert"]["input_bam"], 
        ref= config["lcr-modules"]["_shared"]["ref_genome"],
        bed= config["lcr-modules"]["artifact_alert"]["target_bed"]
    output:
        pileup= os.path.join(OUTDIR, "01-pileup", "{sample_id}.pileup")
    params:
        min_mapq= config["lcr-modules"]["artifact_alert"]["min_mapq"],
        min_baseq= config["lcr-modules"]["artifact_alert"]["min_baseq"]
    conda:
        "envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb=3000
    log:
        os.path.join(LOGDIR, "01-pileup", "generate_pileup_{sample_id}.log")
    shell:
        """
        samtools mpileup \
            -f {input.ref} \
            -l {input.bed} \
            -Q {params.min_baseq} \
            -q {params.min_mapq} \
            --no-BAQ \
            {input.bam} > {output.pileup} 2> {log}
        """

def calc_mut_rate_mem(wildcards, attempt, input):
    if attempt == 1:
        return input.size_mb + 1000
    elif attempt == 2:
        return input.size_mb + 4000
    else:
        return input.size_mb + 8000

rule calculate_mutation_rates:
    input:
        pileup=rules.generate_pileup.output.pileup
    output:
        mutation_rates= os.path.join(OUTDIR, "02-mutation_rates", "{sample_id}_mutation_rates.tsv")
    params:
        min_depth=config["lcr-modules"]["artifact_alert"]["min_depth"],
    conda:
        "envs/python_scripts.yaml"
    threads: 1
    resources:
        mem_mb= calc_mut_rate_mem
    log:
        os.path.join(LOGDIR, "02-mutation_rates", "calculate_mutation_rates_{sample_id}.log")
    script:
        "scripts/calculate_mutation_rates.py"

def get_mutation_rate_files(wildcards):
    """Input function to collect all mutation rate files from previous step"""
    return expand(rules.calculate_mutation_rates.output.mutation_rates,
                  sample_id = MUT_SAMPLES["sample_id"])

rule aggregate_mutation_rates:
    input:
        mutation_files = get_mutation_rate_files
    output:
        aggregated= os.path.join(OUTDIR, "03-aggregated", "background_mutation_rates.tsv")
    conda:
        "daily"
    threads: 1
    resources:
        mem_mb=5000
    log:
        os.path.join(LOGDIR, "03-aggregated", "aggregate_mutation_rates.log"),
    script:
        "scripts/aggregate_mutation_rates.py"

rule FetchMutationRate:
    input:
        rules.aggregate_mutation_rates.output.aggregated
