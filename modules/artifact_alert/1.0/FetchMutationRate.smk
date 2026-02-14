import os

################################ setup functions ################################
def get_chromosomes(genome_version):
    if genome_version in ["hg19", "GRCh37"]:
        
        return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    elif genome_version in ["hg38", "GRCh38"]:
        return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    else:
        raise ValueError(f"Unsupported genome version: {genome_version}")

################################ global variables ################################
PANEL_NAME = config["lcr-modules"]["artifact_alert"]["panel_name"]
OUTDIR = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "artifact_alert", "1.0", PANEL_NAME)
LOGDIR = os.path.join(OUTDIR, "logs")
MUT_SAMPLES = config["lcr-modules"]["_shared"]["samples"]
SAMPLE_TRACKING_FILE = os.path.join(OUTDIR, f"{PANEL_NAME}_MutRateIndex_sampletracker.txt")
CHROMOSOMES = get_chromosomes(config["lcr-modules"]["_shared"]["ref_genome_ver"])
SCRIPTS_DIR = os.path.join(config["lcr-modules"]["_shared"]["lcr-modules"], "modules", "artifact_alert", "1.0","scripts")
BED_FILE =config["lcr-modules"]["artifact_alert"].get("target_bed", "")

################################ reset final outputs ################################

# Check if user wants to reset the aggregated index
RESET_INDEX = config["lcr-modules"]["artifact_alert"].get("reset_mutation_index", False)

# If resetting, remove existing aggregated files before workflow starts
if RESET_INDEX:
    aggregated_file = os.path.join(OUTDIR, "03-aggregated", "background_mutation_rates.tsv.gz")
    index_file = f"{aggregated_file}.tbi"
    
    if os.path.exists(aggregated_file):
        os.remove(aggregated_file)
        print(f"Reset: Removed {aggregated_file}")
    
    if os.path.exists(index_file):
        os.remove(index_file)
        print(f"Reset: Removed {index_file}")

    # Also clear the sample tracker
    if os.path.exists(SAMPLE_TRACKING_FILE):
        os.remove(SAMPLE_TRACKING_FILE)
        print(f"Reset: Removed {SAMPLE_TRACKING_FILE}")

################### Filter for only new samples ################################

# Sample filtering logic
def get_new_samples():
    """Filter samples: deduplicate and remove already-processed samples"""
    # Load processed samples as a set
    processed = set()
    if os.path.exists(SAMPLE_TRACKING_FILE):
        with open(SAMPLE_TRACKING_FILE, 'r') as f:
            processed = {line.strip() for line in f if line.strip()}
    
    # Get new samples: input samples minus already processed
    all_samples = set(config["lcr-modules"]["_shared"]["samples"]["sample_id"])
    return list(all_samples - processed)

NEW_SAMPLES = get_new_samples()
# if no new samples exit snakemake workflow
if len(NEW_SAMPLES) == 0:
    # print wiht purple text

    print('\033[0;31m' + "No new samples to process. Exiting workflow.")
    exit(0)

######################## rules ################################
rule generate_pileup:
    input:
        bam= config["lcr-modules"]["artifact_alert"]["input_bam"], 
    output:
        pileup= temp(os.path.join(OUTDIR, "01-pileup", "{sample_id}_{chrom}.pileup"))
    params:
        ref= config["lcr-modules"]["_shared"]["ref_genome"],
        bed= (f"-l {BED_FILE}" if BED_FILE else ""),
        min_mapq= config["lcr-modules"]["artifact_alert"]["min_mapq"],
        min_baseq= config["lcr-modules"]["artifact_alert"]["min_baseq"]
    conda:
        "envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb=3000
    log:
        os.path.join(LOGDIR, "01-pileup", "generate_pileup_{sample_id}_{chrom}.log")
    shell:
        """
        samtools mpileup \
            -f {params.ref} {params.bed} \
            -Q {params.min_baseq} \
            -q {params.min_mapq} \
            -r {wildcards.chrom} \
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
        mutation_rates= temp(os.path.join(OUTDIR, "02-mutation_rates", "{sample_id}_{chrom}_mutation_rates.tsv.gz"))
    params:
        min_depth=config["lcr-modules"]["artifact_alert"]["min_depth"],
    conda:
        "envs/samtools.yaml"
    threads: 1
    resources:
        mem_mb= calc_mut_rate_mem
    log:
        os.path.join(LOGDIR, "02-mutation_rates", "calculate_mutation_rates_{sample_id}_{chrom}.log")
    script:
        "scripts/calculate_mutation_rates.py"

def get_mutation_rate_files(wildcards):
    """Input function to collect mutation rate files from NEW samples only"""
    if len(NEW_SAMPLES) == 0:
        return []
    return expand(rules.calculate_mutation_rates.output.mutation_rates,
                  sample_id=NEW_SAMPLES, chrom=CHROMOSOMES)

rule aggregate_mutation_rates:
    input:
        mutation_files = get_mutation_rate_files
    output:
        sentinel = touch(os.path.join(OUTDIR, "03-aggregated", ".aggregate_complete"))
    params:
        sample_ids = NEW_SAMPLES,
        chromosomes = CHROMOSOMES,
        sample_tracker = SAMPLE_TRACKING_FILE,
        script = os.path.join(SCRIPTS_DIR, "aggregate_mutation_rates.py"),
        aggregated = os.path.join(OUTDIR, "03-aggregated", "background_mutation_rates.tsv.gz"),
        index = os.path.join(OUTDIR, "03-aggregated", "background_mutation_rates.tsv.gz.tbi")
    conda:
        "envs/samtools.yaml"
    threads: config["lcr-modules"]["artifact_alert"]["agg_threads"]
    resources:
        mem_mb = lambda wildcards, threads: threads * int(config["lcr-modules"]["artifact_alert"]["agg_mem_per_thread"])
    log:
        os.path.join(LOGDIR, "03-aggregated", "aggregate_mutation_rates.log")
    shell:
        """
        python {params.script} \
            -i {input.mutation_files} \
            -o {params.aggregated} \
            -s {params.sample_ids} \
            -c {params.chromosomes} \
            -t {params.sample_tracker} \
            --threads {threads} \
            --log {log} >> {log} 2>&1
        """

rule all_artifact_alert:
    input:
        rules.aggregate_mutation_rates.output.sentinel
