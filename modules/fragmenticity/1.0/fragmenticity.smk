import pandas as pd
import os
import sys

FRAGSAMPLES = config["lcr-modules"]["_shared"]["samples"]
WD = config["lcr-modules"]["_shared"]["working_dir"]
UTILSDIR = os.path.join(config["lcr-modules"]["_shared"]["lcr-modules"], "/modules/fragmenticity/1.0/scripts")
LOGDIR = os.path.join(WD, "logs")

rule Generate_FS_Regional:
    input:
        bam = config["lcr-modules"]["fragmenticity"]["bam"],
    output:
        fs = os.path.join(WD, "99-fragmenticity", "{sample_id}_FragmentScore.tsv"),
        regional = os.path.join(WD, "99-fragmenticity", "{sample_id}_RegionalFragmentScores.tsv"),
        histogram = os.path.join(WD, "99-fragmenticity", "{sample_id}_regional_fragment_distribution.png"),
    params:
        sample_id = lambda wildcards: wildcards.sample_id,
        script = os.path.join(UTILSDIR, "ScoreRegionalFrags.py"),
        frag_ref = config["lcr-modules"]["fragmenticity"]["fragment_score_reference"],
        read_count = config["lcr-modules"]["fragmenticity"]["read_count"],
        min_length = config["lcr-modules"]["fragmenticity"]["min_length"],
        max_length = config["lcr-modules"]["fragmenticity"]["max_length"],
        top_k_regions = config["lcr-modules"]["fragmenticity"]["top_k_regions"],
        min_reads_per_region = config["lcr-modules"]["fragmenticity"]["min_reads_per_region"],
        bed = config["lcr-modules"]["fragmenticity"]["bed_file"],
        seed = config["lcr-modules"]["fragmenticity"]["seed"],
        output_dir = os.path.join(WD, "99-fragmenticity"),
    log:
        os.path.join(LOGDIR, "fragmenticity", "{sample_id}_FragmentScore_Regional.log"),
    conda:
        "analysis",
    threads: 1
    resources:
        mem_mb=8000,  # Increased memory for regional analysis
        time="02:00:00",  # Increased time for regional processing
    shell:
        """
        python {params.script} \
    --bam_file {input.bam} \
    --bed_file {params.bed} \
    --fragment_score_reference {params.frag_ref} \
    --sample_name {params.sample_id} \
    --output {params.output_dir} \
    --read_count {params.read_count} \
    --min_length {params.min_length} \
    --max_length {params.max_length} \
    --top_k_regions {params.top_k_regions} \
    --min_reads_per_region {params.min_reads_per_region} \
    --seed {params.seed} \
    > {log} 2>&1
        """

rule Generate_FS_all:
    input:
        expand(rules.Generate_FS_Regional.output.fs, 
               sample_id=FRAGSAMPLES["sample_id"])