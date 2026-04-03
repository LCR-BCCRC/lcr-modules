"""Workflow for creating FS scores for each fragment length.

Get length distributions of mutated reads from input tumor samples.

Then get length distribution of reads from the same gene targets from
samples of healthy individuals.

Then sample the same number of reads from both distrubutions and create
a density function for that each length of fragment occuring in either group.
Then take the log2 of the ratio of the two to create fragmentation score
for that particular length of fragment. Any length with fewer than 20 reads
can be given a score of 0. The distrubution can be created for a specified
length range. This process is bootstrapped to smooth out the distibution
and reduce sampling error.

wildcards is sample_id which will be in the bams and maf files, and output and log paths.

The input functions will cooridinate making all the profiles for all the samples, as the output
from this workflow is a single file.

"""
import os
import glob
import pandas as pd

SCRIPTS = config['lcr-modules']["fragmenticity"]["scripts_dir"]
WD = os.path.join(config["lcr-modules"]["_shared"]["working_dir"], "FragmentScoreAtlas")

## input functions

def MutatedProfiles(wildcards):
    """Get the mutated profiles for all samples."""
    mutated_samples = config['lcr-modules']["fragmenticity"]["tumor_samples"]

    return expand(rules.ProfileMutatedFragments.output.mut_profile, sample_id = mutated_samples)

def HealthyProfiles(wildcards):
    """Get the healthy profiles for all samples."""
    healthy_samples = config['lcr-modules']["fragmenticity"]["healthy_samples"]

    return expand(rules.ProfileHealthyFragments.output.healthy_profile, sample_id = healthy_samples)

## rules

rule ProfileMutatedFragments:
    input:
        tum_bam = config['lcr-modules']["fragmenticity"]["bam"],
        tum_maf = config['lcr-modules']["fragmenticity"]["maf"],
    output:
        mut_profile = os.path.join(WD, "01-MutateReadProfiles", "{sample_id}_MFP.tsv"),
    params:
        script = os.path.join(SCRIPTS, "ProfileMutatedReads.py"),
    log:
        os.path.join(WD, "logs", "01-MutateReadProfiles", "{sample_id}_mutated_fragments_profile.log"),
    threads: 1
    resources:
        mem_mb = 4000
    conda:
        "analysis"
    shell:
        """
        python {params.script} \
        --bam {input.tum_bam} \
        --maf {input.tum_maf} \
        --exclude_chip_genes \
        --output {output.mut_profile} \
        > {log} 2>&1
        """

rule ProfileHealthyFragments:
    input:
        healthy_bam = config['lcr-modules']["fragmenticity"]["bam"],
    output:
        healthy_profile = os.path.join(WD, "02-HealthyReadProfiles", "{sample_id}_FP.tsv"),
    params:
        script = os.path.join(SCRIPTS, "ProfileHealthyReads.py"),
    log:
        os.path.join(WD, "logs", "02-HealthyReadProfiles", "{sample_id}_healthy_fragments_profile.log"),
    threads: 1
    resources:
        mem_mb = 4000
    conda:
        "analysis"
    shell:
        """
        python {params.script} \
        --bam {input.healthy_bam} \
        --output {output.healthy_profile} \
        > {log} 2>&1
        """

rule CaclulateFS_distribution:
    input:
        mut_profile = MutatedProfiles,
        healthy_profile = HealthyProfiles,
    output:
        fs_distribution = os.path.join(WD, "99-FSDistribution", "FS_distribution.tsv"),
    params:
        script = os.path.join(SCRIPTS, "CalculateFSDistribution.py"),
    log:
        os.path.join(WD, "logs", "99-FSDistribution", "FS_distribution.log"),
    threads: 1
    resources:
        mem_mb = 4000
    conda:
        "analysis"
    shell:
        """
        python {params.script} \
        --mutated {input.mut_profile} \
        --healthy {input.healthy_profile} \
        --output {output.fs_distribution} \
        > {log} 2>&1
        """

rule all_CreateFS:
    input:
        fs_distribution = rules.CaclulateFS_distribution.output.fs_distribution,

