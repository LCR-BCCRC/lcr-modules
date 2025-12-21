#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A

##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import glob
import pandas as pd
import os

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    print(f"ERROR: oncopipe version installed: {current_version}")
    print(f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment")
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section    

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["battenberg"]`
CFG = op.setup_module(
    name = "battenberg",
    version = "1.3",
    subdirectories = ["inputs", "infer_sex", "battenberg", "convert_coordinates", "fill_regions", "normalize", "plink", "outputs"],
)

#set variable for prepending to PATH based on config
SCRIPT_PATH = CFG['inputs']['src_dir']
#this is used in place of the shell.prefix() because that was not working consistently. This is not ideal. 

#this preserves the variable when using lambda functions
_battenberg_CFG = CFG

# Provide backward-compatible defaults for fit-specific settings when older
# project configs haven't defined them yet.
if "battenberg_fit" not in CFG["threads"]:
    CFG["threads"]["battenberg_fit"] = CFG["threads"]["battenberg"]
if "battenberg_fit" not in CFG["resources"]:
    CFG["resources"]["battenberg_fit"] = CFG["resources"]["battenberg"]

# Basic sanity checks for the `runs` configuration to avoid empty-expands
# which can make the `_battenberg_all` rule have no dependencies and therefore
# run without scheduling the real rules. Fail early with a helpful message.
_required_run_keys = [
    "tumour_sample_id",
    "normal_sample_id",
    "tumour_seq_type",
    "tumour_genome_build",
    "pair_status",
]

# Attempt to find the `runs` block either in the module's CFG (preferred) or
# in the top-level `config` merged by Snakemake. Older workflows sometimes
# provide `runs` only via the project config, so fail-fast raising here caused
# surprising errors. Fall back to a warning if runs are not provided, and only
# perform strict validation when a runs dict is available.
_module_runs = None
if "runs" in CFG and isinstance(CFG["runs"], dict):
    _module_runs = CFG["runs"]
else:
    _module_runs = config.get("lcr-modules", {}).get("battenberg", {}).get("runs")

# Determine whether we actually have runs provided and validate according to type.
_has_runs = False
if _module_runs is None:
    _has_runs = False
elif isinstance(_module_runs, dict):
    _has_runs = bool(_module_runs)
elif hasattr(_module_runs, 'empty'):
    # likely a pandas DataFrame
    try:
        _has_runs = (not bool(getattr(_module_runs, 'empty')))
    except Exception:
        _has_runs = True
else:
    # Unknown but truthy object
    _has_runs = True

if not _has_runs:
    # Keep behaviour permissive: warn and set an empty runs dict so other
    # parts of the module that expect CFG["runs"] won't crash immediately.
    print("WARNING: 'runs' missing under lcr-modules.battenberg in config; some module targets may have no inputs.")
    CFG["runs"] = {}
else:
    # If runs is a dict (legacy style), validate keys & lengths
    if isinstance(_module_runs, dict):
        for _k in _required_run_keys:
            _v = _module_runs.get(_k)
            if not _v:
                raise Exception(f"Configuration error: 'runs.{_k}' is missing or empty in module config (value={_v})")

        _lens = [len(_module_runs[k]) for k in ["tumour_sample_id", "normal_sample_id", "tumour_seq_type", "tumour_genome_build", "pair_status"]]
        if len(set(_lens)) != 1:
            raise Exception(f"Configuration error: inconsistent lengths in runs lists: tumour_sample_id/normal_sample_id/tumour_seq_type/tumour_genome_build/pair_status -> lengths={_lens}")

        CFG["runs"] = _module_runs

    # If runs is a pandas DataFrame, ensure required columns exist and are non-empty
    elif hasattr(_module_runs, 'columns'):
        missing_cols = [c for c in _required_run_keys if c not in list(_module_runs.columns)]
        if missing_cols:
            raise Exception(f"Configuration error: runs DataFrame missing required columns: {missing_cols}")
        # Ensure no required column is entirely empty
        empty_cols = [c for c in _required_run_keys if _module_runs[c].dropna().shape[0] == 0]
        if empty_cols:
            raise Exception(f"Configuration error: runs DataFrame has empty required columns: {empty_cols}")
        CFG["runs"] = _module_runs

    else:
        # Unknown runs type but present; try to pass it through
        CFG["runs"] = _module_runs

# Return the canonical ploidy (first entry in config ploidy_runs) or a sensible default
def _canonical_ploidy():
    CFG = config["lcr-modules"]["battenberg"]
    pr = CFG["options"].get("ploidy_runs", [])
    if not pr:
        return "1.6-4.8"
    return pr[0]


def _ploidy_runs_or_default():
    pr = _battenberg_CFG["options"].get("ploidy_runs")
    if pr:
        return pr
    return [_canonical_ploidy()]

# Define rules to be run locally when using a compute cluster
localrules:
    _battenberg_all


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _battenberg_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    group: "setup_run"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bam + ".bai", output.bai)
        op.absolute_symlink(input.bam + ".bai", output.crai)

# this process is very fast on bam files and painfully slow on cram files. 
# The result of calc_sex_status.sh is stored in a file to avoid having to rerun it unnecessarily
rule _infer_patient_sex:
    input: 
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: sex_result = CFG["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex"
    resources:
        **CFG["resources"]["infer_sex"]
    log:
        stderr = CFG["logs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}_infer_sex_stderr.log"
    group: "setup_run"
    conda:
        CFG["conda_envs"]["samtools"]
    threads: 8
    shell:
        op.as_one_line("""
        PATH={SCRIPT_PATH}:$PATH; 
        echo "running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr} ;
        calc_sex_status.sh {input.normal_bam} {input.fasta} {wildcards.normal_id} > {output.sex_result} 2>> {log.stderr} &&
        echo "DONE running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" >> {log.stderr} 
        """)


# This rule runs the entire Battenberg pipeline. Eventually we may want to set this rule up to allow re-starting
# of partially completed jobs (e.g. if they run out of RAM and are killed by the cluster, they can automatically retry)
# TODO: this rule needs to be modified to rely on reference_files and allow setup (downloading) of the Battenberg references
rule _run_battenberg:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        sex_result = CFG["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        # Preprocess outputs (preserve these so multiple fits can reuse them)
        ac=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleCounts.tab",
        mb=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantBAF.tab",
        mlrg=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR_gcCorrected.tab",
        mlr=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR.tab",
        nlr=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalLogR.tab",
        nb=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalBAF.tab",
        hap=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_heterozygousMutBAFs_haplotyped.txt"
    log:
        stdout = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stdout.log",
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stderr.log"
    params:
        reference_path = lambda w: _battenberg_CFG["reference_path"][w.genome_build],
        script = CFG["inputs"]["battenberg_script"],
        chr_prefixed = lambda w: _battenberg_CFG["options"]["chr_prefixed_reference"][w.genome_build],
        out_dir = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}"
    conda:
        CFG["conda_envs"]["battenberg"]
    resources:
        **CFG["resources"]["battenberg"]
    threads:
        CFG["threads"]["battenberg"]
    shell:
       op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stdout};
        sex=$(cut -f 4 {input.sex_result}| tail -n 1); 
        echo "setting sex as $sex";
        Rscript --vanilla {params.script} -t {wildcards.tumour_id} 
        -n {wildcards.normal_id} --tb {input.tumour_bam} --nb {input.normal_bam} -f {input.fasta}
        -o {params.out_dir} --sex $sex --reference {params.reference_path} {params.chr_prefixed} --cpu {threads} >> {log.stdout} 2>> {log.stderr} &&  
        echo "DONE {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" >> {log.stdout}; 
        """)


# Fit rule: run Battenberg fitting using precomputed .tab files. This rule can be run multiple
# times per pair with different `ploidy_constraint` wildcard values (format: MIN-MAX).
rule _run_battenberg_fit:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        sex_result = CFG["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        # precomputed per-pair tables
        ac = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleCounts.tab",
        mb = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantBAF.tab",
        mlrg = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR_gcCorrected.tab",
        mlr = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR.tab",
        nlr = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalLogR.tab",
        nb = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalBAF.tab",
        hap = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_heterozygousMutBAFs_haplotyped.txt"

    output:
        refit = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/ploidy_{ploidy_constraint}/{tumour_id}_refit_suggestion.txt",
        sub = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/ploidy_{ploidy_constraint}/{tumour_id}_subclones.txt",
        cp = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/ploidy_{ploidy_constraint}/{tumour_id}_cellularity_ploidy.txt"
    log:
        stdout = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/ploidy_{ploidy_constraint}/{tumour_id}_battenberg.stdout.log",
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/ploidy_{ploidy_constraint}/{tumour_id}_battenberg.stderr.log"
    params:
        script = CFG["inputs"]["battenberg_script"],
        preprocess_dir = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/",
        out_dir = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/ploidy_{ploidy_constraint}",
        reference_path = lambda w: _battenberg_CFG["reference_path"][w.genome_build],
        ploidy_min = lambda w: w.ploidy_constraint.split('-')[0],
        ploidy_max = lambda w: w.ploidy_constraint.split('-')[1]
    conda:
        CFG["conda_envs"]["battenberg"]
    resources:
        **CFG["resources"]["battenberg_fit"]
    threads:
        CFG["threads"]["battenberg_fit"]
    shell:
        op.as_one_line("""
        mkdir -p "{params.out_dir}";
        cp -al {input.ac} {input.mb} {input.mlrg} {input.mlr} {input.nlr} {input.nb} {input.hap} {params.out_dir}/;
        for f in {params.preprocess_dir}/{wildcards.tumour_id}_alleleFrequencies_chr*.txt; do
          if [[ -e "$f" ]]; then
            ln -sf "$f" {params.out_dir}/;
          fi
        done;
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} ploidy {wildcards.ploidy_constraint} on $(hostname) at $(date)" > {log.stdout};
        if [[ $(head -c 4 {input.fasta}) == ">chr" ]]; then chr_prefixed='true'; else chr_prefixed='false'; fi;
        sex=$(cut -f 4 {input.sex_result}| tail -n 1);
        Rscript --vanilla {params.script} \
        -t {wildcards.tumour_id} \
        -n {wildcards.normal_id} \
        --tb $(readlink -f {input.tumour_bam}) \
        --nb $(readlink -f {input.normal_bam}) \
        -f {input.fasta} --reference $(readlink -f {params.reference_path}) \
        -o {params.out_dir} --chr_prefixed_genome $chr_prefixed \
        --sex $sex --ploidy_constraint {wildcards.ploidy_constraint} \
        --skip_allelecount --skip_preprocessing --skip_phasing \
        --cpu {threads} >> {log.stdout} 2>> {log.stderr} && \
        echo "DONE {rule} for {wildcards.tumour_id}--{wildcards.normal_id} ploidy {wildcards.ploidy_constraint} on $(hostname) at $(date)" >> {log.stdout};
        """)


# Convert the subclones.txt (best fit) to igv-friendly SEG files. 
rule _battenberg_to_igv_seg:
    input:
        # use the canonical ploidy run's subclones file
        sub = lambda w: str(rules._run_battenberg_fit.output.sub).replace("{ploidy_constraint}", _canonical_ploidy()),
        cnv2igv = ancient(CFG["inputs"]["cnv2igv"])
    output:
        seg = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.igv.seg"
    log:
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_seg2igv.stderr.log"
    threads: 1
    group: "battenberg_post_process"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        python {input.cnv2igv} --mode battenberg --sample {wildcards.tumour_id} 
        {input.sub} > {output.seg} 2>> {log.stderr}
        """)


# Fill subclones.txt with empty regions for compatibility with downstream tools
rule _battenberg_fill_subclones:
    input:
        sub = lambda w: str(rules._run_battenberg_fit.output.sub).replace("{ploidy_constraint}", _canonical_ploidy())
    output:
        sub = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.filled.txt"
    log:
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_fill_subclones.stderr.log"
    threads: 1
    group: "battenberg_post_process"
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"],
        script = "fill_segments.sh",
        arm_file = lambda w: "src/chromArm.hg38.bed" if "38" in str({w.genome_build}) else "src/chromArm.grch37.bed",
        blacklist_file = lambda w: "src/blacklisted.hg38.bed" if "38" in str({w.genome_build}) else "src/blacklisted.grch37.bed"
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.path}{params.script}
        {params.path}{params.arm_file}
        {input.sub}
        {params.path}{params.blacklist_file}
        {output.sub}
        {wildcards.tumour_id}
        subclones
        2>> {log.stderr}
        """)


#due to the large number of files (several per chromosome) that are not explicit outputs, do some glob-based cleaning in the output directory
rule _battenberg_cleanup:
    input:
        rules._battenberg_to_igv_seg.output.seg,
        # ensure cleanup waits for all ploidy fit outputs for this pair so tabs can be safely removed
        lambda w: expand(
            str(rules._run_battenberg_fit.output.sub),
            seq_type=[w.seq_type],
            genome_build=[w.genome_build],
            tumour_id=[w.tumour_id],
            normal_id=[w.normal_id],
            ploidy_constraint=_ploidy_runs_or_default()
        )
    output:
        complete = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cleanup_complete.txt"
    group: "battenberg_post_process"
    shell:
        op.as_one_line("""
        d=$(dirname {output});
        for dd in "$d"/ploidy_*; do
          if [[ -d "$dd" ]]; then
            rm -f "$dd"/*impute_input* &&
            rm -f "$dd"/*alleleFrequencies* &&
            rm -f "$dd"/*aplotype* &&
            rm -f "$dd"/*BAFsegmented*;
          fi
        done &&
        touch {output.complete}
        """)


rule _battenberg_normal_fake_vcf:
    input:
        ac = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleCounts.tab"
    output:
        vcf_gz = temp(CFG["dirs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}.normal.fake.vcf.gz"),
        tbi = temp(CFG["dirs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}.normal.fake.vcf.gz.tbi")
    log:
        stderr = CFG["logs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}_fake_vcf.stderr.log"
    params:
        script = CFG["inputs"]["src_dir"] + "battenberg_normal_to_vcf.pl",
        min_dp = lambda w: _battenberg_CFG["options"].get("roh", {}).get("min_dp", 12),
        hom_max_minor = lambda w: _battenberg_CFG["options"].get("roh", {}).get("hom_max_minor", 1),
        min_minor_count = lambda w: _battenberg_CFG["options"].get("roh", {}).get("min_minor_count", 3),
        het_lo = lambda w: _battenberg_CFG["options"].get("roh", {}).get("het_lo", 0.20),
        het_hi = lambda w: _battenberg_CFG["options"].get("roh", {}).get("het_hi", 0.80)
    conda:
        CFG["conda_envs"]["battenberg"]
    threads: 1
    group: "plink"
    shell:
        op.as_one_line("""
        mkdir -p $(dirname {output.vcf_gz});
        perl {params.script} \
        --in {input.ac} \
        --sample {wildcards.normal_id} \
        --out - \
        --minDP {params.min_dp} \
        --homMaxMinor {params.hom_max_minor} \
        --minMinorCount {params.min_minor_count} \
        --hetLo {params.het_lo} \
        --hetHi {params.het_hi} \
        2> {log.stderr} | bgzip -c > {output.vcf_gz} &&
        tabix -p vcf {output.vcf_gz} 2>> {log.stderr}
        """)


rule _battenberg_plink_make_bed:
    input:
        vcf_gz = rules._battenberg_normal_fake_vcf.output.vcf_gz,
        tbi = rules._battenberg_normal_fake_vcf.output.tbi
    output:
        bed = temp(CFG["dirs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}.bed"),
        bim = temp(CFG["dirs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}.bim"),
        fam = temp(CFG["dirs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}.fam")
    log:
        stderr = CFG["logs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}_plink_make_bed.stderr.log"
    params:
        out_prefix = CFG["dirs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}"
    conda:
        CFG["conda_envs"]["plink"]
    threads: 1
    group: "plink"
    shell:
        op.as_one_line("""
        mkdir -p $(dirname {params.out_prefix});
        plink --vcf {input.vcf_gz} --make-bed --out {params.out_prefix} 2> {log.stderr}
        """)


rule _battenberg_plink_roh:
    input:
        bed = rules._battenberg_plink_make_bed.output.bed,
        bim = rules._battenberg_plink_make_bed.output.bim,
        fam = rules._battenberg_plink_make_bed.output.fam
    output:
        hom = CFG["dirs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}.roh.hom"
    log:
        stderr = CFG["logs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}_plink_roh.stderr.log"
    params:
        bfile_prefix = CFG["dirs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}",
        out_prefix = CFG["dirs"]["plink"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{normal_id}.roh",
        homozyg_kb = lambda w: _battenberg_CFG["options"].get("roh", {}).get("homozyg_kb", 1000),
        homozyg_snp = lambda w: _battenberg_CFG["options"].get("roh", {}).get("homozyg_snp", 50)
    conda:
        CFG["conda_envs"]["plink"]
    threads: 1
    group: "plink"
    shell:
        op.as_one_line("""
        plink --bfile {params.bfile_prefix} \
        --homozyg --homozyg-kb {params.homozyg_kb} --homozyg-snp {params.homozyg_snp} \
        --out {params.out_prefix} 2> {log.stderr}
        """)


rule _battenberg_output_roh:
    input:
        hom = rules._battenberg_plink_roh.output.hom
    output:
        hom = CFG["output"]["txt"]["roh"]
    threads: 1
    group: "plink"
    run:
        op.relative_symlink(input.hom, output.hom, in_module=True)


def _battenberg_get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")

# Convert the coordinates of seg file to a different genome build
rule _battenberg_convert_coordinates:
    input:
        battenberg_native = str(rules._battenberg_to_igv_seg.output.seg),
        battenberg_chain = _battenberg_get_chain
    output:
        battenberg_lifted = CFG["dirs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.stderr.log"
    threads: 1
    group: "battenberg_post_process"
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        liftover_minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.liftover_script}
        SEG
        {input.battenberg_native}
        {output.battenberg_lifted}
        {input.battenberg_chain}
        YES
        {params.liftover_minmatch}
        2>> {log.stderr}
        """)

# ensure to request the correct files for each projection and drop wildcards that won't be used downstream
def _battenberg_prepare_projection(wildcards):
    CFG = config["lcr-modules"]["battenberg"]
    tbl = CFG["runs"]
    this_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"].tolist()
    
    prefixed_projections = CFG["options"]["prefixed_projections"]
    non_prefixed_projections = CFG["options"]["non_prefixed_projections"]

    if any(substring in this_genome_build[0] for substring in prefixed_projections):
        hg38_projection = str(rules._battenberg_to_igv_seg.output.seg).replace("{genome_build}", this_genome_build[0])
        grch37_projection = str(rules._battenberg_convert_coordinates.output.battenberg_lifted).replace("{genome_build}", this_genome_build[0])
        # handle the hg19 (prefixed) separately
        if "38" in str(this_genome_build[0]):
            grch37_projection = grch37_projection.replace("{chain}", "hg38ToHg19")
        else:
            grch37_projection = grch37_projection.replace("{chain}", "hg19ToHg38")

    elif any(substring in this_genome_build[0] for substring in non_prefixed_projections):
        grch37_projection = str(rules._battenberg_to_igv_seg.output.seg).replace("{genome_build}", this_genome_build[0])
        hg38_projection = str(rules._battenberg_convert_coordinates.output.battenberg_lifted).replace("{genome_build}", this_genome_build[0])
        # handle the grch38 (non-prefixed) separately
        if "38" in str(this_genome_build[0]):
            hg38_projection = hg38_projection.replace("{chain}", "hg38ToHg19")
        else:
            hg38_projection = hg38_projection.replace("{chain}", "hg19ToHg38")
    else:
        raise AttributeError(f"The specified genome build {this_genome_build[0]} is not specified in the config under options to indicate its chr prefixing.")

    return{
        "grch37_projection": grch37_projection,
        "hg38_projection": hg38_projection
    }


# Fill segments of both native and filled file
rule _battenberg_fill_segments:
    input:
        unpack(_battenberg_prepare_projection)
    output:
        grch37_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"),
        hg38_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg")
    log:
        stderr = CFG["logs"]["fill_regions"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}_fill_segments.stderr.log"
    threads: 1
    group: "battenberg_post_process"
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"]
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        echo "Filling grch37 projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.grch37.bed
        {input.grch37_projection}
        {params.path}src/blacklisted.grch37.bed
        {output.grch37_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        echo "Filling hg38 projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.hg38.bed
        {input.hg38_projection}
        {params.path}src/blacklisted.hg38.bed
        {output.hg38_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        """)


def _battenberg_determine_projection(wildcards):
    CFG = config["lcr-modules"]["battenberg"]
    if any(substring in wildcards.projection for substring in ["hg19", "grch37", "hs37d5"]):
        this_file = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"
    elif any(substring in wildcards.projection for substring in ["hg38", "grch38"]):
        this_file = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg"
    return (this_file)


# Normalize chr prefix of the output file
rule _battenberg_normalize_projection:
    input:
        filled = _battenberg_determine_projection,
        chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
    output:
        projection = CFG["dirs"]["normalize"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    resources:
        **CFG["resources"]["post_battenberg"]
    threads: 1
    group: "battenberg_post_process"
    run:
        # read the main chromosomes file of the projection
        chromosomes = pd.read_csv(input.chrom_file, sep = "\t", names=["chromosome"], header=None)
        # handle chr prefix
        if "chr" in chromosomes["chromosome"][0]:
            seg_open = pd.read_csv(input.filled, sep = "\t")
            chrom = list(seg_open['chrom'])
            # avoid cases of chrchr1 if the prefix already there
            for i in range(len(chrom)):
                if 'chr' not in str(chrom[i]):
                    chrom[i]='chr'+str(chrom[i])
            seg_open.loc[:, 'chrom']=chrom
            seg_open.to_csv(output.projection, sep="\t", index=False)
        else:
            # remove chr prefix
            seg_open = pd.read_csv(input.filled, sep = "\t")
            seg_open["chrom"] = seg_open["chrom"].astype(str).str.replace('chr', '')
            seg_open.to_csv(output.projection, sep="\t", index=False)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _battenberg_output_projection:
    input:
        projection = str(rules._battenberg_normalize_projection.output.projection)
    output:
        projection = CFG["output"]["seg"]["projection"]
    threads: 1
    group: "battenberg_post_process"
    run:
        op.relative_symlink(input.projection, output.projection, in_module = True)

# Symlinks the final output files into the module results directory (under '99-outputs/')
# All plots generated by Battenberg are symlinked using a glob for convenience

rule _battenberg_output_seg:
    input:
        seg = rules._battenberg_to_igv_seg.output.seg,
        sub = rules._battenberg_fill_subclones.output.sub,
        cp = lambda w: str(rules._run_battenberg_fit.output.cp).replace("{ploidy_constraint}", _canonical_ploidy())
    output:
        seg = CFG["output"]["seg"]["original"],
        sub = CFG["output"]["txt"]["subclones"],
        cp = CFG["output"]["txt"]["cell_ploidy"]
    params: 
        batt_dir = CFG["dirs"]["battenberg"] + "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}/ploidy_" + _canonical_ploidy(),
        png_dir = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}"
    group: "battenberg_post_process"
    run:
        plots = glob.glob(params.batt_dir + "/*.png")
        for png in plots:
            bn = os.path.basename(png)
            op.relative_symlink(png, params.png_dir + "/" + bn,in_module=True)
        op.relative_symlink(input.seg, output.seg,in_module=True)
        op.relative_symlink(input.sub, output.sub,in_module=True)
        op.relative_symlink(input.cp, output.cp,in_module=True)

# Generates the target sentinels for each run, which generate the symlinks
rule _battenberg_all:
    input:
        expand(
            [
                rules._battenberg_output_seg.output.sub,
                rules._battenberg_output_seg.output.seg,
                rules._battenberg_cleanup.output.complete,
                rules._battenberg_output_roh.output.hom
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]),
        expand(
            expand(
            [
                str(rules._battenberg_output_projection.output.projection)
            ],
            zip,  # Run expand() with zip(), not product()
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            seq_type=CFG["runs"]["tumour_seq_type"],
            pair_status=CFG["runs"]["pair_status"],
            #repeat the tool name N times in expand so each pair in run is used
            tool=["battenberg"] * len(CFG["runs"]["tumour_sample_id"]),
            allow_missing=True),
            projection=CFG["output"]["requested_projections"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
