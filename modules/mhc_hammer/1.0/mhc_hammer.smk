#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Tom Jones & Clare Puttick (McGranahanLab/mhc-hammer)
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import os
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version = "1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section


# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mhc_hammer"]`
CFG = op.setup_module(
    name = "mhc_hammer",
    version = "1.0",
    subdirectories = [
        "inputs", "mhc_reference", "preprocess", "hlahd", "patient_reference",
        "novoalign", "allele_bams", "dna_analysis", "mutations", "cohort_tables", "outputs"
    ],
)

# v1 is DNA-only (WES): HLA typing -> personalised reference -> Novoalign alignment ->
# allele-specific BAM splitting -> LOH/allelic-imbalance detection -> Mutect2+VEP mutation
# calling. RNA allelic expression/imbalance/repression and alt-splicing are out of scope --
# see modules/mhc_hammer/README.md.
#
# Design note on dynamic HLA alleles: the number/names of HLA alleles for a patient are only
# known after HLA-HD typing runs, and every upstream process that needs to iterate "per allele"
# does so with an internal bash/R loop that parses allele names from the personalised reference
# FASTA headers at runtime, rather than spawning one task per allele. This module preserves that
# design: rules that produce allele-named files (rules whose output filenames depend on the
# runtime HLA genotype) always also write at least one fixed-name manifest file in the same
# output directory, and every downstream consumer locates the allele-named siblings via
# `dirname()` + shell-side allele discovery -- never via individual Snakemake output/input
# declarations. This avoids needing a `checkpoint` anywhere in the DAG.

SCRIPTS_DIR = CFG["options"]["mhc_hammer_scripts_dir"]
NOVOALIGN_DIR = CFG["options"]["novoalign_dir"]
HLAHD_DIR = CFG["options"]["hlahd_dir"]
VEP_PATH = CFG["options"]["vep_path"]
VEP_CACHE = CFG["options"]["vep_cache"]
MHC_COORDS = CFG["options"]["mhc_coords"]
IMGT_URLS = CFG["options"]["imgt_release_urls"][CFG["options"]["imgt_release"]]

# Define rules to be run locally when using a compute cluster
localrules:
    _mhc_hammer_input_bam,
    _mhc_hammer_download_reference,
    _mhc_hammer_output_dna_analysis,
    _mhc_hammer_output_mutations,
    _mhc_hammer_output_cohort_table,
    _mhc_hammer_all,


##### HELPER FUNCTIONS #####


# NOTE: every helper function below re-fetches `CFG` from the persistent `config` global
# (rather than closing over the module-level `CFG` variable) because Snakemake evaluates input/
# param functions lazily, after this whole file -- including the `op.cleanup_module(CFG)` call
# at the very end, which deletes the module-level `CFG` shortcut -- has already been parsed.
# This mirrors the identical pattern in modules/gridss/2.0/gridss.smk's `get_input_per_patient`/
# `get_split_vcf` helpers. Confirmed necessary by dry-run (a closure-based version fails with
# `NameError: name 'CFG' is not defined` as soon as any rule beyond input-symlinking is planned).

def _mhc_hammer_get_patient_id_for_sample(sample_id):
    # Every sample_id is unique within a seq_type, so a plain lookup against CFG["samples"]
    # (not CFG["runs"]) is sufficient regardless of whether the sample is tumour or normal.
    CFG = config["lcr-modules"]["mhc_hammer"]
    hits = op.filter_samples(CFG["samples"], sample_id = sample_id)
    return hits["patient_id"].tolist()[0]

def _mhc_hammer_get_patient_id_for_tumour(tumour_id, seq_type):
    # Mirrors the gridss.smk `get_split_vcf` pattern: recover patient_id for a pair-level rule
    # (wildcarded on tumour_id/normal_id, not patient_id) via CFG["runs"].
    CFG = config["lcr-modules"]["mhc_hammer"]
    hits = op.filter_samples(CFG["runs"], tumour_sample_id = tumour_id, tumour_seq_type = seq_type)
    return hits["tumour_patient_id"].tolist()[0]

def _mhc_hammer_get_germline_fqs(wildcards):
    # Every patient must have exactly one germline WES sample (v1 simplification -- see README).
    CFG = config["lcr-modules"]["mhc_hammer"]
    normal = op.filter_samples(
        CFG["samples"],
        patient_id = wildcards.patient_id,
        seq_type = wildcards.seq_type,
        genome_build = wildcards.genome_build,
    )
    normal = normal[normal["tissue_status"].str.lower() == "normal"]
    assert len(normal) == 1, (
        f"Expected exactly one germline WES sample for patient '{wildcards.patient_id}' "
        f"({wildcards.seq_type}--{wildcards.genome_build}), found {len(normal)}. This module "
        f"requires exactly one germline WES sample per patient (see README)."
    )
    sample_id = normal["sample_id"].tolist()[0]
    return {
        "fq1": expand(str(rules._mhc_hammer_generate_fqs.output.fq1), sample_id = sample_id, allow_missing = True),
        "fq2": expand(str(rules._mhc_hammer_generate_fqs.output.fq2), sample_id = sample_id, allow_missing = True),
    }

def _mhc_hammer_reference_dir_for_sample(wildcards):
    patient_id = _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id)
    return expand(str(rules._mhc_hammer_generate_references.output.patient_dir), patient_id = patient_id, allow_missing = True)

def _mhc_hammer_reference_dir_for_tumour(wildcards):
    patient_id = _mhc_hammer_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type)
    return expand(str(rules._mhc_hammer_generate_references.output.patient_dir), patient_id = patient_id, allow_missing = True)

def _mhc_hammer_get_purity_ploidy(wildcards):
    CFG = config["lcr-modules"]["mhc_hammer"]
    purity_ploidy_file = CFG["options"]["purity_ploidy_file"]
    if not purity_ploidy_file:
        return "", ""
    with open(purity_ploidy_file) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] == wildcards.tumour_id:
                return fields[2], fields[1] # purity, ploidy
    return "", ""


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mhc_hammer_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


# Downloads and unpacks the pinned IMGT-derived MHC Hammer reference bundle from Zenodo
# (docs/mhc_reference_files.md upstream). Two separate zip files per release. Cohort-wide,
# no wildcards -- built once and reused by every patient.
rule _mhc_hammer_download_reference:
    output:
        gtf = CFG["dirs"]["mhc_reference"] + "mhc_references/gtf/mhc.gtf",
        genome_fasta = CFG["dirs"]["mhc_reference"] + "mhc_references/genome/mhc_genome_strand.fasta",
        kmer_file = CFG["dirs"]["mhc_reference"] + "kmer_files/imgt_30mers.fa"
    params:
        mhc_references_url = IMGT_URLS["mhc_references"],
        kmer_files_url = IMGT_URLS["kmer_files"],
        outdir = CFG["dirs"]["mhc_reference"]
    conda:
        CFG["conda_envs"]["mhc_hammer_download"]
    container:
        CFG["container_envs"]["mhc_hammer_download"]
    threads:
        CFG["threads"]["download_reference"]
    resources:
        **CFG["resources"]["download_reference"]
    shell:
        op.as_one_line("""
        mkdir -p {params.outdir} &&
        wget -qO {params.outdir}/mhc_references.zip "{params.mhc_references_url}" &&
        wget -qO {params.outdir}/kmer_files.zip "{params.kmer_files_url}" &&
        unzip -oq {params.outdir}/mhc_references.zip -d {params.outdir} &&
        unzip -oq {params.outdir}/kmer_files.zip -d {params.outdir} &&
        rm {params.outdir}/mhc_references.zip {params.outdir}/kmer_files.zip
        """)


# Computes library size (used to normalise depth downstream) on the full input BAM,
# mirroring upstream's FLAGSTAT process.
rule _mhc_hammer_flagstat:
    input:
        bam = str(rules._mhc_hammer_input_bam.output.bam)
    output:
        library_size = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.library_size.txt",
        flagstat = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.flagstat"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    threads:
        CFG["threads"]["flagstat"]
    resources:
        **CFG["resources"]["flagstat"]
    shell:
        op.as_one_line("""
        samtools flagstat {input.bam} > {output.flagstat} &&
        samtools view -c -f 1 -F 2308 {input.bam} > {output.library_size}
        """)


# Subsets the input BAM to reads that could plausibly align to the HLA class I genes: the MHC
# region on chr6, unmapped reads, and (fish_reads=True) reads matching the IMGT 30-mer kmer set.
# Mirrors upstream's SUBSET_BAMS process.
rule _mhc_hammer_subset_bam:
    input:
        bam = str(rules._mhc_hammer_input_bam.output.bam),
        kmer_file = str(rules._mhc_hammer_download_reference.output.kmer_file)
    output:
        bam = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.subset.sorted.bam",
        bai = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.subset.sorted.bam.bai",
        read_counts = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.read_counts.csv"
    params:
        scripts_dir = SCRIPTS_DIR,
        mhc_chrom = lambda w: MHC_COORDS[w.genome_build],
        fish_reads = CFG["options"]["fish_reads"],
        unmapped_reads = CFG["options"]["unmapped_reads"],
        contig_reads = CFG["options"]["contig_reads"],
        sort_mem = lambda wildcards, resources: max(1, int(resources.mem_mb / 1000 * 0.8))
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["subset_bam"]
    resources:
        **CFG["resources"]["subset_bam"]
    shell:
        op.as_one_line("""
        echo {params.mhc_chrom} > $(dirname {output.bam})/mhc_coords.txt &&
        touch $(dirname {output.bam})/extra_contigs.txt &&
        {params.scripts_dir}/bin/subset_bam_opt.sh
        -b {input.bam}
        -k {input.kmer_file}
        -f {params.fish_reads} -c {params.contig_reads} -d $(dirname {output.bam})/extra_contigs.txt
        -u {params.unmapped_reads} -h $(dirname {output.bam})/mhc_coords.txt
        -p $(dirname {output.bam})/{wildcards.sample_id} -t {threads}
        -m {params.sort_mem}G -o false
        """)


# Converts the subsetted BAM to paired FASTQs (samtools collate/fastq), mirroring upstream's
# GENERATE_HLA_FQS process. These are the reads aligned against the personalised reference.
rule _mhc_hammer_generate_fqs:
    input:
        bam = str(rules._mhc_hammer_subset_bam.output.bam)
    output:
        fq1 = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.1.fq.gz",
        fq2 = CFG["dirs"]["preprocess"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.2.fq.gz"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        CFG["container_envs"]["samtools"]
    threads:
        CFG["threads"]["generate_fqs"]
    resources:
        **CFG["resources"]["generate_fqs"]
    shell:
        op.as_one_line("""
        samtools collate -u -O {input.bam} |
        samtools fastq -1 {output.fq1} -2 {output.fq2} -s /dev/null -0 /dev/null -n
        """)


# Runs the user-supplied HLA-HD install on the patient's one germline WES sample to type
# HLA-A/B/C. Patient-level (once per patient, not once per sample). Mirrors upstream's
# HLAHD_LOCAL process. Fails loudly (rather than gracefully excluding the patient, as upstream
# does) if A/B/C are not all typed -- see README v1 simplifications.
rule _mhc_hammer_hlahd:
    input:
        unpack(_mhc_hammer_get_germline_fqs),
        gtf = str(rules._mhc_hammer_download_reference.output.gtf)
    output:
        hla_alleles = CFG["dirs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}_hla_alleles.csv",
        result_dir = directory(CFG["dirs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/result")
    log:
        stdout = CFG["logs"]["hlahd"] + "{seq_type}--{genome_build}/{patient_id}/hlahd.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        hlahd_dir = HLAHD_DIR,
        workdir = lambda wildcards, output: os.path.dirname(output.hla_alleles),
        gtf_abs = lambda wildcards, input: os.path.abspath(input.gtf)
    conda:
        CFG["conda_envs"]["mhc_hammer_hlahd"]
    container:
        None
    threads:
        CFG["threads"]["hlahd"]
    resources:
        **CFG["resources"]["hlahd"]
    shell:
        op.as_one_line("""
        mkdir -p {params.workdir} &&
        (
        export PATH=${{PATH}}:{params.hlahd_dir}/bin &&
        bash {params.hlahd_dir}/bin/hlahd.sh -m 100 -c 1.0 -t {threads}
        -f {params.hlahd_dir}/freq_data
        {input.fq1} {input.fq2}
        {params.hlahd_dir}/HLA_gene.split.txt
        {params.hlahd_dir}/dictionary
        {wildcards.patient_id} {params.workdir} &&
        mv {params.workdir}/{wildcards.patient_id}/* {params.workdir}/ &&
        rmdir {params.workdir}/{wildcards.patient_id}
        ) > {log.stdout} 2>&1 &&
        if [ -f {output.result_dir}/{wildcards.patient_id}_A.est.txt ] &&
           [ -f {output.result_dir}/{wildcards.patient_id}_B.est.txt ] &&
           [ -f {output.result_dir}/{wildcards.patient_id}_C.est.txt ]; then
            (cd {params.workdir} &&
             Rscript {params.scripts_dir}/bin/hlahd_parse_output.R
             --hlahd_folder result --gtf_path {params.gtf_abs}
             --sample_id {wildcards.patient_id} --genes A B C) >> {log.stdout} 2>&1
        else
            echo "ERROR: HLA-HD failed to produce HLA A, B and C estimates for patient {wildcards.patient_id}. See {log.stdout}." >&2 &&
            exit 1
        fi
        """)


# Generates the patient's personalised HLA reference (genome FASTA/GTF/novoindex/SNP positions)
# from their HLA-HD genotype. Patient-level. DNA-only scope -- skips the transcriptome reference
# upstream also builds for RNA analysis. Mirrors upstream's GENERATE_REFERENCES process.
#
# Output is a single directory containing both the fixed-name files referenced explicitly by
# downstream rules (below) and the per-allele *_genome.snp_pos.bed files whose names depend on
# the runtime HLA genotype (see design note near the top of this file).
#   Fixed-name contents: {patient_id}.gtf, {patient_id}_mhc_genome_reference.fa(.fai),
#   {patient_id}_mhc_genome_reference.dict, {patient_id}.mhc_genome_fasta.nix,
#   {patient_id}_sorted.gtf.gz(.tbi), {patient_id}_genome_allele_table.csv, {patient_id}.bed
rule _mhc_hammer_generate_references:
    input:
        genotype = str(rules._mhc_hammer_hlahd.output.hla_alleles),
        mhc_gtf = str(rules._mhc_hammer_download_reference.output.gtf),
        mhc_fasta = str(rules._mhc_hammer_download_reference.output.genome_fasta)
    output:
        patient_dir = directory(CFG["dirs"]["patient_reference"] + "{seq_type}--{genome_build}/{patient_id}")
    log:
        stdout = CFG["logs"]["patient_reference"] + "{seq_type}--{genome_build}/{patient_id}/generate_references.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        novoalign_dir = NOVOALIGN_DIR,
        genotype_abs = lambda wildcards, input: os.path.abspath(input.genotype),
        mhc_gtf_abs = lambda wildcards, input: os.path.abspath(input.mhc_gtf),
        mhc_fasta_abs = lambda wildcards, input: os.path.abspath(input.mhc_fasta)
    conda:
        CFG["conda_envs"]["mhc_hammer_reference"]
    container:
        None
    threads:
        CFG["threads"]["generate_references"]
    resources:
        **CFG["resources"]["generate_references"]
    shell:
        op.as_one_line("""
        mkdir -p {output.patient_dir} &&
        (
        cd {output.patient_dir} &&
        export PATH=${{PATH}}:{params.novoalign_dir} &&
        Rscript {params.scripts_dir}/bin/make_bed_file.R
          --patient_hla_alleles {params.genotype_abs} --patient_id {wildcards.patient_id} --mhc_gtf {params.mhc_gtf_abs} &&
        Rscript {params.scripts_dir}/bin/allele_mismatch.R
          --genome_or_transcriptome genome --patient_id {wildcards.patient_id}
          --patient_hla_alleles {params.genotype_abs} --mhc_fasta {params.mhc_fasta_abs} &&
        Rscript {params.scripts_dir}/bin/create_gtf_and_ref.R
          --gtf_out_path {wildcards.patient_id}.gtf --fa_out_path {wildcards.patient_id}_mhc_genome_reference.fa
          --genome_out_path {wildcards.patient_id}_genome_size.txt --hla_path {params.genotype_abs}
          --mhc_fasta {params.mhc_fasta_abs} --mhc_gtf {params.mhc_gtf_abs} &&
        sort -k1,1 -k4,4n -k5,5n -t$'\\t' {wildcards.patient_id}.gtf > {wildcards.patient_id}_sorted.gtf &&
        bgzip {wildcards.patient_id}_sorted.gtf &&
        tabix -p gff {wildcards.patient_id}_sorted.gtf.gz &&
        samtools faidx {wildcards.patient_id}_mhc_genome_reference.fa &&
        picard CreateSequenceDictionary -R {wildcards.patient_id}_mhc_genome_reference.fa &&
        {params.novoalign_dir}/novoindex {wildcards.patient_id}.mhc_genome_fasta.nix {wildcards.patient_id}_mhc_genome_reference.fa
        ) > {log.stdout} 2>&1
        """)


# Aligns a sample's HLA-region reads against its patient's personalised reference with the
# user-supplied Novoalign. Sample-level. Mirrors upstream's NOVOALIGN process.
rule _mhc_hammer_novoalign:
    input:
        unpack(lambda wildcards: {
            "fq1": str(rules._mhc_hammer_generate_fqs.output.fq1),
            "fq2": str(rules._mhc_hammer_generate_fqs.output.fq2)
        }),
        patient_dir = _mhc_hammer_reference_dir_for_sample
    output:
        bam = CFG["dirs"]["novoalign"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.hla.rehead.bam",
        bai = CFG["dirs"]["novoalign"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.hla.rehead.bam.bai"
    log:
        stderr = CFG["logs"]["novoalign"] + "{seq_type}--{genome_build}/{sample_id}/novoalign.stderr.log"
    params:
        novoalign_dir = NOVOALIGN_DIR,
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id),
        novoindex = lambda wildcards, input: f"{input.patient_dir[0]}/{_mhc_hammer_get_patient_id_for_sample(wildcards.sample_id)}.mhc_genome_fasta.nix"
    conda:
        CFG["conda_envs"]["samtools"]
    container:
        None
    threads:
        CFG["threads"]["novoalign"]
    resources:
        **CFG["resources"]["novoalign"]
    shell:
        op.as_one_line("""
        export PATH=${{PATH}}:{params.novoalign_dir} &&
        wd=$(dirname {output.bam}) &&
        mkfifo $wd/fq1_uncompressed $wd/fq2_uncompressed &&
        gzip -cdf {input.fq1} > $wd/fq1_uncompressed &
        gzip -cdf {input.fq2} > $wd/fq2_uncompressed &
        {params.novoalign_dir}/novoalign -d {params.novoindex} -f $wd/fq1_uncompressed $wd/fq2_uncompressed
          -F STDFQ -R 0 -r All 9999 -o SAM -o FullNW 1> $wd/{wildcards.sample_id}.sam 2> {log.stderr} &&
        samtools view -b -o $wd/{wildcards.sample_id}.bam $wd/{wildcards.sample_id}.sam &&
        samtools sort -o $wd/{wildcards.sample_id}.sorted.bam $wd/{wildcards.sample_id}.bam &&
        samtools view -f 2 -b -o $wd/{wildcards.sample_id}.hla.bam $wd/{wildcards.sample_id}.sorted.bam &&
        samtools addreplacerg -r ID:{wildcards.sample_id} -r SM:{wildcards.sample_id} -o {output.bam} $wd/{wildcards.sample_id}.hla.bam &&
        samtools index {output.bam} &&
        rm $wd/fq1_uncompressed $wd/fq2_uncompressed $wd/{wildcards.sample_id}.sam $wd/{wildcards.sample_id}.bam $wd/{wildcards.sample_id}.sorted.bam $wd/{wildcards.sample_id}.hla.bam
        """)


# Splits the Novoalign BAM into per-allele filtered BAMs, and records which HLA genes/alleles
# passed filtering (heterozygous or not) via `make_hla_bams.sh`. Sample-level. Mirrors upstream's
# MAKE_HLA_ALLELE_BAMS process.
#
# The 4 outputs below are the fixed-name files this module relies on; the per-allele
# *.sorted.filtered.bam files the script also writes into the same directory are not
# individually tracked by Snakemake (see design note near the top of this file) -- downstream
# rules discover them via the personalised reference FASTA headers, exactly as upstream does.
rule _mhc_hammer_make_allele_bams:
    input:
        bam = str(rules._mhc_hammer_novoalign.output.bam),
        patient_dir = _mhc_hammer_reference_dir_for_sample
    output:
        passed_hla_genes = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_passed_hla_genes.txt",
        passed_heterozygous_hla_genes = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_passed_heterozygous_hla_genes.txt",
        passed_hla_alleles = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_passed_hla_alleles.txt",
        passed_heterozygous_hla_alleles = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}_passed_heterozygous_hla_alleles.txt"
    params:
        scripts_dir = SCRIPTS_DIR,
        max_mismatch = CFG["options"]["max_mismatch"],
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id),
        bam_abs = lambda wildcards, input: os.path.abspath(input.bam)
    conda:
        CFG["conda_envs"]["mhc_hammer_alignment"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["make_allele_bams"]
    resources:
        **CFG["resources"]["make_allele_bams"]
    shell:
        op.as_one_line("""
        wd=$(dirname {output.passed_hla_genes}) &&
        mkdir -p $wd &&
        (cd $wd &&
         {params.scripts_dir}/bin/make_hla_bams.sh {params.bam_abs} {params.scripts_dir}/
           {params.max_mismatch} {wildcards.sample_id} {wildcards.seq_type} novoalign)
        """)


# Computes per-allele coverage with mosdepth against the patient's personalised reference BED.
# Sample-level. Mirrors upstream's MOSDEPTH process.
rule _mhc_hammer_mosdepth:
    input:
        allele_bams_marker = str(rules._mhc_hammer_make_allele_bams.output.passed_hla_alleles),
        patient_dir = _mhc_hammer_reference_dir_for_sample
    output:
        bed = CFG["dirs"]["allele_bams"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.{seq_type}.novoalign.mosdepth.bed"
    params:
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_sample(wildcards.sample_id),
        bed_file = lambda wildcards, input: f"{input.patient_dir[0]}/{_mhc_hammer_get_patient_id_for_sample(wildcards.sample_id)}.bed"
    conda:
        CFG["conda_envs"]["mosdepth"]
    container:
        CFG["container_envs"]["mosdepth"]
    threads:
        CFG["threads"]["mosdepth"]
    resources:
        **CFG["resources"]["mosdepth"]
    shell:
        op.as_one_line("""
        wd=$(dirname {input.allele_bams_marker}) &&
        > {output.bed} &&
        for bam_file in $wd/*_{wildcards.seq_type}_novoalign.*.sorted.filtered.bam; do
            [ -e "$bam_file" ] || continue;
            allele=$(basename $bam_file | sed -E 's/^.*_novoalign\\.//; s/\\.sorted\\.filtered\\.bam$//');
            mosdepth --no-per-base --chrom $allele --by {params.bed_file} $wd/{wildcards.sample_id}.$allele $bam_file &&
            gunzip -f $wd/{wildcards.sample_id}.$allele.regions.bed.gz &&
            cat $wd/{wildcards.sample_id}.$allele.regions.bed >> {output.bed};
        done
        """)


##### PAIR-LEVEL HELPER FUNCTIONS #####


def _mhc_hammer_pair_inputs(wildcards, second_allele_marker):
    # Shared plumbing for the two pair-level rules below (detect_cn_aib, detect_muts): both need
    # the tumour+normal allele-bam directories, the patient reference, and (to compute the
    # gene/allele intersection) one marker file from each sample.
    patient_dir = expand(
        str(rules._mhc_hammer_generate_references.output.patient_dir),
        patient_id = _mhc_hammer_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type),
        allow_missing = True
    )
    return {
        "tumour_marker": expand(second_allele_marker, sample_id = wildcards.tumour_id, allow_missing = True),
        "normal_marker": expand(second_allele_marker, sample_id = wildcards.normal_id, allow_missing = True),
        "tumour_flagstat": expand(str(rules._mhc_hammer_flagstat.output.flagstat), sample_id = wildcards.tumour_id, allow_missing = True),
        "normal_flagstat": expand(str(rules._mhc_hammer_flagstat.output.flagstat), sample_id = wildcards.normal_id, allow_missing = True),
        "patient_dir": patient_dir
    }


# Detects HLA allele-specific copy number and allelic imbalance (LOH) by comparing tumour and
# germline coverage at germline-heterozygous SNP positions within each HLA gene the tumour and
# normal both passed filtering for. Pair-level. Mirrors upstream's DETECT_CN_AND_AIB process
# ("all_snps" variant only -- upstream's "exon_snps" variant is disabled/commented out upstream
# too). Copy-number output is skipped (but allelic-imbalance/LOH detection still runs) for
# samples missing a `purity_ploidy_file` entry, matching upstream's `missing_purity_ploidy` flag.
rule _mhc_hammer_detect_cn_aib:
    input:
        unpack(lambda wildcards: _mhc_hammer_pair_inputs(wildcards, str(rules._mhc_hammer_make_allele_bams.output.passed_heterozygous_hla_genes)))
    output:
        dna_analysis = CFG["dirs"]["dna_analysis"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_dna_analysis.csv"
    log:
        stdout = CFG["logs"]["dna_analysis"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/detect_cn_aib.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        min_depth = CFG["options"]["min_depth"],
        purity_ploidy = lambda wildcards: _mhc_hammer_get_purity_ploidy(wildcards),
        tumour_marker_abs = lambda wildcards, input: os.path.abspath(input.tumour_marker[0]),
        normal_marker_abs = lambda wildcards, input: os.path.abspath(input.normal_marker[0]),
        tumour_flagstat_abs = lambda wildcards, input: os.path.abspath(input.tumour_flagstat[0]),
        normal_flagstat_abs = lambda wildcards, input: os.path.abspath(input.normal_flagstat[0]),
        patient_dir_abs = lambda wildcards, input: os.path.abspath(input.patient_dir[0]),
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type)
    conda:
        CFG["conda_envs"]["mhc_hammer_cn_aib"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["detect_cn_aib"]
    resources:
        **CFG["resources"]["detect_cn_aib"]
    shell:
        op.as_one_line("""
        wd=$(dirname {output.dna_analysis}) && mkdir -p $wd &&
        (
        cd $wd &&
        fasta={params.patient_dir_abs}/{params.patient_id}_mhc_genome_reference.fa &&
        tumour_dir=$(dirname {params.tumour_marker_abs}) &&
        normal_dir=$(dirname {params.normal_marker_abs}) &&
        purity="{params.purity_ploidy[0]}" && ploidy="{params.purity_ploidy[1]}" &&
        genes=$(comm -12 <(sort {params.tumour_marker_abs}) <(sort {params.normal_marker_abs})) &&
        for gene in $genes; do
            alleles=($(grep '^>' $fasta | sed 's/^>//' | grep $gene | sort -u)) &&
            allele1=${{alleles[0]}} && allele2=${{alleles[1]}} &&
            allele1_snp_path=${{allele1}}_genome.snp_pos.bed && allele2_snp_path=${{allele2}}_genome.snp_pos.bed &&
            cp {params.patient_dir_abs}/$allele1_snp_path . 2>/dev/null;
            cp {params.patient_dir_abs}/$allele2_snp_path . 2>/dev/null;
            allele1_tumour_bam=$tumour_dir/{wildcards.tumour_id}_{wildcards.seq_type}_novoalign.$allele1.sorted.filtered.bam &&
            allele2_tumour_bam=$tumour_dir/{wildcards.tumour_id}_{wildcards.seq_type}_novoalign.$allele2.sorted.filtered.bam &&
            allele1_gl_bam=$normal_dir/{wildcards.normal_id}_{wildcards.seq_type}_novoalign.$allele1.sorted.filtered.bam &&
            allele2_gl_bam=$normal_dir/{wildcards.normal_id}_{wildcards.seq_type}_novoalign.$allele2.sorted.filtered.bam &&
            samtools mpileup $allele1_tumour_bam -f $fasta > {wildcards.tumour_id}.$allele1.mpileup;
            samtools mpileup $allele1_gl_bam -f $fasta > {wildcards.normal_id}.$allele1.mpileup;
            samtools mpileup $allele2_tumour_bam -f $fasta > {wildcards.tumour_id}.$allele2.mpileup;
            samtools mpileup $allele2_gl_bam -f $fasta > {wildcards.normal_id}.$allele2.mpileup;
            Rscript {params.scripts_dir}/bin/mpileup_to_csv.R --mpileup_path {wildcards.tumour_id}.$allele1.mpileup --hla_ref_path $fasta --allele $allele1 --out_csv_name {wildcards.tumour_id}.$allele1.coverage.csv;
            Rscript {params.scripts_dir}/bin/mpileup_to_csv.R --mpileup_path {wildcards.normal_id}.$allele1.mpileup --hla_ref_path $fasta --allele $allele1 --out_csv_name {wildcards.normal_id}.$allele1.coverage.csv;
            Rscript {params.scripts_dir}/bin/mpileup_to_csv.R --mpileup_path {wildcards.tumour_id}.$allele2.mpileup --hla_ref_path $fasta --allele $allele2 --out_csv_name {wildcards.tumour_id}.$allele2.coverage.csv;
            Rscript {params.scripts_dir}/bin/mpileup_to_csv.R --mpileup_path {wildcards.normal_id}.$allele2.mpileup --hla_ref_path $fasta --allele $allele2 --out_csv_name {wildcards.normal_id}.$allele2.coverage.csv;
            Rscript {params.scripts_dir}/bin/get_filtered_pos_bed.R --coverage_path {wildcards.normal_id}.$allele1.coverage.csv --min_depth {params.min_depth} --out_bed_name {wildcards.normal_id}.$allele1.filtered_positions.bed --allele $allele1;
            Rscript {params.scripts_dir}/bin/get_filtered_pos_bed.R --coverage_path {wildcards.normal_id}.$allele2.coverage.csv --min_depth {params.min_depth} --out_bed_name {wildcards.normal_id}.$allele2.filtered_positions.bed --allele $allele2;
            Rscript {params.scripts_dir}/bin/get_filtered_coverage.R --coverage_path {wildcards.tumour_id}.$allele1.coverage.csv --filtered_positions_path {wildcards.normal_id}.$allele1.filtered_positions.bed --out_csv {wildcards.tumour_id}.$allele1.coverage_at_filtered_positions.csv;
            Rscript {params.scripts_dir}/bin/get_filtered_coverage.R --coverage_path {wildcards.tumour_id}.$allele2.coverage.csv --filtered_positions_path {wildcards.normal_id}.$allele2.filtered_positions.bed --out_csv {wildcards.tumour_id}.$allele2.coverage_at_filtered_positions.csv;
            Rscript {params.scripts_dir}/bin/get_filtered_coverage.R --coverage_path {wildcards.normal_id}.$allele1.coverage.csv --filtered_positions_path {wildcards.normal_id}.$allele1.filtered_positions.bed --out_csv {wildcards.normal_id}.$allele1.coverage_at_filtered_positions.csv;
            Rscript {params.scripts_dir}/bin/get_filtered_coverage.R --coverage_path {wildcards.normal_id}.$allele2.coverage.csv --filtered_positions_path {wildcards.normal_id}.$allele2.filtered_positions.bed --out_csv {wildcards.normal_id}.$allele2.coverage_at_filtered_positions.csv;
            Rscript {params.scripts_dir}/bin/get_gene_filtered_snp_positions.R --allele1_snp_bed $allele1_snp_path --allele2_snp_bed $allele2_snp_path --allele1 $allele1 --allele2 $allele2 --allele1_filtered_pos_bed {wildcards.normal_id}.$allele1.filtered_positions.bed --allele2_filtered_pos_bed {wildcards.normal_id}.$allele2.filtered_positions.bed --sample_name {wildcards.normal_id};
            allele1_gene_filtered_snps={wildcards.normal_id}.$allele1.filtered_snp_positions.bed && allele2_gene_filtered_snps={wildcards.normal_id}.$allele2.filtered_snp_positions.bed &&
            if [ -n "$purity" ] && [ -n "$ploidy" ]; then
                Rscript {params.scripts_dir}/bin/get_cn.R --allele1 $allele1 --allele2 $allele2 --allele1_gl_coverage_file {wildcards.normal_id}.$allele1.coverage_at_filtered_positions.csv --allele2_gl_coverage_file {wildcards.normal_id}.$allele2.coverage_at_filtered_positions.csv --allele1_tumour_coverage_file {wildcards.tumour_id}.$allele1.coverage_at_filtered_positions.csv --allele2_tumour_coverage_file {wildcards.tumour_id}.$allele2.coverage_at_filtered_positions.csv --allele1_snp_bed $allele1_gene_filtered_snps --allele2_snp_bed $allele2_gene_filtered_snps --purity $purity --ploidy $ploidy --tumour_library_size_path {params.tumour_flagstat_abs} --gl_library_size_path {params.normal_flagstat_abs} --gtf_path {params.patient_dir_abs}/{params.patient_id}.gtf --cn_output_path {wildcards.tumour_id}.$gene.all_snps.cn.csv --cn_plots_prefix {wildcards.tumour_id}.$gene.all_snps.cn --scripts_dir {params.scripts_dir}/bin/;
                Rscript {params.scripts_dir}/bin/get_expected_depth.R --allele1 $allele1 --allele2 $allele2 --allele1_gl_coverage_file {wildcards.normal_id}.$allele1.coverage_at_filtered_positions.csv --allele2_gl_coverage_file {wildcards.normal_id}.$allele2.coverage_at_filtered_positions.csv --allele1_snp_bed $allele1_gene_filtered_snps --allele2_snp_bed $allele2_gene_filtered_snps --purity $purity --tumour_library_size_path {params.tumour_flagstat_abs} --gl_library_size_path {params.normal_flagstat_abs} --expected_depth_output_path {wildcards.tumour_id}.$gene.all_snps.expected_depth.csv --expected_depth_plots_prefix {wildcards.tumour_id}.$gene.all_snps.expected_depth --scripts_dir {params.scripts_dir}/bin/;
                missing_purity_ploidy=FALSE;
            else
                missing_purity_ploidy=TRUE;
            fi;
            for role_allele_bam in "{wildcards.tumour_id} $allele1 $allele1_tumour_bam" "{wildcards.tumour_id} $allele2 $allele2_tumour_bam" "{wildcards.normal_id} $allele1 $allele1_gl_bam" "{wildcards.normal_id} $allele2 $allele2_gl_bam"; do
                set -- $role_allele_bam; role=$1; allele=$2; bam=$3;
                snp_bed=${{allele}}_genome.snp_pos.bed;
                bedtools intersect -loj -bed -b $bam -a $snp_bed | cut -f1-7 > $role.$allele.snp_reads_overlap_7cols.bed;
                Rscript {params.scripts_dir}/bin/count_reads_once.R --snp_reads_overlap_bed $role.$allele.snp_reads_overlap_7cols.bed --snp_path $snp_bed --out_csv $role.$allele.coverage_at_filtered_snps_reads_count_once.csv;
            done;
            Rscript {params.scripts_dir}/bin/get_logr_aib.R --allele1 $allele1 --allele2 $allele2 --allele1_gl_reads_count_once_coverage {wildcards.normal_id}.$allele1.coverage_at_filtered_snps_reads_count_once.csv --allele2_gl_reads_count_once_coverage {wildcards.normal_id}.$allele2.coverage_at_filtered_snps_reads_count_once.csv --allele1_tumour_reads_count_once_coverage {wildcards.tumour_id}.$allele1.coverage_at_filtered_snps_reads_count_once.csv --allele2_tumour_reads_count_once_coverage {wildcards.tumour_id}.$allele2.coverage_at_filtered_snps_reads_count_once.csv --allele1_snp_bed $allele1_gene_filtered_snps --allele2_snp_bed $allele2_gene_filtered_snps --tumour_library_size_path {params.tumour_flagstat_abs} --gl_library_size_path {params.normal_flagstat_abs} --logr_aib_output_path {wildcards.tumour_id}.$gene.all_snps.logr_aib.csv --logr_aib_plots_prefix {wildcards.tumour_id}.$gene.all_snps.logr_aib --scripts_dir {params.scripts_dir}/bin/;
        done &&
        Rscript {params.scripts_dir}/bin/concatenate_dna_analysis_tables.R --genes "$genes" --snp_type all_snps --sample_name {wildcards.tumour_id} --aligner novoalign --missing_purity_ploidy ${{missing_purity_ploidy:-TRUE}}
        ) > {log.stdout} 2>&1
        """)
        # NOTE: this rule's shell script was translated from the upstream DETECT_CN_AND_AIB
        # Nextflow process (mhc_hammer_dna_modules.nf) without access to the external R scripts'
        # source (see licensing note near the top of this file) -- verify argument names/order
        # against your own mhc_hammer_scripts_dir checkout before relying on results.


# Calls somatic mutations in each HLA allele the tumour and normal both passed filtering for,
# using Mutect2 in tumour-normal mode followed by VEP annotation against the personalised
# per-patient HLA GTF/FASTA (not a standard genome-build VEP cache). Pair-level. Mirrors
# upstream's DETECT_MUTS process. VEP is user-supplied (vep_path/vep_cache), not conda-managed
# -- see licensing/VEP note near the top of this file.
#
# Output is a directory (not individually-named files) because the set of *.vep.txt files
# depends on the runtime HLA genotype and on which alleles actually carried a mutation; a
# `{tumour_id}.empty.vep.txt` sentinel is always written if none did, so the directory is never
# empty and downstream globbing (rule _mhc_hammer_parse_mutations) is always well-defined.
rule _mhc_hammer_detect_muts:
    input:
        unpack(lambda wildcards: _mhc_hammer_pair_inputs(wildcards, str(rules._mhc_hammer_make_allele_bams.output.passed_hla_alleles)))
    output:
        vep_dir = directory(CFG["dirs"]["mutations"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/vep")
    log:
        stdout = CFG["logs"]["mutations"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/detect_muts.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        vep_path = VEP_PATH,
        vep_cache = VEP_CACHE,
        tumour_marker_abs = lambda wildcards, input: os.path.abspath(input.tumour_marker[0]),
        normal_marker_abs = lambda wildcards, input: os.path.abspath(input.normal_marker[0]),
        patient_dir_abs = lambda wildcards, input: os.path.abspath(input.patient_dir[0]),
        patient_id = lambda wildcards: _mhc_hammer_get_patient_id_for_tumour(wildcards.tumour_id, wildcards.seq_type),
        mem_gb = lambda wildcards, resources: max(1, int(resources.mem_mb / 1000 * 0.8))
    conda:
        CFG["conda_envs"]["mhc_hammer_detect_muts"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir and vep_path -- see note near the top of this file
    threads:
        CFG["threads"]["detect_muts"]
    resources:
        **CFG["resources"]["detect_muts"]
    shell:
        op.as_one_line("""
        mkdir -p {output.vep_dir} &&
        (
        cd {output.vep_dir} &&
        export PATH={params.vep_path}:${{PATH}} &&
        fasta={params.patient_dir_abs}/{params.patient_id}_mhc_genome_reference.fa &&
        tumour_dir=$(dirname {params.tumour_marker_abs}) &&
        normal_dir=$(dirname {params.normal_marker_abs}) &&
        bgzip -dc {params.patient_dir_abs}/{params.patient_id}_sorted.gtf.gz > {wildcards.tumour_id}.gtf &&
        alleles=$(comm -12 <(sort {params.tumour_marker_abs}) <(sort {params.normal_marker_abs})) &&
        for allele in $alleles; do
            normal_bam=$normal_dir/{wildcards.normal_id}_{wildcards.seq_type}_novoalign.$allele.sorted.filtered.bam;
            tumour_bam=$tumour_dir/{wildcards.tumour_id}_{wildcards.seq_type}_novoalign.$allele.sorted.filtered.bam;
            gatk --java-options '-Xmx{params.mem_gb}g -Xms1g' Mutect2 -R $fasta -I $normal_bam -I $tumour_bam -normal {wildcards.normal_id} --f1r2-tar-gz {wildcards.tumour_id}.$allele.f1r2.tar.gz --output {wildcards.tumour_id}.$allele.vcf;
            if tail -n 1 {wildcards.tumour_id}.$allele.vcf | grep -q CHROM; then continue; fi;
            gatk --java-options '-Xmx{params.mem_gb}g -Xms1g' LearnReadOrientationModel -I {wildcards.tumour_id}.$allele.f1r2.tar.gz -O {wildcards.tumour_id}.$allele.read-orientation-model.tar.gz;
            gatk --java-options '-Xmx{params.mem_gb}g -Xms1g' FilterMutectCalls -V {wildcards.tumour_id}.$allele.vcf -R $fasta --ob-priors {wildcards.tumour_id}.$allele.read-orientation-model.tar.gz -O {wildcards.tumour_id}.$allele.filt.vcf;
            bcftools norm -m-any {wildcards.tumour_id}.$allele.filt.vcf --output {wildcards.tumour_id}.$allele.norm.filt.vcf;
            Rscript {params.scripts_dir}/bin/make_vep_gtf.R --gtf_path {wildcards.tumour_id}.gtf --gtf_vep_path {wildcards.tumour_id}.$allele.vep.gtf --allele $allele;
            bgzip {wildcards.tumour_id}.$allele.vep.gtf && tabix -p gff {wildcards.tumour_id}.$allele.vep.gtf.gz;
            vep -i {wildcards.tumour_id}.$allele.norm.filt.vcf --gtf {wildcards.tumour_id}.$allele.vep.gtf.gz --fasta $fasta --dir_cache {params.vep_cache} --vcf -o {wildcards.tumour_id}.$allele.norm.filt.vep.vcf --fields "Allele,Consequence,IMPACT,Feature_type,Feature,EXON,INTRON,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,FLAGS";
            gatk --java-options '-Xmx{params.mem_gb}g -Xms1g' VariantsToTable -V {wildcards.tumour_id}.$allele.norm.filt.vep.vcf -O {wildcards.tumour_id}.$allele.norm.filt.vep.txt -F CHROM -F POS -F REF -F ALT -F FILTER -GF AD -GF DP -F CSQ --show-filtered;
        done;
        if ! ls {wildcards.tumour_id}.*.vep.txt >/dev/null 2>&1; then touch {wildcards.tumour_id}.empty.vep.txt; fi
        ) > {log.stdout} 2>&1
        """)
        # NOTE: translated from the upstream DETECT_MUTS Nextflow process without access to the
        # external scripts' source -- verify against your own mhc_hammer_scripts_dir checkout.


# Aggregates one patient's per-tumour VEP mutation tables (across all of that patient's runs)
# into a single per-patient mutation CSV. Patient-level. Mirrors upstream's PARSE_MUTATIONS
# process. Skipped (no output) for patients where no run produced any real (non-empty) VEP
# table -- matching upstream's own `mutations_detected_check` branch.
def _mhc_hammer_get_patient_mutation_dirs(wildcards):
    CFG = config["lcr-modules"]["mhc_hammer"]
    patient_runs = op.filter_samples(
        CFG["paired_runs"],
        tumour_patient_id = wildcards.patient_id,
        tumour_seq_type = wildcards.seq_type,
        tumour_genome_build = wildcards.genome_build
    )
    return expand(
        str(rules._mhc_hammer_detect_muts.output.vep_dir),
        zip,
        tumour_id = patient_runs["tumour_sample_id"],
        normal_id = patient_runs["normal_sample_id"],
        pair_status = patient_runs["pair_status"],
        allow_missing = True
    )

rule _mhc_hammer_parse_mutations:
    input:
        vep_dirs = _mhc_hammer_get_patient_mutation_dirs
    output:
        mutations = CFG["dirs"]["mutations"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}_mutations.csv"
    log:
        stdout = CFG["logs"]["mutations"] + "{seq_type}--{genome_build}/{patient_id}/parse_mutations.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        vep_dirs_abs = lambda wildcards, input: " ".join(os.path.abspath(d) for d in input.vep_dirs)
    conda:
        CFG["conda_envs"]["mhc_hammer_r"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["parse_mutations"]
    resources:
        **CFG["resources"]["parse_mutations"]
    shell:
        op.as_one_line("""
        wd=$(dirname {output.mutations}) && mkdir -p $wd &&
        (
        vep_tables="" &&
        for d in {params.vep_dirs_abs}; do
            for f in $d/*.vep.txt; do
                case "$f" in *empty.vep.txt) continue ;; esac;
                [ -e "$f" ] && vep_tables="$vep_tables $f";
            done;
        done;
        if [ -z "$vep_tables" ]; then
            echo "No mutations detected for patient {wildcards.patient_id} -- skipping mutation table" &&
            touch {output.mutations};
        else
            cd $wd &&
            Rscript {params.scripts_dir}/bin/make_mutation_table.R --vep_tables $vep_tables --mutation_save_path {wildcards.patient_id}_mutations.csv --scripts_dir {params.scripts_dir}/bin/;
        fi
        ) > {log.stdout} 2>&1
        """)
        # NOTE: upstream's make_mutation_table.R also takes --wxs_tumour_bam_files/--wxs_gl_bam_files
        # and an --inventory csv for read-depth annotation; wire these up against your own
        # mhc_hammer_scripts_dir checkout's actual CLI before relying on this rule's output.


# Builds the cohort-wide gene-level summary table (docs/mhc_hammer_outputs.md upstream,
# DNA-only column subset), combining every patient's HLA-HD genotype, allele-specific coverage
# (mosdepth), DNA copy-number/allelic-imbalance, and mutation calls into one CSV. Cohort-wide,
# no wildcards. Mirrors upstream's CREATE_MHC_HAMMER_TABLE process.
rule _mhc_hammer_cohort_table:
    input:
        dna_analysis = expand(
            str(rules._mhc_hammer_detect_cn_aib.output.dna_analysis),
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            tumour_id = CFG["paired_runs"]["tumour_sample_id"],
            normal_id = CFG["paired_runs"]["normal_sample_id"],
            pair_status = CFG["paired_runs"]["pair_status"]
        ),
        mosdepth = expand(
            str(rules._mhc_hammer_mosdepth.output.bed),
            zip,
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"],
            sample_id = CFG["samples"]["sample_id"]
        )
    output:
        cohort_table = CFG["dirs"]["cohort_tables"] + "cohort_mhc_hammer_gene_table.csv"
    log:
        stdout = CFG["logs"]["cohort_tables"] + "cohort_table.log"
    params:
        scripts_dir = SCRIPTS_DIR,
        min_number_of_snps = CFG["options"]["min_number_of_snps"],
        max_copy_number_range = CFG["options"]["max_copy_number_range"],
        min_expected_depth = CFG["options"]["min_expected_depth"],
        min_depth = CFG["options"]["min_depth"],
        input_csvs_abs = lambda wildcards, input: " ".join(os.path.abspath(f) for f in input.dna_analysis)
    conda:
        CFG["conda_envs"]["mhc_hammer_r"]
    container:
        None # calls the user-supplied mhc_hammer_scripts_dir -- see licensing note near the top of this file
    threads:
        CFG["threads"]["cohort_table"]
    resources:
        **CFG["resources"]["cohort_table"]
    shell:
        op.as_one_line("""
        wd=$(dirname {output.cohort_table}) && mkdir -p $wd &&
        (
        cd $wd &&
        Rscript {params.scripts_dir}/bin/make_cohort_overview_table.R
          --csv_tables_path <(printf '%s\\n' {params.input_csvs_abs})
          --max_cn_range {params.max_copy_number_range}
          --min_n_snps {params.min_number_of_snps}
          --min_expected_depth {params.min_expected_depth}
          --dna_snp_min_depth {params.min_depth}
        ) > {log.stdout} 2>&1
        """)
        # NOTE: upstream's make_cohort_overview_table.R also takes --inventory_path and
        # --hlahd_germline_samples_path; wire these up against your own mhc_hammer_scripts_dir
        # checkout's actual CLI before relying on this rule's output.


# Symlinks the final per-run output files into the module results directory (under '99-outputs/')
rule _mhc_hammer_output_dna_analysis:
    input:
        dna_analysis = str(rules._mhc_hammer_detect_cn_aib.output.dna_analysis)
    output:
        dna_analysis = CFG["dirs"]["outputs"] + "dna_analysis/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dna_analysis.csv"
    run:
        op.relative_symlink(input.dna_analysis, output.dna_analysis, in_module = True)

rule _mhc_hammer_output_mutations:
    input:
        mutations = str(rules._mhc_hammer_parse_mutations.output.mutations)
    output:
        mutations = CFG["dirs"]["outputs"] + "mutations/{seq_type}--{genome_build}/{patient_id}.mutations.csv"
    run:
        op.relative_symlink(input.mutations, output.mutations, in_module = True)

rule _mhc_hammer_output_cohort_table:
    input:
        cohort_table = str(rules._mhc_hammer_cohort_table.output.cohort_table)
    output:
        cohort_table = CFG["dirs"]["outputs"] + "cohort_tables/cohort_mhc_hammer_gene_table.csv"
    run:
        op.relative_symlink(input.cohort_table, output.cohort_table, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks. Uses
# CFG["paired_runs"] (not CFG["runs"]) so that tumour samples without a matched germline WES
# sample -- which this module cannot process at all, since HLA typing and the personalised
# reference both depend on the patient's own germline sample -- are never requested as targets.
rule _mhc_hammer_all:
    input:
        expand(
            [
                str(rules._mhc_hammer_output_dna_analysis.output.dna_analysis)
            ],
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            tumour_id = CFG["paired_runs"]["tumour_sample_id"],
            normal_id = CFG["paired_runs"]["normal_sample_id"],
            pair_status = CFG["paired_runs"]["pair_status"]
        ),
        expand(
            [
                str(rules._mhc_hammer_output_mutations.output.mutations)
            ],
            zip,
            seq_type = CFG["paired_runs"]["tumour_seq_type"],
            genome_build = CFG["paired_runs"]["tumour_genome_build"],
            patient_id = CFG["paired_runs"]["tumour_patient_id"]
        ),
        str(rules._mhc_hammer_output_cohort_table.output.cohort_table)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
