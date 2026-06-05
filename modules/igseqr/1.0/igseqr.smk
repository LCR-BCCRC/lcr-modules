#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["igseqr"]`
CFG = op.setup_module(
    name = "igseqr",
    version = "1.0",
    subdirectories = ["inputs", "hisat2", "bam2fastq", "trinity", "blastn", "kallisto", "reports", "outputs"],
)

assert len(CFG["inputs"]["hisat_ref_url"]) == 1, (
    "Config 'hisat_ref_url' must contain exactly one entry (version_string: url). "
    "To use a different index, replace the key and URL."
)

HISAT_REF_VERSION = list(CFG["inputs"]["hisat_ref_url"].keys())[0]
HISAT_REF_URL     = CFG["inputs"]["hisat_ref_url"][HISAT_REF_VERSION]
HISAT_REF_DIR     = CFG["dirs"]["inputs"] + "hisat_ref"

CHAINS = CFG["options"]["chains"]
if isinstance(CHAINS, str):
    CHAINS = CHAINS.split()

assert all(c in {"IGH", "IGKL"} for c in CHAINS), (
    "Config 'options.chains' must contain 'IGH', 'IGKL', or both."
)

# IG loci for read extraction (Ensembl-style chromosome names, no 'chr' prefix).
# All three loci are always extracted; BLAST separates chains in a downstream rule.
IG_REGIONS = "14:105550000-106900000 2:88697000-92240000 22:22005000-23590000"

# Define rules to be run locally when using a compute cluster
localrules:
    _igseqr_input_fastq,
    _igseqr_get_imgt_db,
    _igseqr_output_files,
    _igseqr_all,

IMGT_DB_DIR = CFG["dirs"]["inputs"] + "imgt_db"


##### RULES #####


# Downloads the HISAT2 splice-aware index specified in hisat_ref_url.
# The tarball is extracted under {inputs_dir}/hisat_ref/ and the index prefix
# (path to the *.ht2 files without the .N.ht2 suffix) is written to the sentinel
# file so downstream rules can discover it regardless of the archive's internal layout.
rule _igseqr_get_hisat_ref:
    params:
        ref_dir = HISAT_REF_DIR,
        version = HISAT_REF_VERSION,
        url     = HISAT_REF_URL,
    output:
        complete = HISAT_REF_DIR + "/hisat2_" + HISAT_REF_VERSION + "_downloaded.success"
    log:
        stdout = CFG["logs"]["hisat2"] + "get_hisat_ref.stdout.log",
        stderr = CFG["logs"]["hisat2"] + "get_hisat_ref.stderr.log",
    resources:
        **CFG["resources"]["get_hisat_ref"]
    shell:
        '''
        wget -O {params.ref_dir}/{params.version}.tar.gz "{params.url}" \
            > {log.stdout} 2> {log.stderr}
        tar -xzf {params.ref_dir}/{params.version}.tar.gz -C {params.ref_dir}/ \
            >> {log.stdout} 2>> {log.stderr}
        rm {params.ref_dir}/{params.version}.tar.gz
        find {params.ref_dir} -name "*.1.ht2" | head -1 | sed 's/\.1\.ht2$//' \
            > {output.complete}
        '''


# Symlinks the input FASTQ pair into the module results directory (under '00-inputs/')
rule _igseqr_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"],
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R2.fastq.gz",
    run:
        op.absolute_symlink(input.fastq_1, output.fastq_1)
        op.absolute_symlink(input.fastq_2, output.fastq_2)


# Aligns paired-end FASTQs to the genome with HISAT2.
# The sorted BAM and its index are used by _igseqr_filter_reads for region extraction.
rule _igseqr_hisat2_align:
    input:
        fastq_1         = str(rules._igseqr_input_fastq.output.fastq_1),
        fastq_2         = str(rules._igseqr_input_fastq.output.fastq_2),
        fastq_1_real    = CFG["inputs"]["sample_fastq_1"],  # Prevent premature deletion of fastqs marked as temp
        fastq_2_real    = CFG["inputs"]["sample_fastq_2"],
        hisat_ref_ready = str(rules._igseqr_get_hisat_ref.output.complete),
    output:
        bam = temp(CFG["dirs"]["hisat2"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}.aligned.sorted.bam"),
        bai = temp(CFG["dirs"]["hisat2"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}.aligned.sorted.bam.bai"),
    log:
        stdout = CFG["logs"]["hisat2"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/hisat2_align.stdout.log",
        stderr = CFG["logs"]["hisat2"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/hisat2_align.stderr.log",
    params:
        opts = CFG["options"]["hisat2_align"],
    threads:
        CFG["threads"]["hisat2_align"]
    resources:
        **CFG["resources"]["hisat2_align"]
    conda:
        CFG["conda_envs"]["igseqr"]
    container:
        CFG["container_envs"]["igseqr"]
    shell:
        '''
        HISAT_PREFIX=$(cat {input.hisat_ref_ready})
        hisat2 -p {threads} --phred33 -t \
            -x $HISAT_PREFIX \
            -1 {input.fastq_1} -2 {input.fastq_2} \
            {params.opts} \
            2> {log.stderr} \
        | samtools view -@ {threads} -bS - 2>> {log.stderr} \
        | samtools sort -@ {threads} -o {output.bam} 2>> {log.stderr}
        hisat2_exit=${{PIPESTATUS[0]}}
        [ "$hisat2_exit" -eq 0 ] || exit "$hisat2_exit"
        samtools index -@ {threads} {output.bam} > {log.stdout} 2>> {log.stderr}
        '''


# Extracts reads overlapping IG loci (IGH, IGK, IGL) plus all unmapped reads,
# then converts to a paired FASTQ for Trinity assembly.
rule _igseqr_filter_reads:
    input:
        bam = str(rules._igseqr_hisat2_align.output.bam),
        bai = str(rules._igseqr_hisat2_align.output.bai),
    output:
        fastq_1 = temp(CFG["dirs"]["bam2fastq"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}.IG.R1.fastq.gz"),
        fastq_2 = temp(CFG["dirs"]["bam2fastq"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}.IG.R2.fastq.gz"),
    log:
        stdout = CFG["logs"]["bam2fastq"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/filter_reads.stdout.log",
        stderr = CFG["logs"]["bam2fastq"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/filter_reads.stderr.log",
    params:
        regions      = IG_REGIONS,
        tmp_unmapped = CFG["dirs"]["bam2fastq"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_tmp_unmapped.bam",
        tmp_ig_loci  = CFG["dirs"]["bam2fastq"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_tmp_ig_loci.bam",
        tmp_merged   = CFG["dirs"]["bam2fastq"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_tmp_merged.bam",
    threads:
        CFG["threads"]["filter_reads"]
    resources:
        **CFG["resources"]["filter_reads"]
    conda:
        CFG["conda_envs"]["igseqr"]
    container:
        CFG["container_envs"]["igseqr"]
    shell:
        '''
        samtools view -@ {threads} -b -f 4 {input.bam} \
            -o {params.tmp_unmapped} > {log.stdout} 2> {log.stderr}
        samtools view -@ {threads} -b {input.bam} {params.regions} \
            -o {params.tmp_ig_loci} >> {log.stdout} 2>> {log.stderr}
        samtools merge -f {params.tmp_merged} \
            {params.tmp_unmapped} {params.tmp_ig_loci} >> {log.stdout} 2>> {log.stderr}
        samtools sort -n -@ {threads} {params.tmp_merged} 2>> {log.stderr} \
        | samtools fastq -@ {threads} -n -c 6 \
            -1 {output.fastq_1} -2 {output.fastq_2} \
            -0 /dev/null -s /dev/null >> {log.stdout} 2>> {log.stderr}
        rm {params.tmp_unmapped} {params.tmp_ig_loci} {params.tmp_merged}
        '''


# De novo assembles IG reads with Trinity. Runs once per sample across all chains;
# BLAST in the next step separates heavy and light chain transcripts.
rule _igseqr_trinity_assemble:
    input:
        fastq_1 = str(rules._igseqr_filter_reads.output.fastq_1),
        fastq_2 = str(rules._igseqr_filter_reads.output.fastq_2),
    output:
        fasta = CFG["dirs"]["trinity"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_Trinity.Trinity.fasta",
    log:
        stdout = CFG["logs"]["trinity"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/trinity.stdout.log",
        stderr = CFG["logs"]["trinity"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/trinity.stderr.log",
    params:
        # Trinity outputs {trinity_dir}.Trinity.fasta when using --full_cleanup
        trinity_dir = CFG["dirs"]["trinity"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_Trinity",
        max_mem     = lambda wildcards, resources: str(int(resources.mem_mb/1000)) + "G",
        opts        = CFG["options"]["trinity_assemble"],
    threads:
        CFG["threads"]["trinity_assemble"]
    resources:
        **CFG["resources"]["trinity_assemble"]
    conda:
        CFG["conda_envs"]["igseqr"]
    container:
        CFG["container_envs"]["trinity"]
    shell:
        op.as_one_line("""
        Trinity
        --seqType fq
        --left {input.fastq_1}
        --right {input.fastq_2}
        --output {params.trinity_dir}
        --max_memory {params.max_mem}
        --min_contig_length 500
        --full_cleanup
        --no_normalize_reads
        --CPU {threads}
        {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)


# Downloads the pre-built IMGT BLAST databases from the IgSeqR GitHub repository.
# Stored under {inputs_dir}/imgt_db/{chain}/ and valid for conda and container runs.
rule _igseqr_get_imgt_db:
    params:
        out_dir = IMGT_DB_DIR,
    output:
        complete = IMGT_DB_DIR + "/imgt_db_downloaded.success"
    log:
        stdout = CFG["logs"]["blastn"] + "get_imgt_db.stdout.log",
        stderr = CFG["logs"]["blastn"] + "get_imgt_db.stderr.log",
    resources:
        **CFG["resources"]["get_imgt_db"]
    shell:
        '''
        wget https://github.com/ForconiLab/IgSeqR/archive/refs/heads/main.tar.gz \
            -O {params.out_dir}/igseqr_main.tar.gz > {log.stdout} 2> {log.stderr}
        tar -xzf {params.out_dir}/igseqr_main.tar.gz \
            -C {params.out_dir} \
            --strip-components=3 \
            IgSeqR-main/data/IMGT >> {log.stdout} 2>> {log.stderr}
        rm {params.out_dir}/igseqr_main.tar.gz
        touch {output.complete}
        '''


# BLASTs the Trinity assembly against the IMGT germline database for each chain.
rule _igseqr_blastn:
    input:
        fasta        = str(rules._igseqr_trinity_assemble.output.fasta),
        imgt_db_ready = str(rules._igseqr_get_imgt_db.output.complete),
    output:
        blast = temp(CFG["dirs"]["blastn"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}_blast.outfmt6"),
    log:
        stdout = CFG["logs"]["blastn"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/blastn.stdout.log",
        stderr = CFG["logs"]["blastn"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/blastn.stderr.log",
    params:
        db   = lambda wildcards: IMGT_DB_DIR + f"/{wildcards.chain}/{wildcards.chain}.fasta",
        opts = CFG["options"]["blastn"],
    wildcard_constraints:
        chain = "|".join(CHAINS),
    threads:
        CFG["threads"]["blastn"]
    resources:
        **CFG["resources"]["blastn"]
    conda:
        CFG["conda_envs"]["igseqr"]
    container:
        CFG["container_envs"]["igseqr"]
    shell:
        '''
        if [ ! -s {input.fasta} ]; then
            touch {output.blast}
        else
            blastn \
                -db {params.db} \
                -query {input.fasta} \
                -outfmt 6 \
                -num_threads {threads} \
                -out {output.blast} \
                {params.opts} \
                > {log.stdout} 2> {log.stderr}
        fi
        '''


# Extracts the Trinity contigs that had BLAST hits against the chain's IMGT database.
rule _igseqr_extract_transcripts:
    input:
        fasta = str(rules._igseqr_trinity_assemble.output.fasta),
        blast = str(rules._igseqr_blastn.output.blast),
    output:
        transcripts = CFG["dirs"]["blastn"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}_transcripts.fasta",
    log:
        stdout = CFG["logs"]["blastn"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/extract_transcripts.stdout.log",
        stderr = CFG["logs"]["blastn"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/extract_transcripts.stderr.log",
    wildcard_constraints:
        chain = "|".join(CHAINS),
    resources:
        **CFG["resources"]["extract_transcripts"]
    conda:
        CFG["conda_envs"]["igseqr"]
    container:
        CFG["container_envs"]["igseqr"]
    shell:
        '''
        samtools faidx {input.fasta} 2> {log.stderr}
        HITS=$(cut -f1 {input.blast} | sort -u)
        if [ -n "$HITS" ]; then
            echo "$HITS" | xargs samtools faidx {input.fasta} \
                > {output.transcripts} 2>> {log.stderr}
        else
            touch {output.transcripts}
        fi
        echo "Done" > {log.stdout}
        '''


# Builds a kallisto index from the chain-specific transcripts FASTA.
rule _igseqr_kallisto_index:
    input:
        transcripts = str(rules._igseqr_extract_transcripts.output.transcripts),
    output:
        idx = temp(CFG["dirs"]["kallisto"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}.kallisto.idx"),
    log:
        stdout = CFG["logs"]["kallisto"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/kallisto_index.stdout.log",
        stderr = CFG["logs"]["kallisto"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/kallisto_index.stderr.log",
    wildcard_constraints:
        chain = "|".join(CHAINS),
    resources:
        **CFG["resources"]["kallisto_index"]
    conda:
        CFG["conda_envs"]["igseqr"]
    container:
        CFG["container_envs"]["igseqr"]
    shell:
        '''
        if [ ! -s {input.transcripts} ]; then
            touch {output.idx}
        else
            kallisto index -i {output.idx} {input.transcripts} > {log.stdout} 2> {log.stderr}
        fi
        '''


# Quantifies transcript abundance with kallisto using the IG-filtered reads.
rule _igseqr_kallisto_quant:
    input:
        idx     = str(rules._igseqr_kallisto_index.output.idx),
        fastq_1 = str(rules._igseqr_filter_reads.output.fastq_1),
        fastq_2 = str(rules._igseqr_filter_reads.output.fastq_2),
    output:
        abundance = CFG["dirs"]["kallisto"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/abundance.tsv",
    log:
        stdout = CFG["logs"]["kallisto"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/kallisto_quant.stdout.log",
        stderr = CFG["logs"]["kallisto"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/kallisto_quant.stderr.log",
    params:
        out_dir = CFG["dirs"]["kallisto"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}",
    wildcard_constraints:
        chain = "|".join(CHAINS),
    threads:
        CFG["threads"]["kallisto_quant"]
    resources:
        **CFG["resources"]["kallisto_quant"]
    conda:
        CFG["conda_envs"]["igseqr"]
    container:
        CFG["container_envs"]["igseqr"]
    shell:
        '''
        if [ ! -s {input.idx} ]; then
            printf 'target_id\tlength\teff_length\test_counts\ttpm\n' > {output.abundance}
        else
            kallisto quant \
                -i {input.idx} \
                -o {params.out_dir} \
                -t {threads} \
                {input.fastq_1} {input.fastq_2} \
                > {log.stdout} 2> {log.stderr}
        fi
        '''


# Generates TSV reports and a TPM-filtered FASTA from kallisto abundance and transcripts.
rule _igseqr_make_report:
    input:
        abundance   = str(rules._igseqr_kallisto_quant.output.abundance),
        transcripts = str(rules._igseqr_extract_transcripts.output.transcripts),
    output:
        report          = CFG["dirs"]["reports"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}_report.tsv",
        dominant_report = CFG["dirs"]["reports"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}_dominant_report.tsv",
        tpm_fasta       = CFG["dirs"]["reports"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}_TPM_filtered.fasta",
    log:
        stdout = CFG["logs"]["reports"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/make_report.stdout.log",
        stderr = CFG["logs"]["reports"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{chain}/make_report.stderr.log",
    params:
        script    = CFG["scripts"]["igseqr_report"],
        sample_id = "{sample_id}",
    wildcard_constraints:
        chain = "|".join(CHAINS),
    resources:
        **CFG["resources"]["make_report"]
    conda:
        CFG["conda_envs"]["igseqr"]
    container:
        CFG["container_envs"]["igseqr"]
    shell:
        op.as_one_line("""
        python {params.script}
        --abundance {input.abundance}
        --transcripts {input.transcripts}
        --sample_id {params.sample_id}
        --report {output.report}
        --dominant_report {output.dominant_report}
        --tpm_fasta {output.tpm_fasta}
        > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _igseqr_output_files:
    input:
        transcripts_fasta = str(rules._igseqr_extract_transcripts.output.transcripts),
        report            = str(rules._igseqr_make_report.output.report),
        dominant_report   = str(rules._igseqr_make_report.output.dominant_report),
        tpm_fasta         = str(rules._igseqr_make_report.output.tpm_fasta),
    output:
        transcripts_fasta = CFG["dirs"]["outputs"] + "fasta/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_transcripts.fasta",
        report            = CFG["dirs"]["outputs"] + "tsv/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_report.tsv",
        dominant_report   = CFG["dirs"]["outputs"] + "tsv/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_dominant_report.tsv",
        tpm_fasta         = CFG["dirs"]["outputs"] + "fasta/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_TPM_filtered.fasta",
    wildcard_constraints:
        chain = "|".join(CHAINS),
    run:
        op.relative_symlink(input.transcripts_fasta, output.transcripts_fasta, in_module=True)
        op.relative_symlink(input.report,            output.report,            in_module=True)
        op.relative_symlink(input.dominant_report,   output.dominant_report,   in_module=True)
        op.relative_symlink(input.tpm_fasta,         output.tpm_fasta,         in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _igseqr_all:
    input:
        expand(
            expand(
                [
                    rules._igseqr_output_files.output.transcripts_fasta,
                    rules._igseqr_output_files.output.report,
                ],
                chain=CHAINS,
                allow_missing=True,
            ),
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            sample_id=CFG["samples"]["sample_id"],
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)