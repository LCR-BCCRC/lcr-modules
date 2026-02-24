#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  sgillis
# Module Author:    sgillis
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

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
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["panel_of_normals"]`
CFG = op.setup_module(
    name = "panel_of_normals",
    version = "1.0",
    subdirectories = ["inputs", "fix_bed", "target_sites", "coverage", "pon_cnn", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _panel_of_normals_input_bam,
    _panel_of_normals_input_capspace,
    _panel_of_normals_canonical_capspace,
    _panel_of_normals_filter_main_chrs,
    _panel_of_normals_output_cnn


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _panel_of_normals_input_bam:
    input:
        bam = ancient(CFG["inputs"]["sample_bam"]),
        bai = ancient(CFG["inputs"]["sample_bai"])
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{capture_space}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)

# Symlinks the capture bed file
rule _panel_of_normals_input_capspace:
    input:
        bed = ancient(reference_files("genomes/{genome_build}/capture_space/{capture_space}.bed"))
    output:
        bed = CFG["dirs"]["inputs"] + "bed/{seq_type}--{genome_build}/{capture_space}.bed"
    run:
        op.absolute_symlink(input.bed, output.bed)

# Removes non-canonical chroms from capture space bed
rule _panel_of_normals_canonical_capspace:
    input:
        bed = str(rules._panel_of_normals_input_capspace.output.bed)
    output:
        bed = CFG["dirs"]["fix_bed"] + "{seq_type}--{genome_build}/{capture_space}.canonical.bed"
    shell:
        op.as_one_line("""
        awk -F"\t" -v OFS="\t" '$1 !~ /(_|M|EBV)/' {input.bed} > {output.bed}
        """)

# Download gene annotation files
rule _panel_of_normals_get_refFlat:
    output:
        refFlat = CFG["dirs"]["inputs"] + "references/{genome_build}/refFlat.final.txt"
    params:
        url = "http://hgdownload.soe.ucsc.edu/goldenPath/",
        build = lambda w: "hg38" if "38" in str({w.genome_build}) else "hg19",
        txt = CFG["dirs"]["inputs"] + "references/{genome_build}/refFlat.txt",
        prefix = lambda w: "chr" if "hg" in str({w.genome_build}) else ""
    threads: 1
    shell:
        op.as_one_line("""
        wget {params.url}{params.build}/database/refFlat.txt.gz -O - | gzip -d
         > {params.txt} &&
         sed 's/chr/{params.prefix}/g' {params.txt} > {output.refFlat}
        """)

rule _panel_of_normals_accessible_regions:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        access = CFG["dirs"]["inputs"] + "references/access.{genome_build}.bed"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["cnvkit"]
    resources:
        **CFG["resources"]["cnvkit"]
    log:
        stdout = CFG["logs"]["inputs"] + "{genome_build}_access.log"
    shell:
        op.as_one_line("""
        cnvkit.py access {input.fasta} -o {output.access} &> {log.stdout}
        """)

# Filters out chrG, chrJ, chrM from bed
rule _panel_of_normals_filter_main_chrs:
    input:
        access = CFG["dirs"]["inputs"] + "references/access.{genome_build}.bed"
    output:
        access_main = CFG["dirs"]["inputs"] + "references/access_main.{genome_build}.bed"
    shell:
        op.as_one_line("""
        grep -v GL {input.access} | grep -v J | grep -v M > {output.access_main}
        """)

# Collects all normals per genome_build--capture_space combo
def _get_normals_per_combo(wildcards):
    CFG = config["lcr-modules"]["panel_of_normals"]
    tbl = CFG["samples"]
    samples = tbl[(tbl.genome_build == wildcards.genome_build) & (tbl.capture_space == wildcards.capture_space)]["sample_id"].tolist()
    normals = expand(
        str(rules._panel_of_normals_input_bam.output.bam),
        sample_id = samples,
        allow_missing = True)
    return normals

# Outputs the sample_ids per capture space combo
rule _panel_of_normals_record_samples:
    input:
        _get_normals_per_combo
    output:
        tsv = CFG["dirs"]["inputs"] + "sample_metadata_tsvs/{seq_type}--{genome_build}/{capture_space}_samples_metadata.tsv"
    params:
        metadata = CFG["samples"]
    threads: 1
    run:
        metadata_for_combo = params.metadata[(params.metadata.genome_build == wildcards.genome_build) & (params.metadata.capture_space == wildcards.capture_space)]
        metadata_for_combo.to_csv(output.tsv, sep="\t", index=False, na_rep='NA')

rule _panel_of_normals_target_sites:
    input:
        access = str(rules._panel_of_normals_filter_main_chrs.output.access_main),
        bam = _get_normals_per_combo,
        bed = str(rules._panel_of_normals_canonical_capspace.output.bed),
        refFlat = str(rules._panel_of_normals_get_refFlat.output.refFlat)
    output:
        target = CFG["dirs"]["target_sites"] + "{seq_type}--{genome_build}/{capture_space}/target_sites.target.bed",
        antitarget = CFG["dirs"]["target_sites"] + "{seq_type}--{genome_build}/{capture_space}/target_sites.antitarget.bed"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["cnvkit"]
    resources:
        **CFG["resources"]["cnvkit"]
    log:
        stdout = CFG["logs"]["target_sites"] + "{seq_type}--{genome_build}/{capture_space}_autobin.log"
    shell:
        op.as_one_line("""
        cnvkit.py autobin {input.bam} -t {input.bed} -g {input.access} --annotate
         {input.refFlat} --short-names --target-output-bed {output.target}
         --antitarget-output-bed {output.antitarget} &> {log.stdout}
        """)

# Get coverage for each sample
rule _panel_of_normals_coverage_target:
    input:
        bam = str(rules._panel_of_normals_input_bam.output.bam),
        bed = str(rules._panel_of_normals_target_sites.output.target),
    output:
        cov = CFG["dirs"]["coverage"] + "target/{seq_type}--{genome_build}/{capture_space}/{sample_id}.targetcoverage.cnn"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["cnvkit"]
    resources:
        **CFG["resources"]["cnvkit"]
    log:
        stdout = CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}_target.log"
    shell:
        op.as_one_line("""
        cnvkit.py coverage {input.bam} {input.bed} -o {output.cov} -p {threads}
         &> {log.stdout}
        """)

# Get coverage of anti-target sample
rule _panel_of_normals_coverage_antitarget:
    input:
        bam = str(rules._panel_of_normals_input_bam.output.bam),
        bed = str(rules._panel_of_normals_target_sites.output.antitarget),
    output:
        cov = CFG["dirs"]["coverage"] + "antitarget/{seq_type}--{genome_build}/{capture_space}/{sample_id}.antitargetcoverage.cnn"
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["cnvkit"]
    resources:
        **CFG["resources"]["cnvkit"]
    log:
        stdout = CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}_antitarget.log"
    shell:
        op.as_one_line("""
        cnvkit.py coverage {input.bam} {input.bed} -o {output.cov} -p {threads}
         &> {log.stdout}
        """)

# Create reference coverage cnn file using all normal cnn files per combo
def _get_coverage_per_combo(wildcards):
    CFG = config["lcr-modules"]["panel_of_normals"]
    tbl = CFG["samples"]
    samples = tbl[(tbl.genome_build == wildcards.genome_build) & (tbl.capture_space == wildcards.capture_space)]["sample_id"].tolist()
    target_cov = expand(
        str(rules._panel_of_normals_coverage_target.output.cov),
        sample_id = samples,
        allow_missing = True)
    antitarget_cov = expand(
        str(rules._panel_of_normals_coverage_antitarget.output.cov),
        sample_id = samples,
        allow_missing = True)
    return{
        "control_target": target_cov,
        "control_antitarget": antitarget_cov
    }

rule _panel_of_normals_create_pon_ref:
    input:
        unpack(_get_coverage_per_combo),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        pon = CFG["dirs"]["pon_cnn"] + "cnn/{seq_type}--{genome_build}/{capture_space}_normal_reference.cnn"
    params:
        male_reference = CFG["options"]["male_ref"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["cnvkit"]
    resources:
        **CFG["resources"]["cnvkit"]
    log:
        stdout = CFG["logs"]["pon_cnn"] + "{seq_type}--{genome_build}/{capture_space}/pon_cnn.log"
    shell:
        op.as_one_line("""
        cnvkit.py reference {input.control_target} {input.control_antitarget}
         --fasta {input.fasta} -o {output.pon} {params.male_reference}
         &> {log.stdout}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _panel_of_normals_output_cnn:
    input:
        cnn = str(rules._panel_of_normals_create_pon_ref.output.pon)
    output:
        cnn = CFG["dirs"]["outputs"] + "cnn/{seq_type}--{genome_build}/{capture_space}_normal_reference.cnn"
    run:
        op.relative_symlink(input.cnn, output.cnn, in_module= True)

rule _panel_of_normals_output_tsv:
    input:
        tsv = str(rules._panel_of_normals_record_samples.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{capture_space}_samples_metadata.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)

# Generates the target sentinels for each run, which generate the symlinks
rule _panel_of_normals_all:
    input:
        expand(
            [
                str(rules._panel_of_normals_output_cnn.output.cnn),
                str(rules._panel_of_normals_output_tsv.output.tsv)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            capture_space = CFG["samples"]["capture_space"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
