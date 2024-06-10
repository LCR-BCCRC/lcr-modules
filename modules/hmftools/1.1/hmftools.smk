#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["hmftools"]`
CFG = op.setup_module(
    name = "hmftools",
    version = "1.1",
    subdirectories = ["inputs", "prepare_slms3", "amber", "cobalt", "purple", "linx", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _hmftools_input_bam,
    _hmftools_input_slms3,
    _hmftools_slms3_sample_names, 
    _hmftools_input_gridss, 
    _hmftools_input_references,
    _hmftools_get_cobalt_gc,
    _hmftools_get_cobalt_bed, 
    _hmftools_get_amber_snps, 
    _hmftools_get_purple_drivers, 
    _hmftools_get_linx_db,
    _hmftools_get_ensembl_cache, 
    _hmftools_purple_output,
    _hmftools_purple_plots,  
    _hmftools_all


VERSION_MAP_HMFTOOLS = CFG["options"]["version_map"]

possible_genome_builds = ", ".join(list(VERSION_MAP_HMFTOOLS.keys()))
for genome_build in CFG["runs"]["tumour_genome_build"]:
    assert genome_build in possible_genome_builds, (
        f"Samples table includes genome builds not yet compatible with this module. "
        f"This module is currently only compatible with {possible_genome_builds}. "
    )


masked_string = "" 
if CFG["options"]["use_masked_ref"]:
    masked_string = "_masked"


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _hmftools_input_bam:
    input:
        bam = ancient(CFG["inputs"]["sample_bam"]), 
        bai = ancient(CFG["inputs"]["sample_bai"]), 
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam", 
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai", 
    group: "input_and_vcf"
    wildcard_constraints:
        genome_build = "|".join(VERSION_MAP_HMFTOOLS.keys()),
        pair_status = "matched|unmatched"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)

rule _hmftools_input_slms3: 
    input: 
        vcf = CFG["inputs"]["slms3_vcf"], 
    output: 
        vcf = CFG["dirs"]["inputs"] + "slms3_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/slms3.vcf.gz" 
    group: "input_and_vcf"
    run: 
        op.relative_symlink(input.vcf, output.vcf)

rule _hmftools_input_gridss: 
    input: 
        gridss_somatic_vcf = CFG["inputs"]["gridss_somatic"], 
        gridss_somatic_tbi = CFG["inputs"]["gridss_somatic_tbi"], 
        gridss_filtered_vcf = CFG["inputs"]["gridss_somatic_filtered"], 
        gridss_filtered_tbi = CFG["inputs"]["gridss_somatic_filtered_tbi"]
    output: 
        gridss_somatic_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz",
        gridss_somatic_tbi = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz.tbi",
        gridss_filtered_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.vcf.gz", 
        gridss_filtered_tbi = CFG["dirs"]["inputs"] + "gridss_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.vcf.gz.tbi"
    group: "input_and_vcf"
    run: 
        op.absolute_symlink(input.gridss_somatic_vcf, output.gridss_somatic_vcf)
        op.absolute_symlink(input.gridss_somatic_tbi, output.gridss_somatic_tbi)
        op.absolute_symlink(input.gridss_filtered_vcf, output.gridss_filtered_vcf)
        op.absolute_symlink(input.gridss_filtered_tbi, output.gridss_filtered_tbi)

# Rules to download and setup reference files

rule _hmftools_input_references: 
    input: 
        genome_fa = reference_files("genomes/{genome_build}" + masked_string + "/genome_fasta/genome.fa"),
        genome_fai = reference_files("genomes/{genome_build}" + masked_string + "/genome_fasta/genome.fa.fai"),
        genome_dict = reference_files("genomes/{genome_build}" + masked_string + "/genome_fasta/genome.dict")
    output: 
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}" + masked_string + "/genome_fa/genome.fa", 
        genome_fai = CFG["dirs"]["inputs"] + "references/{genome_build}" + masked_string + "/genome_fa/genome.fa.fai", 
        genome_dict = CFG["dirs"]["inputs"] + "references/{genome_build}" + masked_string + "/genome_fa/genome.dict"
    shell: 
        op.as_one_line("""
        ln -s {input.genome_fa} {output.genome_fa} &&
        ln -s {input.genome_fai} {output.genome_fai} &&
        ln -s {input.genome_dict} {output.genome_dict}
        """)

rule _hmftools_get_cobalt_gc: 
    output: 
        gc = CFG["dirs"]["inputs"] + "references/{genome_build}/cobalt/GC_profile.1000bp.cnp"
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/cobalt",
        alt_build = lambda w: VERSION_MAP_HMFTOOLS[w.genome_build]
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -O {output.gc} {params.url}/GC_profile.1000bp.{params.alt_build}.cnp'

rule _hmftools_get_cobalt_bed: 
    output: 
        bed = CFG["dirs"]["inputs"] + "references/{genome_build}/cobalt/DiploidRegions.bed"
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/cobalt",
        alt_build = lambda w: VERSION_MAP_HMFTOOLS[w.genome_build]
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -O {output.bed} {params.url}/DiploidRegions.{params.alt_build}.bed'

rule _hmftools_get_amber_snps: 
    output: 
        vcf = CFG["dirs"]["inputs"] + "references/{genome_build}/amber/GermlineHetPon.vcf.gz", 
        snpcheck = CFG["dirs"]["inputs"] + "references/{genome_build}/amber/Amber.snpcheck.vcf.gz"
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/amber",
        alt_build = lambda w: VERSION_MAP_HMFTOOLS[w.genome_build]
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -O {output.vcf} {params.url}/GermlineHetPon.{params.alt_build}.vcf.gz; '
        'wget -O {output.snpcheck} {params.url}/Amber.snpcheck.{params.alt_build}.vcf'

rule _hmftools_get_purple_drivers: 
    output: 
        hotspots = CFG["dirs"]["inputs"] + "references/{genome_build}/purple/KnownHotspots.vcf.gz", 
        gene_panel = CFG["dirs"]["inputs"] + "references/{genome_build}/purple/DriverGenePanel.tsv"
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/purple",
        alt_build = lambda w: VERSION_MAP_HMFTOOLS[w.genome_build]
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -O {output.hotspots} {params.url}/KnownHotspots.somatic.{params.alt_build}.vcf.gz && '
        'wget -O {output.hotspots}.tbi {params.url}/KnownHotspots.somatic.{params.alt_build}.vcf.gz.tbi && '
        'wget -O {output.gene_panel} {params.url}/DriverGenePanel.{params.alt_build}.tsv'

rule _hmftools_get_linx_db: 
    output: 
        directory(CFG["dirs"]["inputs"] + "references/{genome_build}/linx_db")
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/linx/Linx", 
        alt_build = lambda w: VERSION_MAP_HMFTOOLS[w.genome_build]
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -r -np -nd -P {output} -A .bed,.csv {params.url}/{params.alt_build}  && '
        'wget -O {output}/viral_host_ref.csv {params.url}/viral_host_ref.csv'

rule _hmftools_get_ensembl_cache: 
    output: 
        cache = directory(CFG["dirs"]["inputs"] + "references/{genome_build}/ensembl_cache/"), 
        complete = touch(CFG["dirs"]["inputs"] + "references/{genome_build}/ensembl_cache/cache.complete")
    params: 
        url = "www.bcgsc.ca/downloads/morinlab/hmftools-references/ensembl_data_cache",
        alt_build = lambda w: VERSION_MAP_HMFTOOLS[w.genome_build] 
    conda: 
        CFG["conda_envs"]["wget"]
    shell: 
        'wget -O {output.cache}/{params.alt_build}.zip {params.url}/{params.alt_build}.zip && '
        'unzip -d {output.cache} {output.cache}/{params.alt_build}.zip'

# Prepare SLMS-3 VCF files for use with PURPLE
# SnpEff annotation enables driver discovery logic

rule _hmftools_slms3_sample_names: 
    input: 
        vcf = rules._hmftools_input_slms3.output.vcf
    output: 
        vcf = temp(CFG["dirs"]["prepare_slms3"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/tmp.slms3.vcf")
    log: CFG["dirs"]["prepare_slms3"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/vcf_sample_names.log"
    conda: 
        CFG["conda_envs"]["bcftools"]
    threads: CFG["threads"]["vcf_sample_names"]
    resources: 
        **CFG["resources"]["vcf_sample_names"]
    group: "input_and_vcf"
    shell: 
        op.as_one_line("""
        bcftools view -Ov {input.vcf} | 
        sed 's/TUMOR/{wildcards.tumour_id}/g' | 
        sed 's/NORMAL/{wildcards.normal_id}/g'  
        > {output.vcf}
        """)

rule _hmftools_snpeff_vcf: 
    input: 
        vcf = str(rules._hmftools_slms3_sample_names.output.vcf)
    output: 
        sample_key = temp(CFG["dirs"]["prepare_slms3"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sample_key.txt"),
        vcf = temp(CFG["dirs"]["prepare_slms3"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/slms3.snpeff.vcf.gz")
    resources: 
        **CFG["resources"]["snpeff"]
    params: 
        snpeff_build = lambda w: {
            "grch37": "GRCh37.75", 
            "hs37d5": "GRCh37.75", 
            "hg38": "hg38"
        }[w.genome_build], 
        config = "$(readlink -e $(which snpEff)).config",
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    log: 
        CFG["logs"]["prepare_slms3"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/snpeff_slms3.log"
    conda: 
        CFG["conda_envs"]["snpeff"]
    threads: 
        CFG["threads"]["snpeff"]
    shell: 
        op.as_one_line("""
        printf "{wildcards.normal_id}\t{wildcards.tumour_id}\n" > {output.sample_key} && 
        snpEff -Xmx{params.mem_mb}m   
        -c {params.config} -noStats
        -cancer -cancerSamples {output.sample_key} 
        {params.snpeff_build} {input.vcf} | 
        bcftools view -Oz -o {output.vcf} - && 
        bcftools index -t {output.vcf}
        """)


# Run AMBER to calculate BAFs
rule _hmftools_amber_matched: 
    input: 
        tumour_bam = ancient(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam"),
        normal_bam = ancient(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"), 
        snps = str(rules._hmftools_get_amber_snps.output.vcf), 
        fasta = str(rules._hmftools_input_references.output.genome_fa)
    output: 
        vcf = CFG["dirs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.amber.baf.vcf.gz"
    resources: 
        **CFG["resources"]["amber"]
    params:
        options = CFG["options"]["amber"], 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8) 
    log: CFG["logs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/amber.log"
    wildcard_constraints: 
        pair_status = "matched"
    conda: 
        CFG["conda_envs"]["amber"]
    threads: 
        CFG["threads"]["amber"]
    shell: 
        op.as_one_line("""
        AMBER -Xmx{params.jvmheap}m
        -reference {wildcards.normal_id} -reference_bam {input.normal_bam}
        -tumor {wildcards.tumour_id} -tumor_bam {input.tumour_bam}
        -output_dir `dirname {output.vcf}`
        -threads {threads}
        -loci {input.snps}
        -ref_genome {input.fasta} 
        {params.options}
        2>&1 | tee -a {log} 
        """)

rule _hmftools_amber_unmatched: 
    input: 
        tumour_bam = ancient(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam"),
        normal_bam = ancient(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"), 
        snps = str(rules._hmftools_get_amber_snps.output.vcf), 
        fasta = str(rules._hmftools_input_references.output.genome_fa)
    output: 
        vcf = CFG["dirs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.amber.baf.vcf.gz"
    resources: 
        **CFG["resources"]["amber"]
    params:
        options = CFG["options"]["amber"], 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    log: CFG["logs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/amber.log"
    wildcard_constraints: 
        pair_status = "unmatched"
    conda: 
        CFG["conda_envs"]["amber"]
    threads: 
        CFG["threads"]["amber"]
    shell: 
        op.as_one_line("""
        AMBER -Xmx{params.jvmheap}m
        -tumor_only 
        -tumor {wildcards.tumour_id} -tumor_bam {input.tumour_bam} 
        -output_dir `dirname {output.vcf}`
        -threads {threads}
        -loci {input.snps}
        -ref_genome {input.fasta}
        {params.options}
        2>&1 | tee -a {log} 
        """)

# Run COBALT to estimate depth across the genome
rule _hmftools_cobalt: 
    input: 
        tumour_bam = ancient(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam"),
        normal_bam = ancient(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"), 
        gc_profile = str(rules._hmftools_get_cobalt_gc.output.gc), 
        fasta = str(rules._hmftools_input_references.output.genome_fa)
    output: 
        tumour_ratio = CFG["dirs"]["cobalt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cobalt.ratio.pcf", 
        normal_ratio = CFG["dirs"]["cobalt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{normal_id}.cobalt.ratio.pcf", 
        tumour_tsv = temp(CFG["dirs"]["cobalt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cobalt.ratio.tsv"), 
    log: ratio = CFG["logs"]["cobalt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cobalt.log"
    resources: 
        **CFG["resources"]["cobalt"]
    params:
        options = CFG["options"]["cobalt"], 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    wildcard_constraints: 
        pair_status = "matched|unmatched"
    conda: 
        CFG["conda_envs"]["cobalt"]
    threads: 
        CFG["threads"]["cobalt"]
    shell: 
        op.as_one_line("""
        COBALT -Xmx{params.jvmheap}m
        -reference {wildcards.normal_id} -reference_bam {input.normal_bam} 
        -tumor {wildcards.tumour_id} -tumor_bam {input.tumour_bam} 
        -ref_genome {input.fasta} 
        -output_dir `dirname {output.tumour_ratio}` 
        -threads {threads} 
        -gc_profile {input.gc_profile} 
        {params.options} 
        2>&1 | tee -a {log} 
        """)


# Run PURPLE for final CNV calling 

# Define variables for output file names
purple_out = [
    "purity.tsv", 
    "purity.range.tsv", 
    "cnv.gene.tsv", 
    "sv.vcf.gz", 
]
purple_plots = [
        "circos",
        "input",
        "map",
        "purity.range",
        "segment"
    ]

rule _hmftools_purple_matched:
    input: 
        amber = CFG["dirs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.amber.baf.vcf.gz",
        cobalt_tumour = str(rules._hmftools_cobalt.output.tumour_ratio),
        cobalt_normal = str(rules._hmftools_cobalt.output.normal_ratio),
        cobalt_tumour_tsv = str(rules._hmftools_cobalt.output.tumour_tsv), 
        slms3_vcf = str(rules._hmftools_snpeff_vcf.output.vcf), 
        gridss_somatic_vcf = str(rules._hmftools_input_gridss.output.gridss_somatic_vcf),
        gridss_filtered_vcf = str(rules._hmftools_input_gridss.output.gridss_filtered_vcf),
        reference_fa = str(rules._hmftools_input_references.output.genome_fa), 
        gene_panel = str(rules._hmftools_get_purple_drivers.output.gene_panel), 
        hotspots = str(rules._hmftools_get_purple_drivers.output.hotspots), 
        gc_profile = str(rules._hmftools_get_cobalt_gc.output.gc)
    output: 
        files = expand(CFG["dirs"]["purple"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/{{tumour_id}}.purple.{out_file}", 
            out_file = purple_out), 
        plots = expand(CFG["dirs"]["purple"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/plot/{{tumour_id}}.{plot_name}.png", 
            plot_name = purple_plots)
    log: CFG["logs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/purple.log"
    resources: 
        **CFG["resources"]["purple"]
    params: 
        outdir = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}", 
        options = CFG["options"]["purple"],
        circos = "`which circos`", 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.9)
    wildcard_constraints: 
        pair_status = "matched|unmatched", 
        out_file = "|".join(purple_out), 
        plot_name = "|".join(purple_plots)
    conda: 
        CFG["conda_envs"]["purple"]
    threads: 
        CFG["threads"]["purple"]
    shell: 
        op.as_one_line("""
        PURPLE -Xmx{params.jvmheap}m -driver_catalog 
            -reference {wildcards.normal_id} 
            -tumor {wildcards.tumour_id} 
            -output_dir {params.outdir} 
            -amber `dirname {input.amber}` 
            -cobalt `dirname {input.cobalt_tumour}` 
            -gc_profile {input.gc_profile} 
            -ref_genome {input.reference_fa}
            -somatic_hotspots {input.hotspots} 
            -driver_gene_panel {input.gene_panel}  
            -somatic_vcf {input.slms3_vcf} 
            -structural_vcf {input.gridss_filtered_vcf} 
            -sv_recovery_vcf {input.gridss_somatic_vcf} 
            -circos {params.circos} 
            {params.options}
            -threads {threads}
            2>&1 | tee -a {log}
        """)




# Run LINX to cluster and visualize CNV and SV data
rule _hmftools_linx: 
    input: 
        purple_vcf = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.purple.sv.vcf.gz", 
        ensembl_cache = str(rules._hmftools_get_ensembl_cache.output.cache), 
        linx_db = str(rules._hmftools_get_linx_db.output)
    output: 
        clusters = CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.linx.vis_sv_data.tsv", 
        svs = CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.linx.svs.tsv"
    log: CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/linx.log"
    resources: 
        **CFG["resources"]["linx"]
    params: 
      ref_genome_version = lambda w: VERSION_MAP_HMFTOOLS[w.genome_build], 
      jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
      options = CFG["options"]["linx"], 
      cache_subdir = lambda w: config["lcr-modules"]["hmftools"]["dirs"]["inputs"] + "references/" + w.genome_build + "/ensembl_cache/" + VERSION_MAP_HMFTOOLS[w.genome_build]
    conda: 
        CFG["conda_envs"]["linx"]
    threads: 
        CFG["threads"]["linx"]
    shell: 
        op.as_one_line("""
        linx -Xmx{params.jvmheap}m 
            -sample {wildcards.tumour_id} 
            -ref_genome_version {params.ref_genome_version} 
            -sv_vcf {input.purple_vcf} 
            -purple_dir `dirname {input.purple_vcf}` 
            -output_dir `dirname {output.clusters}`  
            -gene_transcripts_dir {params.cache_subdir} 
            -fragile_site_file {input.linx_db}/fragile_sites_hmf.{params.ref_genome_version}.csv 
            -line_element_file {input.linx_db}/line_elements.{params.ref_genome_version}.csv 
            -viral_hosts_file {input.linx_db}/viral_host_ref.csv 
            -known_fusion_file {input.linx_db}/known_fusion_data.{params.ref_genome_version}.csv 
            -check_fusions 
            -check_drivers 
            -write_vis_data 
            {params.options} 
            2>&1 | tee -a {log}
        """)

rule _hmftools_linx_viz: 
    input:
        clusters = rules._hmftools_linx.output.clusters,
        svs = rules._hmftools_linx.output.svs,
        ensembl_cache = str(rules._hmftools_get_ensembl_cache.output.cache)
    output:
        plots = directory(CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot"),
        data = directory(CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/data")
    log: CFG["logs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/linx_viz.log"
    resources:
        **CFG["resources"]["linx_viz"]
    params: 
        linx_jar = "$(ls $(dirname $(readlink -e $(which linx)))/*.jar)", 
        circos = "$(which circos)", 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
        options = CFG["options"]["linx_viz"], 
        cache_subdir = lambda w: config["lcr-modules"]["hmftools"]["dirs"]["inputs"] + "references/" + w.genome_build + "/ensembl_cache/" + VERSION_MAP_HMFTOOLS[w.genome_build], 
        alt_build = lambda w: VERSION_MAP_HMFTOOLS[w.genome_build]
    conda:
        CFG["conda_envs"]["linx"]
    threads: 
        CFG["threads"]["linx_viz"]
    
    shell:
        op.as_one_line("""
        to_plot=$(dirname {input.svs})/to_plot.tsv;
        tail -n +2 {input.svs} | awk '{{FS=OFS="\\t"}} $4 != "" {{print $3}}' | sort | uniq > $to_plot; 
        if [[ $(cat $to_plot | wc -l) -lt 50 ]]; then
            cat $to_plot | while read cluster; do 
                java -Xmx{params.jvmheap}m -cp {params.linx_jar} com.hartwig.hmftools.linx.visualiser.SvVisualiser 
                    -sample {wildcards.tumour_id} 
                    -ref_genome_version V{params.alt_build}
                    -gene_transcripts_dir {params.cache_subdir} 
                    -plot_out {output.plots} 
                    -data_out {output.data} 
                    -vis_file_dir $(dirname {input.clusters})
                    -circos {params.circos} 
                    -threads {threads}  
                    -clusterId $cluster
                    -plot_cluster_genes         
                    2>&1 | tee -a {log};
            done; 
        else 
            echo "Too many clusters to plot for {wildcards.tumour_id}--{wildcards.normal_id}--{wildcards.pair_status}. See chromosome outputs and consider manually selecting clusters to plot. " 2>&1 | tee -a {log}; 
        fi;
        for chrom in $(tail -n +2 {input.clusters} | cut -f8 | sort | uniq); do 
            java -Xmx{params.jvmheap}m -cp {params.linx_jar} com.hartwig.hmftools.linx.visualiser.SvVisualiser 
                -sample {wildcards.tumour_id} 
                -ref_genome_version V{params.alt_build}
                -gene_transcripts_dir {params.cache_subdir} 
                -plot_out {output.plots} 
                -data_out {output.data} 
                -vis_file_dir $(dirname {input.clusters})
                -circos {params.circos} 
                -threads {threads}  
                -chromosome ${{chrom}}
                2>&1 | tee -a {log}; 
        done
        """) 
        



# Symlinks the final output files into the module results directory (under '99-outputs/')

rule _hmftools_purple_output:
    input:
        files = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.purple.{out_file}" 
    output:
        files = CFG["dirs"]["outputs"] + "purple_output/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.purple.{out_file}" 
    wildcard_constraints: 
        out_file = "|".join(purple_out) 
    run:
        op.relative_symlink(input.files, output.files, in_module=True)

rule _hmftools_purple_plots: 
    input:
        plots = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot/{tumour_id}.{plot_name}.png"
    output: 
        plots = CFG["dirs"]["outputs"] + "purple_plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{plot_name}.png"
    wildcard_constraints: 
        plot_name = "|".join(purple_plots)
    run: 
        op.relative_symlink(input.plots, output.plots, in_module=True)


rule _hmftools_linx_plots: 
    input:
        plots = CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot", 
    output: 
        plots = CFG["dirs"]["outputs"] + "linx_plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.symlinked"
    shell: 
        op.as_one_line("""
        workdir=$PWD &&
        cd `dirname {output.plots}` && 
        find $workdir/{input.plots} -type f -name "*.png" -exec cp -s {{}} . \; && 
        touch $workdir/{output.plots} && 
        cd $workdir
        """)

rule _hmftools_dispatch: 
    input: 
        files = expand(CFG["dirs"]["outputs"] + "purple_output/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}.purple.{out_file}", 
            out_file = purple_out), 
        plots = expand(CFG["dirs"]["outputs"] + "purple_plots/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}.{plot_name}.png", 
            plot_name = purple_plots),
        linx = rules._hmftools_linx_plots.output.plots
    output: 
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _hmftools_all:
    input:
        expand(
            [
                str(rules._hmftools_dispatch.output.dispatched),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
