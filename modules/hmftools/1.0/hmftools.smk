#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["hmftools"]`
CFG = op.setup_module(
    name = "hmftools",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "prepare_strelka", "amber", "cobalt", "purple", "linx", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _hmftools_input,
    _hmftools_input_references,
    _hmftools_linx_prepare_ensembl,
    _hmftools_output_vcf,
    _hmftools_all,

wildcard_constraints: 
    genome_build = "|".join(CFG["switches"]["het_snps"].keys())

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _hmftools_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"], 
        bai = CFG["inputs"]["sample_bai"], 
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam", 
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai", 
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)

rule _hmftools_input_strelka: 
    input: 
        strelka_vcf = CFG["inputs"]["strelka_vcf"]
    output: 
        strelka_vcf = CFG["dirs"]["inputs"] + "strelka_vcf/{tumour_id}--{normal_id}--matched/{var_type}.passed.vcf"
    run: 
        op.relative_symlink(input.strelka_vcf, output.strelka_vcf)

rule _hmftools_input_gridss: 
    input: 
        gridss_somatic_vcf = CFG["inputs"]["gridss_somatic"], 
        gridss_filtered_vcf = CFG["inputs"]["gridss_somatic_filtered"]
    output: 
        gridss_somatic_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz",
        gridss_filtered_vcf = CFG["dirs"]["inputs"] + "gridss_vcf/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filtered.vcf.gz"
    run: 
        op.relative_symlink(input.gridss_vcf, output.gridss_vcf)
        op.relative_symlink(input.gridss_full_vcf, output.gridss_full_vcf)


rule _hmftools_input_references: 
    input: 
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
    output: 
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.fa"
    conda: 
        CFG["envs"]["samtools"]
    shell: 
        op.as_one_line("""
        ln -s {input.genome_fa} {output.genome_fa} &&
        ln -s {input.genome_fa}.fai {output.genome_fa}.fai &&
        samtools dict {output.genome_fa} -o {output.genome_fa}.dict
        """)

# Prepare Strelka VCF files for use with PURPLE
rule _hmftools_merge_strelka: 
    input: 
        snv_vcf = expand(rules._hmftools_input_strelka.output.strelka_vcf, var_type = "somatic.snvs"), 
        indel_vcf = expand(rules._hmftools_input_strelka.output.strelka_vcf, var_type = "somatic.indels") 
    output: 
        vcf = CFG["dirs"]["prepare_strelka"] + "genome--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.merged.vcf"
    log: CFG["logs"]["prepare_strelka"] + "genome--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merge_strelka.log"
    conda: 
        CFG["envs"]["bcftools"]
    threads: 
        CFG["threads"]["merge_strelka"]
    resources: 
        mem_mb = CFG["mem_mb"]["merge_strelka"]
    shell: 
        op.as_one_line("""
        bcftools concat -a -O v {input.snv_vcf} {input.indel_vcf} | 
        bcftools view -f PASS -O v | 
        bcftools sort - | 
        sed "s/TUMOR/{wildcards.tumour_id}/g" | sed "s/NORMAL/{wildcards.normal_id}/g" 
        > {output.vcf} 2>{log}
        """)

rule _hmftools_annotate_strelka: 
    input: 
        vcf = rules._hmftools_merge_strelka.output.vcf
    output: 
        vcf = vcf = CFG["dirs"]["prepare_strelka"] + "genome--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.merged.annotated.vcf"
    params: 
        purple_jar = "$(dirname $(readlink -e $(which PURPLE)))/purple.jar"
    log: CFG["logs"]["prepare_strelka"] + "genome--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/annotate_strelka.log"
    conda: CFG["envs"]["purple"]
    threads: 
        CFG["threads"]["annotate_strelka"]
    resources: 
        mem_mb = CFG["mem_mb"]["annotate_strelka"]
    shell: 
        op.as_one_line("""
        java -Xmx{resources.mem_mb}m 
            -cp {params.purple_jar} com.hartwig.hmftools.purple.tools.AnnotateStrelkaWithAllelicDepth
            -in {input.vcf} -out {output.vcf}
        """)

# Run AMBER to calculate BAFs
rule _hmftools_amber_paired: 
    input: 
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
    output: 
        vcf = CFG["dirs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched/{tumour_id}.amber.baf.vcf.gz"
    params:
        snps = op.switch_on_wildcard("genome_build", CFG["switches"]["het_snps"]), 
        options = CFG["options"]["amber"]
    log: CFG["logs"]["amber"] + "genome--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/amber.log"
    conda: 
        CFG["envs"]["amber"]
    threads: 
        CFG["threads"]["amber"]
    resources: 
        mem_mb = CFG["mem_mb"]["amber"]
    shell: 
        op.as_one_line("""
        AMBER 
        -reference {wildcards.normal_id} -reference_bam {input.normal_bam}
        -tumor {wildcards.tumour_id} -tumor_bam {input.tumour_bam}
        -output_dir `dirname {output.vcf}`
        -threads {params.threads}
        -loci {params.snps}
        {params.options}
        2>&1 | tee -a {log} 
        """)
    
# rule _hmftools_amber_unpaired: 
#     input: 
#         tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam"
#     output: 
#         out_dir = directory(CFG["dirs"]["amber"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--unmatched")
#     params:
#         snps = op.switch_on_wildcard("genome_build", CFG["switches"]["het_snps"]), 
#         options = CFG["options"]["amber"]
#     log: CFG["logs"]["amber"] + "genome--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/amber.log"
#     conda: 
#         CFG["envs"]["amber"]
#     threads: 
#         CFG["threads"]["amber"]
#     resources: 
#         mem_mb = CFG["mem_mb"]["amber"]
#     shell: 
#         op.as_one_line("""
#         AMBER 
#         -tumor {wildcards.tumour_id} -tumor_bam {input.tumour_bam}
#         -output_dir {output.out_dir}
#         -threads {params.threads}
#         -loci {params.snps}
#         -tumor_only
#         {params.options}
#         2>&1 | tee -a {log} 
#         """)

# Run COBALT to estimate depth across the genome
rule _hmftools_cobalt: 
    input: 
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
    output: 
        ratio = CFG["dirs"]["cobalt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cobalt.ratio.tsv"
    params:
        gc_profile = op.switch_on_wildcard("genome_build", CFG["switches"]["gc_profile"]) 
        options = CFG["options"]["cobalt"]
    conda: 
        CFG["envs"]["cobalt"]
    threads: 
        CFG["threads"]["cobalt"]
    resources: 
        mem_mb = CFG["mem_mb"]["cobalt"]
    shell: 
        op.as_one_line("""
        COBALT 
        -reference {wildcards.normal_id} -reference_bam {input.normal_bam} 
        -tumor {wildcards.tumour_id} -tumor_bam {input.tumour_bam} 
        -output_dir `dirname {output.ratio}` 
        -threads {params.threads} 
        -gc_profile {params.gc_profile} 
        {params.options} 
        2>&1 | tee -a {log} 
        """)

# Run PURPLE for final CNV calling 
rule _hmftools_purple_paired:
    input: 
        amber = rules._hmftools_amber_paired.output.out_dir,
        cobalt = rules._hmftools_cobalt.output.out_dir,
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
        strelka_vcf = rules._hmftools_annotate_strelka.output.vcf, 
        gridss_somatic_vcf = rules._hmftools_input_gridss.output.gridss_somatic_vcf,
        gridss_filtered_vcf = rules._hmftools_input_gridss.output.gridss_filtered_vcf,
        reference_fa = rules._hmftools_input_references.output.genome_fa
    output: 
        files = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched/{tumour_id}.{purple_output}", 
        circos = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched/plot/{tumour_id}.circos.png"
    params: 
        gc_profile = op.switch_on_wildcard("genome_build", CFG["switches"]["gc_profile"]) ,
        options = CFG["options"]["purple"]
        circos = "`which circos`"
    conda: 
        CFG["envs"]["purple"]
    threads: 
        CFG["threads"]["purple"]
    resources: 
        mem_mb = CFG["mem_mb"]["purple"]
    shell: 
        op.as_one_line("""
        PURPLE \
            -reference {wildcards.normal_id} 
            -tumor {wildcards.tumour_id} 
            -output_dir `dirname {output.sv_vcf}` 
            -amber {input.amber} 
            -cobalt {input.cobalt} 
            -gc_profile {params.gc_profile} 
            -ref_genome {input.reference_fa} 
            -somatic_vcf {input.strelka_vcf} 
            -structural_vcf {input.gridss_filtered_vcf} 
            -sv_recovery_vcf {input.gridss_somatic_vcf} 
            -circos {params.circos} 
            {params.options}
        """)
    
# rule _hmftools_purple_unpaired:
#     input: 
#         amber = rules._hmftools_amber_paired.output.out_dir,
#         cobalt = rules._hmftools_cobalt.output.out_dir,
#         tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
#         normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam", 
#         gridss_vcf = rules._hmftools_input_gridss.output.gridss_vcf,
#         gridss_full_vcf = rules._hmftools_input_gridss.output.gridss_full_vcf,
#         reference_fa = rules._hmftools_input_references.output.genome_fa
#     output: 
#         files = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--unmatched/{tumour_id}.{purple_output}", 
#         plots = directory(CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--unmatched/plot") 
#     params: 
#         gc_profile = op.switch_on_wildcard("genome_build", CFG["switches"]["gc_profile"]),
#         options = CFG["options"]["purple"]
#         circos = "`which circos`"
#     conda: 
#         CFG["envs"]["purple"]
#     threads: 
#         CFG["threads"]["purple"]
#     resources: 
#         mem_mb = CFG["mem_mb"]["purple"]
#     shell: 
#         op.as_one_line("""
#         PURPLE 
#             -reference {wildcards.normal_id} 
#             -tumor {wildcards.tumour_id} 
#             -output_dir `dirname {output.sv_vcf}` 
#             -amber {input.amber} 
#             -cobalt {input.cobalt} 
#             -gc_profile {params.gc_profile} 
#             -ref_genome {input.reference_fa} 
#             -structural_vcf {input.gridss_vcf} 
#             -sv_recovery_vcf {input.gridss_full_vcf} 
#             -circos {params.circos} 
#             -tumor_only
#             {params.options}
#         """)


# Prepare LINX Ensembl data cache

rule _hmftools_linx_prepare_ensembl: 
    output: 
        ensembl_cache = directory(CFG["dirs"]["inputs"] + "references/ensembl/HG{ensembl_build}")
    params: 
        url = op.switch_on_wildcard("ensembl_build", CFG["switches"]["ensembl_url"])
        linx_jar = "$(dirname $(readlink -e $(which linx)))/sv-linx.jar"
    conda: 
        CFG["envs"]["linx"]
    shell: 
        op.as_one_line("""
        java -cp {params.jar} com.hartwig.hmftools.linx.gene.GenerateEnsemblDataCache
            -ensembldb {params.url} -ensembl_user "anonymous" -ensembl_pass ""
            -output_dir {output.ensembl_cache} -ref_genome_version {wildcards.ensembl_build}
        """)
    
 def _hmftools_get_ensembl_build(wildcards): 
    CFG = config["lcr-modules"]["hmftools"]
    if wildcards.genome_build == "hg38":
        ensembl_build = "38"
    else:
        ensembl_build = "37"
    ensembl_cache = expand(str(rules._hmftools_linx_prepare_ensembl.output.ensembl_cache), ensembl_build = ensembl_build)
    return ensembl_cache



# Run LINX to cluster and visualize CNV and SV data
rule _hmftools_linx: 
    input: 
        purple_vcf = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.purple.sv.vcf.gz", 
        ensembl_cache = _hmftools_get_ensembl_build
    output: 
        linx_dir = directory(CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}")
    log: CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/linx.log"
    params: 
      linx_db = CFG["references"]["linx_db"], 
      options = CFG["options"]["linx"]
    conda: 
        CFG["envs"]["linx"]
    threads: 
        CFG["threads"]["linx"]
    resources: 
        mem_mb = CFG["mem_mb"]["linx"]
    shell: 
        op.as_one_line("""
        linx 
            -sample {wildcards.tumour_id} 
            -ref_genome_version {wildcards.ensembl_build} 
            -sv_vcf {input.purple_vcf} 
            -purple_dir `dirname {input.purple_vcf}` 
            -output_dir {output.linx_dir} 
            -gene_transcripts_dir {input.ensembl_cache} 
            -fragile_site_file {params.linx_db}/fragile_sites_hmf.csv 
            -line_element_file {params.linx_db}/line_elements.csv 
            -replication_origins_file {params.linx_db}/heli_rep_origins.bed 
            -viral_hosts_file {params.linx_db}/viral_host_ref.csv 
            -check_fusions 
            -fusion_pairs_csv {params.linx_db}/knownFusionPairs.csv 
            -promiscuous_five_csv {params.linx_db}/knownPromiscuousFive.csv 
            -promiscuous_three_csv {params.linx_db}/knownPromiscuousThree.csv  
            -check_drivers 
            -write_vis_data 
            {params.options}
        """)

checkpoint _hmftools_linx_viz: 
    input:
        linx_dir = rules._hmftools_linx.output.linx_dir,
        ensembl_cache = _hmftools_get_ensembl_build
    output:
        plots = directory(CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot"),
        data = directory(CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/data")
    params: 
        linx_jar = "$(ls $(dirname $(readlink -e $(which linx)))/*.jar)", 
        circos = "$(which circos)", 
        options = CFG["options"]["linx_viz"]
    conda:
        CFG["envs"]["linx"]
    threads: 
        CFG["threads"]["linx_viz"]
    resources:
        CFG["mem_mb"]["linx_viz"]
    shell:
        op.as_one_line("""
        java -cp {params.linx_jar} com.hartwig.hmftools.linx.visualiser.SvVisualiser 
            -sample {wildcards.tumour_bam} 
            -gene_transcripts_dir {input.ensembl_cache} 
            -plot_out {output.plots} 
            -data_out {output.data} 
            -segment {input.linx_dir}/{wildcards.tumour_id}.linx.vis_segments.tsv 
            -link {input.linx_dir}/{wildcards.tumour_id}.linx.vis_sv_data.tsv 
            -exon {input.linx_dir}/{wildcards.tumour_id}.linx.vis_gene_exon.tsv 
            -cna {input.linx_dir}/{wildcards.tumour_id}.linx.vis_copy_number.tsv 
            -protein_domain {input.linx_dir}/{wildcards.tumour_id}.linx.vis_protein_domain.tsv 
            -fusion {input.linx_dir}/{wildcards.tumour_id}.linx.fusions.tsv 
            -circos {params.circos} 
            -threads {threads} 
        """) 

# Symlinks the final output files into the module results directory (under '99-outputs/')

rule _hmftools_purple_output:
    input:
        files = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.{purple_output}"
    output:
        files = CFG["dirs"]["outputs"] + "purple_output/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{purple_output}"
    run:
        op.relative_symlink(input.files, output.files)

def _hmftools_predict_purple_output(wildcards): 

    # Define standard purple output file extensions
    purple_files = 
        [
            "purple.purity.tsv", 
            "purple.purity.range.tsv", 
            "purple.cnv.somatic.tsv", 
            "purple.cnv.gene.tsv", 
            "driver.catalog.tsv", 
            "purple.sv.vcf.gz"
        ]

    # Append the somatic vcf output if the tumour is paired
    if wildcards.pair_status == "matched":
        purple_files = purple_files.append("purple.somatic.vcf.gz") 

    # Expand the outputs
    purple_outputs = expand(
        CFG["dirs"]["outputs"] + "purple_output/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{purple_output}",
        purple_output = purple_files,
        **wildcards
    )

    return purple_outputs

rule _hmftools_purple_plots: 
    input:
        plots = CFG["dirs"]["purple"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot/{tumour_id}.{plot_name}"
    output: 
        plots = CFG["dirs"]["outputs"] + "purple_plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{plot_name}"
    run: 
        op.relative_symlink(input.plots, output.plots)

def _hmftools_predict_purple_plots(wildcards): 

    # Define the extensions of all purple plots
    purple_plot_names = [
        "circos.png",
        "copynumber.png",
        "input.png",
        "map.png",
        "purity.range.png",
        "segment.png",
        "somatic.clonality.png",
        "somatic.png",
        "somatic.rainfall.png"
    ]
    
    # Expand the output files
    purple_plots = expand(
        CFG["dirs"]["outputs"] + "purple_plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{plot_name}", 
        plot_name = purple_plot_names, 
        **wildcards
    )

    return purple_plots

rule _hmftools_linx_plots: 
    input:
        plots = CFG["dirs"]["linx"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/plot/{tumour_id}.{plot_name}"
    output: 
        plots = CFG["dirs"]["outputs"] + "linx_plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{plot_name}"
    run: 
        op.relative_symlink(input.plots, output.plots)

def _hmftools_predict_linx_plots(wildcards):
    checkpoint_output = checkpoints.clustering.get(**wildcards).output.plots
    linx_plots = expand(
        CFG["dirs"]["outputs"] + "linx_plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{plot_name}", 
        plot_name = glob_wildcards(os.path.join(checkpoint_output, "{plot_name}.png").plot_name), 
        **wildcards
    )
    return linx_plots


rule _hmftools_dispatch: 
    input: 
        _hmftools_predict_purple_output,
        _hmftools_predict_purple_plots,  
        _hmftools_predict_linx_plots
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
