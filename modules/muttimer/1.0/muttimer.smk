#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Shaghayegh Soudi
# Module Author:    Shaghayegh Soudi
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["muttimer"]`
CFG = op.setup_module(
    name = "muttimer",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "processed_inputs", "muttimer", "outputs"],
)


# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _muttimer_input_maf,
    _muttimer_all,


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _muttimer_input_maf:
    input:       
        maf = CFG["inputs"]["maf"] 
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched_final_augmented.maf"  
    run:
        op.absolute_symlink(input.maf, output.maf)  


rule _muttimer_input_battenberg:
    input:      
        cellularity = CFG["inputs"]["cellularity"],
        subclones = CFG["inputs"]["subclones"]
    output:
        cellularity = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}.cellularity_ploidy.txt",
        subclones = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}.subclones.txt"
    run:
        op.absolute_symlink(input.subclones, output.subclones)
        op.absolute_symlink(input.cellularity, output.cellularity)


#dpclust_extension2 =  str(CFG["options"]["dpclust_run"]["iters"]) + "iters_" + str(CFG["options"]["dpclust_run"]["burnin"]) + "burnin"
rule _muttimer_input_dpclust:
    input: 
        clust = CFG["inputs"]["clust"]
    output:
        clust= CFG["dirs"]["inputs"] + "dpclust/{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_DPoutput_2000iters_1000burnin_seed123/{tumour_id}_2000iters_1000burnin_bestClusterInfo.txt" 
    run:
        op.absolute_symlink(input.clust, output.clust)
        

rule _muttimer_formatted_inputs: 
    input:
        maf = str(rules._muttimer_input_maf.output.maf), 
        cnv = str(rules._muttimer_input_battenberg.output.subclones),
        cellularity = str(rules._muttimer_input_battenberg.output.cellularity),   
        clust = str(rules._muttimer_input_dpclust.output.clust),
    output: 
        processed_vcf = CFG["dirs"]["processed_inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_maf_convertedTo.vcf",
        processed_cnv = CFG["dirs"]["processed_inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cnv_mtr.txt",
        processed_clust = CFG["dirs"]["processed_inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cluster_mtr.txt"
    log:
        stderr = CFG["logs"]["processed_inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}.create_processed_mutationtimer_inputs.stderr.log",
        stdout = CFG["logs"]["processed_inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}.create_processed_mutationtimer_inputs.stdout.log"
    params:      
        out_dir = CFG["dirs"]["processed_inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        script = CFG["scripts"]["mtr_make_inputs"]
    conda:
        CFG["conda_envs"]["muttimer"]    
    shell:
        op.as_one_line("""
        Rscript {params.script} -s {wildcards.tumour_id} -m {input.maf} -p {input.cellularity} -c {input.cnv} -l {input.clust} -o {params.out_dir}
        2> {log.stderr} > {log.stdout} 
        """)


rule _install_muttimer:
    output:
        complete = CFG["dirs"]["inputs"] + "mutationtimer_dependencies_installed.success"
    conda:
        CFG["conda_envs"]["muttimer"]
    log:
        input = CFG["logs"]["inputs"] + "input.log"
    shell:
        """
        R -q -e 'BiocManager::install("rtracklayer",force = TRUE)' &&
        R -q -e 'BiocManager::install(c("optparse","GenomicFeatures","VariantAnnotation","GenomicRanges","IRanges"))' &&
        R -q -e 'devtools::install_github("mg14/mg14")' >> {log.input} && 
        R -q -e 'devtools::install_github("gerstung-lab/MutationTimeR")' >> {log.input} &&
        touch {output.complete}"""


rule _muttimer_run: 
    input:
        vcf = str(rules._muttimer_formatted_inputs.output.processed_vcf), 
        cnv = str(rules._muttimer_formatted_inputs.output.processed_cnv),
        clust = str(rules._muttimer_formatted_inputs.output.processed_clust),
        installed = CFG["dirs"]["inputs"] + "mutationtimer_dependencies_installed.success"
    output: 
        clust_summary_mtr = CFG["dirs"]["muttimer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_copy_number_segments_molecular_time.table",
        variant_annotated_mtr = CFG["dirs"]["muttimer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_variants_annotatted.table",
        varianttype_count_mtr = CFG["dirs"]["muttimer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_variant_type_counts.table",
        molecular_time_mtr = CFG["dirs"]["muttimer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_copy_number_segments_molecular_time.png"
    log:
        stderr = CFG["logs"]["muttimer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}.run_mutationtimer_inputs.stderr.log",
        stdout = CFG["logs"]["muttimer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}.run_mutationtimer_inputs.stdout.log" 
    params:      
        script = CFG["scripts"]["mtr_main"],
        main_out_dir = CFG["dirs"]["muttimer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}"
    conda:
        CFG["conda_envs"]["muttimer"]    
    shell:
        op.as_one_line("""
        Rscript {params.script} -s {wildcards.tumour_id} -v {input.vcf} -c {input.cnv} -l {input.clust} -o {params.main_out_dir}
        2> {log.stderr} > {log.stdout}
        """)
        
# Symlinks the final output files into the module results directory (und er '99-outputs/')
# All plots generated by Battenberg are symlinked using a glob for convenience
rule _muttimer_output:
    input:
        clust_summary = rules._muttimer_run.output.clust_summary_mtr,
        variant_annotated = rules._muttimer_run.output.variant_annotated_mtr,
        varianttype_count = rules._muttimer_run.output.varianttype_count_mtr,
        molecular_time = rules._muttimer_run.output.molecular_time_mtr
    output:
        clust_summary = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}_copy_number_segments_molecular_time.table",
        variant_annotated = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}_variants_annotatted.table", 
        varianttype_count = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}_variant_type_counts.table",
        molecular_time = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}_copy_number_segments_molecular_time.png",
    run:
        op.relative_symlink(input.clust_summary, output.clust_summary,in_module=True)
        op.relative_symlink(input.variant_annotated, output.variant_annotated,in_module=True)
        op.relative_symlink(input.varianttype_count, output.varianttype_count,in_module=True)
        op.relative_symlink(input.molecular_time, output.molecular_time,in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _muttimer_all:
    input:
        expand(
            [
                rules._muttimer_output.output.clust_summary,
                rules._muttimer_output.output.variant_annotated,
                rules._muttimer_output.output.varianttype_count,
                rules._muttimer_output.output.molecular_time
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
