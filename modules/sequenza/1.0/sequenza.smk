#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["sequenza"]`
CFG = op.setup_module(
    name = "sequenza",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "sequenza","seqz_filtered", "outputs"]
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _sequenza_input_bam,
#    _sequenza_output_seg,
    _sequenza_all,


##### RULES #####
#print("Sequenza CFG")

#print(CFG)
print("Sequenza CFG inputs")
print(CFG["inputs"])
# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _sequenza_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)


rule _sequenza_bam2seqz:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        sequenza_gc = reference_files("genomes/{genome_build}/annotations/{genome_build}.gc50Base.txt.gz"),
        genome = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        seqz = CFG["dirs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_{chr}.binned.out.seqz.gz"
    threads:
        CFG["threads"]["bam2seqz"]
    log:
        stdout = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_bam2seqz.stdout_{chr}.log",
        stderr = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_bam2seqz.stderr_{chr}.log"
    resources:
        mem_mb = CFG["mem_mb"]["bam2seqz"]
    conda:
        CFG["conda_envs"]["sequenza"]
    shell:
        "sequenza-utils bam2seqz --qlimit 30 -gc {input.sequenza_gc} --fasta {input.genome} "
        "-n {input.normal_bam} -t {input.tumour_bam} --chromosome {wildcards.chr} | sequenza-utils "
        "seqz_binning -w 300 -s - | gzip > {output}"

merged_string = CFG["dirs"]["sequenza"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/{{tumour_id}}_{chr}.binned.out.seqz.gz"

rule _sequenza_merge_seqz:
    input:
        expand(merged_string,chr=CFG["chroms"])
    output:
        merged_seqz = CFG["dirs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.merged.seqz.gz"
    params:
        merge = CFG["options"]["merge_seqz"]
    log:
        stdout = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_bam2seqz.stdout.log",
        stderr = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_bam2seqz.stderr.log"
    resources: mem_mb = 10000
    shell:
        "bash {params.merge} {input} | gzip > {output.merged_seqz}"


rule _sequenza_filter_seqz:
    input:
        merged_seqz = rules._sequenza_merge_seqz.output.merged_seqz
    output:
        filtered_seqz = CFG["dirs"]["seqz_filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.filtered.seqz.gz"
    log:
        stdout = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_filter_seqz.stdout.log",
        stderr = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_filter_seqz.stderr.log"
    resources: mem_mb = 20000
    params:
        dbsnp_pos = "reference/genomes/{genome_build}/annotations/{genome_build}.dbsnp.pos.sort.C", #fix this to use the path to the reference files based on the config
        filter_seqz = CFG["options"]["filter_seqz"]
    shell:
        "{params.filter_seqz} {input.merged_seqz} {params.dbsnp_pos} | gzip > {output.filtered_seqz}"

rule _sequenza_unfiltered_analysis:
    input:
        merged_seqz = rules._sequenza_merge_seqz.output.merged_seqz,
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"
    output:
        segments = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_unfiltered_segments.txt"
    log:
        stdout = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_unfiltered_analysis.stdout.log",
        stderr = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_unfiltered_analysis.stderr.log"
    resources: mem_mb = 100000
    conda:
        CFG["conda_envs"]["sequenza"]
    threads:
        CFG["threads"]["sequenza"]
    params:
        sequenza_analysis = CFG["options"]["sequenza_analysis"],
        outdir = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}",
        sample_id = "{tumour_id}_unfiltered",
        calc_sex_status = CFG["options"]["calc_sex_status"]
    shell:
    	"""mkdir -p {params.outdir} && sex=$({params.calc_sex_status} {input.normal_bam} a | tail -1 | awk '{{if( $3/$2 > 0.1) print "male"; else print "female"}}') ; """
        "Rscript {params.sequenza_analysis} {input.merged_seqz} {params.outdir} $sex {threads} {params.sample_id} {wildcards.genome_build}"

rule _sequenza_filtered_analysis:
    input:
        filtered_seqz = rules._sequenza_filter_seqz.output.filtered_seqz,
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"
    output:
        segments = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_filtered_segments.txt"
    log:
        stdout = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_filtered_analysis.stdout.log",
        stderr = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_sequenza_filtered_analysis.stderr.log"
    resources: mem_mb = 100000
    conda:
        CFG["conda_envs"]["sequenza"]
    threads:
        CFG["threads"]["sequenza"]
    params:
        sequenza_analysis = CFG["options"]["sequenza_analysis"],
        outdir = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}",
        sample_id = "{tumour_id}_filtered",
        calc_sex_status = CFG["options"]["calc_sex_status"]
    shell:
    	"""mkdir -p {params.outdir} && sex=$({params.calc_sex_status} {input.normal_bam} a | tail -1 | awk '{{if( $3/$2 > 0.1) print "male"; else print "female"}}') ; """
        "Rscript {params.sequenza_analysis} {input.filtered_seqz} {params.outdir} $sex {threads} {params.sample_id} {wildcards.genome_build}"	

#I am not fond of having to duplicate these rules but could not figure out how to expand the targets for another wildcard (e.g. filter-mode) so this is an inelegant workaround
rule _sequenza_unfiltered_igv_segments:
    input:
        rules._sequenza_unfiltered_analysis.output.segments
    output:
        igv = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_unfiltered_segments.igv.seg"
    log:
        stderr = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_unfiltered_cnv2igv.stderr.log"
    conda:
        CFG["conda_envs"]["sequenza"]
    params:
        cnv2igv =  CFG["options"]["cnv2igv"],
        sample_id = "{tumour_id}"
    shell:
        "python "
        "{params.cnv2igv} --mode sequenza --sample {params.sample_id} "
        "{input} > {output} 2> {log.stderr}"

rule _sequenza_filtered_igv_segments:
    input:
        rules._sequenza_filtered_analysis.output.segments
    output:
        igv = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_filtered_segments.igv.seg"
    log:
        stderr = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_filtered_cnv2igv.stderr.log"
    conda:
        CFG["conda_envs"]["sequenza"]
    params:
        cnv2igv =  CFG["options"]["cnv2igv"],
        sample_id = "{tumour_id}"
    shell:
        "python "
        "{params.cnv2igv} --mode sequenza --sample {params.sample_id} "
        "{input} > {output} 2> {log.stderr}"

# Generates the target sentinels for each run, which generate the symlinks
rule _sequenza_all:
    input:
        expand(
            [
                rules._sequenza_filter_seqz.output.filtered_seqz,
                rules._sequenza_merge_seqz.output.merged_seqz,
                rules._sequenza_filtered_analysis.output.segments,
                rules._sequenza_unfiltered_analysis.output.segments,
                rules._sequenza_filtered_igv_segments.output.igv,
                rules._sequenza_unfiltered_igv_segments.output.igv
                # TODO: If applicable, add other output rules here
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
