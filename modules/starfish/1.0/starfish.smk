#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["starfish"]`
CFG = op.setup_module(
    name = "starfish",
    version = "1.0",
    subdirectories = ["inputs", "starfish", "vcf_to_bed", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _starfish_input_vcf,
    _starfish_output_vcf,
    _starfish_output_bed,
    _starfish_all,


##### GLOBAL VARIABLES TO HELP READABILITY IN RULES #####

tool1 = CFG["inputs"]["sample_tool"][0]
tool2 = CFG["inputs"]["sample_tool"][1]
output_base_vcf = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/vcf/{tumour_id}--{normal_id}--{pair_status}."
output_base_bed = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/bed/{tumour_id}--{normal_id}--{pair_status}."
run_starfish_base = CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/"
vcf_to_bed_base = CFG["dirs"]["vcf_to_bed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/"

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _starfish_input_vcf:
    input:
        vcf1 = CFG["inputs"]["sample_vcf"][0],
        vcf2 = CFG["inputs"]["sample_vcf"][1]
    output:
        vcf1 = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}." + tool1 + ".vcf.gz",
        vcf2 = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}." + tool2 + ".vcf.gz"
    run:
        op.relative_symlink(input.vcf1, output.vcf1),
        op.relative_symlink(input.vcf2, output.vcf2),
        op.relative_symlink(input.vcf1 + ".tbi" , output.vcf1 + ".tbi"),
        op.relative_symlink(input.vcf2 + ".tbi" , output.vcf2 + ".tbi")

rule _starfish_run:
    input:
        vcf1 = str(rules._starfish_input_vcf.output.vcf1),
        vcf2 = str(rules._starfish_input_vcf.output.vcf2),
        reference = CFG["inputs"]["reference"],
        starfish_script = CFG["inputs"]["starfish_script"]
    output:
        venn = run_starfish_base + "venn.pdf",
        tool1_only = run_starfish_base + "A.vcf.gz",
        tool1_only_tbi = run_starfish_base + "A.vcf.gz.tbi",
        tool2_only = run_starfish_base + "B.vcf.gz",
        tool2_only_tbi = run_starfish_base + "B.vcf.gz.tbi",
        intersect = run_starfish_base + "A_and_B.vcf.gz",
        intersect_tbi = run_starfish_base + "A_and_B.vcf.gz.tbi",
        completed = run_starfish_base + "starfish_run.complete"
    log:
        stdout = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
        stderr = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
    params:
        tool1 = tool1,
        tool2 = tool2,
        vcf_dir = run_starfish_base
    conda:
        CFG["conda_envs"]["starfish"]
    threads:
        CFG["threads"]["starfish"]
    resources:
        mem_mb = CFG["mem_mb"]["starfish"]
    shell:
       op.as_one_line("""
        {input.starfish_script} --sdf {input.reference} -O {params.vcf_dir}
        --names {params.tool1} {params.tool2}
        --sample ALT --squash-ploidy --vennout {output.venn}
        -V {input.vcf1} {input.vcf2} > {log.stdout} 2> {log.stderr} && touch {output.completed}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _starfish_output_vcf:
    input:
        tool1_only = str(rules._starfish_run.output.tool1_only),
        tool2_only = str(rules._starfish_run.output.tool2_only),
        intersect = str(rules._starfish_run.output.intersect),
        starfish_venn = str(rules._starfish_run.output.venn),
        completed = str(rules._starfish_run.output.completed),
    output:
        t1 = output_base_vcf + tool1 + "-unique.vcf.gz",
        t2 = output_base_vcf + tool2 + "-unique.vcf.gz",
        isec = output_base_vcf + tool1 + "-and-" + tool2 + ".vcf.gz",
        completed = output_base_vcf + "starfish_complete"
    run:
        op.relative_symlink(input.tool1_only, output.t1, in_module = True),
        op.relative_symlink(input.tool1_only + ".tbi", output.t1 + ".tbi", in_module = True),
        op.relative_symlink(input.tool2_only, output.t2, in_module = True),
        op.relative_symlink(input.tool2_only + ".tbi", output.t2 + ".tbi", in_module = True),
        op.relative_symlink(input.intersect, output.isec, in_module = True),
        op.relative_symlink(input.intersect + ".tbi", output.isec + ".tbi", in_module = True),
        op.relative_symlink(input.completed, output.completed, in_module = True),

#should generalize for all VCFs to avoid redundancy. Note the need for the Strelka indels so there are additional outputs here.
#This rule keeps all indels from both tools (i.e not just Strelka)
rule _starfish_vcf_to_bed:
    input:
        tool1_only = str(rules._starfish_run.output.tool1_only),
        tool2_only = str(rules._starfish_run.output.tool2_only),
        intersect = str(rules._starfish_run.output.intersect),
        completed = str(rules._starfish_run.output.completed),
        vcf1 = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}." + tool1 + ".vcf.gz",
        vcf2 = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}." + tool2 + ".vcf.gz"
    output:
        tool1_only = vcf_to_bed_base + "A.bed",
        tool2_only = vcf_to_bed_base + "B.bed",
        intersect = vcf_to_bed_base + "A_and_B.bed",
        tool1_only_indel_bed = vcf_to_bed_base + "A.indels.bed",
        tool2_only_indel_bed = vcf_to_bed_base + "B.indels.bed",
        intersect_plus_indels = vcf_to_bed_base + "A_and_B_plus_indels.bed",
        union = vcf_to_bed_base + "A_plus_B.bed"
    params: base_dir = vcf_to_bed_base
    log:
        stderr = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/vcf_to_bed.stderr.log"
    conda:
        CFG["conda_envs"]["bedops"]
    threads:
        CFG["threads"]["vcf_to_bed"]
    resources:
        mem_mb = CFG["mem_mb"]["vcf_to_bed"]
    shell:
        op.as_one_line("""
        sort -u -S {resources.mem_mb}M  <(zcat {input.vcf1} | vcf2bed | cut -f 1-3) <(zcat {input.vcf2} | vcf2bed | cut -f 1-3) | sort -S {resources.mem_mb}M -k1,1 -k2,2n > {output.union} ;
        zcat {input.tool1_only} | vcf2bed | cut -f 1-3 > {output.tool1_only} 2>> {log.stderr};
        zcat {input.tool2_only} | vcf2bed | cut -f 1-3 > {output.tool2_only} 2>> {log.stderr};
        zcat {input.intersect} | vcf2bed | cut -f 1-3 > {output.intersect} 2>> {log.stderr};
        zcat {input.tool1_only} | awk 'length($4)>1 || length($5)>1' | vcf2bed | cut -f 1-3 > {output.tool1_only_indel_bed} 2>> {log.stderr};
        zcat {input.tool2_only} | awk 'length($4)>1 || length($5)>1' | vcf2bed | cut -f 1-3 > {output.tool2_only_indel_bed} 2>> {log.stderr};
        cat {output.tool1_only_indel_bed} {output.tool2_only_indel_bed} {output.intersect} | sort -S {resources.mem_mb}M -k1,1 -k2,2n > {output.intersect_plus_indels} 2>> {log.stderr};
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _starfish_output_bed:
    input:
        tool1_only = str(rules._starfish_vcf_to_bed.output.tool1_only),
        tool2_only = str(rules._starfish_vcf_to_bed.output.tool2_only),
        intersect = str(rules._starfish_vcf_to_bed.output.intersect),
        tool1_only_indel_bed = str(rules._starfish_vcf_to_bed.output.tool1_only_indel_bed),
        tool2_only_indel_bed = str(rules._starfish_vcf_to_bed.output.tool2_only_indel_bed),
        intersect_plus_indels = str(rules._starfish_vcf_to_bed.output.intersect_plus_indels)
    output:
        tool1_only = output_base_bed + tool1 + "-unique.bed",
        tool2_only = output_base_bed + tool2 + "-unique.bed",
        intersect = output_base_bed + tool1 + "-and-" + tool2 + ".bed",
        tool1_only_indel_bed = output_base_bed + tool1 + "-unique.indels.bed",
        tool2_only_indel_bed = output_base_bed + tool2 + "-unique.indels.bed",
        intersect_plus_indels = output_base_bed + tool1 + "-and-" + tool2 + "-and-all-indels.bed"
    run:
        op.relative_symlink(input.tool1_only, output.tool1_only),
        op.relative_symlink(input.tool2_only, output.tool2_only),
        op.relative_symlink(input.intersect, output.intersect),
        op.relative_symlink(input.tool1_only_indel_bed, output.tool1_only_indel_bed),
        op.relative_symlink(input.tool2_only_indel_bed, output.tool2_only_indel_bed),
        op.relative_symlink(input.intersect_plus_indels, output.intersect_plus_indels)

# Generates the target sentinels for each run, which generate the symlinks
rule _starfish_all:
    input:
        expand(
            [
                str(rules._starfish_output_vcf.output.completed),
                str(rules._starfish_output_bed.output.intersect_plus_indels)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            normal_id = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"])

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
