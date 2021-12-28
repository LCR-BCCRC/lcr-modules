#!/usr/bin/env snakemake

##### ATTRIBUTION #####


# Original Author:  Bruno Grande
# Module Author:    Bruno Grande
# Contributors:     Kostia Dreval, Ryan Morin


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
# `CFG` is a shortcut to `config["lcr-modules"]["lofreq"]`
CFG = op.setup_module(
    name = "lofreq",
    version = "1.0",
    subdirectories = ["inputs", "lofreq", "combined", "filtered", "outputs"],
)
#set variable for prepending to PATH based on config
SCRIPT_PATH = CFG['inputs']['src_dir']
#this is used in place of the shell.prefix() because that was not working consistently. This is not ideal. 

# Define rules to be run locally when using a compute cluster
localrules:
    _lofreq_input_bam,
    _lofreq_output_vcf,
    _lofreq_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _lofreq_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Run LoFreq in somatic variant calling mode
rule _lofreq_run:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dbsnp = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz"), #in our experience, this filter doesn't remove as many SNPs as one would expect
        bed = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.bed")
    output:
        vcf_snvs_filtered = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.snvs.vcf.gz",
        vcf_indels_filtered = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.indels.vcf.gz",
        vcf_snvs_all = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.snvs.vcf.gz",
        vcf_indels_all = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.indels.vcf.gz"
    log:
        stdout = CFG["logs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.stdout.log",
        stderr = CFG["logs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.stderr.log"
    params:
        opts = CFG["options"]["lofreq"],
        regions = op.switch_on_wildcard("seq_type", CFG["switches"]["regions_bed"]),
        rm_files = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/normal_relaxed.vcf.gz"
    conda:
        CFG["conda_envs"]["lofreq"]
    threads:
        CFG["threads"]["lofreq"]
    resources:
        **CFG["resources"]["lofreq"]
    shell:
        op.as_one_line("""
        SCRIPT_PATH={SCRIPT_PATH};
        PATH=$SCRIPT_PATH:$PATH;
        SCRIPT="$SCRIPT_PATH/lofreq2_call_pparallel.py";
        if [[ $(which lofreq2_call_pparallel.py) =~ $SCRIPT ]]; then 
            echo "using bundled patched script $SCRIPT";
            if [ ! -e {output.vcf_snvs_filtered} ] && [ -e {params.rm_files} ]; then rm $(dirname {output.vcf_snvs_all})/*; fi
            &&
            lofreq somatic {params.opts} --threads {threads} -t {input.tumour_bam} -n {input.normal_bam}
            -f {input.fasta} -o $(dirname {output.vcf_snvs_filtered})/ -d {input.dbsnp} {params.regions} --bed {input.bed}
            > {log.stdout} 2> {log.stderr}
            &&
            rm {params.rm_files};
        else echo "WARNING: PATH is not set properly, using $(which lofreq2_call_pparallel.py)"; fi
        """)

# indels are not yet called but this rule merges the empty indels file with the snvs file to produce the consistently named "combined" vcf. 
rule _lofreq_combine_vcf:
    input:
        vcf_all = expand(CFG["dirs"]["lofreq"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/somatic_final.{var_type}.vcf.gz",
                    var_type = ["indels", "snvs"]),
        vcf_all_filtered = expand(CFG["dirs"]["lofreq"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/somatic_final_minus-dbsnp.{var_type}.vcf.gz",
                    var_type = ["indels", "snvs"])
    output:
        vcf_all = temp(CFG["dirs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.combined.vcf.gz"),
        vcf_all_filtered = temp(CFG["dirs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.combined.vcf.gz"),
    resources:
        **CFG["resources"]["bcftools_sort"]
    conda:
        CFG["conda_envs"]["lofreq"]
    log:
        stdout_all = CFG["logs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_final.combined.stdout.log",
        stderr_all = CFG["logs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_final.combined.stderr.log",
        stdout_all_filtered = CFG["logs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_final_minus-dbsnp.combined.stdout.log",
        stderr_all_filtered = CFG["logs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_final_minus-dbsnp.combined.stderr.log"
    shell:
                op.as_one_line("""
        bcftools concat -a {input.vcf_all} |
        bcftools sort --max-mem {resources.mem_mb}M -Oz -o {output.vcf_all}
        > {log.stdout_all} 2> {log.stderr_all}
            &&
        bcftools concat -a {input.vcf_all_filtered} |
        bcftools sort --max-mem {resources.mem_mb}M -Oz -o {output.vcf_all_filtered}
        > {log.stdout_all_filtered} 2> {log.stderr_all_filtered}
        """)


rule _lofreq_filter_vcf:
    input:
        vcf_all = rules._lofreq_combine_vcf.output.vcf_all,
        vcf_all_filtered = rules._lofreq_combine_vcf.output.vcf_all_filtered
    output:
        vcf_all_clean = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.combined.filtered.vcf.gz",
        vcf_all_filtered_clean = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.combined.filtered.vcf.gz",
        vcf_all_clean_tbi = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.combined.filtered.vcf.gz.tbi",
        vcf_all_filtered_clean_tbi = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.combined.filtered.vcf.gz.tbi"
    resources:
        **CFG["resources"]["bcftools_sort"]
    conda:
        CFG["conda_envs"]["lofreq"]
    shell:
        op.as_one_line("""
        PATH={SCRIPT_PATH}:$PATH;
        SCRIPT=$(which lofreq_filter.sh); 
        echo "using bundled custom filtering script $SCRIPT";
        lofreq_filter.sh {input.vcf_all} | bgzip > {output.vcf_all_clean}
          && tabix -p vcf {output.vcf_all_clean}
              &&
        lofreq_filter.sh {input.vcf_all_filtered} | bgzip > {output.vcf_all_filtered_clean}
          && tabix -p vcf {output.vcf_all_filtered_clean}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _lofreq_output_vcf:
    input:
        vcf_all = rules._lofreq_filter_vcf.output.vcf_all_clean,
        vcf_all_filtered = rules._lofreq_filter_vcf.output.vcf_all_filtered_clean
    output:
        vcf_all = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lofreq.snvs.vcf.gz",
        vcf_all_filtered = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_minus-dbsnp.lofreq.snvs.vcf.gz"
    run:
        op.relative_symlink(input.vcf_all, output.vcf_all, in_module=True)
        op.relative_symlink(input.vcf_all + ".tbi", output.vcf_all + ".tbi", in_module=True)
        op.relative_symlink(input.vcf_all_filtered, output.vcf_all_filtered, in_module=True)
        op.relative_symlink(input.vcf_all_filtered + ".tbi", output.vcf_all_filtered + ".tbi", in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _lofreq_all:
    input:
        expand(
            [
                rules._lofreq_output_vcf.output.vcf_all,
                rules._lofreq_output_vcf.output.vcf_all_filtered,
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
