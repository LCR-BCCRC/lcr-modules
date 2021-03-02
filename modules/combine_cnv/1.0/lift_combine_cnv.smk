#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
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
# `CFG` is a shortcut to `config["lcr-modules"]["combine_cnv"]`
CFG = op.setup_module(
    name = "combine_cnv",
    version = "1.0",
    subdirectories = ["inputs", "filter_cnv", "combine_cnv", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _combine_cnv_input,
    _combine_cnv_filteR,
    _combine_cnv_seg2bed,
    _combine_cnv_fill_segments,
    _combine_cnv_output_seg,
    _combine_cnv_all,


CFG["inputs"]["seg"] = dict(zip(CFG["inputs"]["names"], CFG["inputs"]["sample_seg"]))
callers = list(CFG["inputs"]["seg"].keys())
callers = [caller.lower() for caller in callers]


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _combine_cnv_input:
    input:
        seg = lambda w: config["lcr-modules"]["combine_cnv"]["inputs"]["seg"][w.caller],
        vcf = CFG["inputs"]["sample_vcf"]
    output:
        seg = CFG["dirs"]["inputs"] + "seg/{caller}/{seq_type}/{tumour_id}--{normal_id}--{pair_status}.seg",
        vcf = CFG["dirs"]["inputs"] + "vcf/for_{caller}/{seq_type}/{tumour_id}--{normal_id}--{pair_status}.vcf.gz"
    run:
        op.absolute_symlink(input.seg, output.seg)
        op.absolute_symlink(input.vcf, output.vcf)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _combine_cnv_filteR:
    input:
        seg_file = str(rules._combine_cnv_input.output.seg),
        vcf_file = str(rules._combine_cnv_input.output.vcf),
        blacklist = reference_files("genomes/grch37/encode/encode-blacklist.grch37.bed")
    output:
        cnv = temp(CFG["dirs"]["filter_cnv"] + "{caller}/{seq_type}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.CNV.tsv"),
        seg = temp(CFG["dirs"]["filter_cnv"] + "{caller}/{seq_type}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.filtered.seg")
    params:
        genome_build="hg19"
    conda:
        CFG["conda_envs"]["CNVfilteR"]
    threads:
        CFG["threads"]["combine_cnv"]
    resources:
        **CFG["resources"]["combine_cnv"]
    script:
        "src/R/CNVfilteR.R"


rule _combine_cnv_seg2bed:
    input:
        seg = str(rules._combine_cnv_filteR.output.seg),
        cnv = str(rules._combine_cnv_filteR.output.cnv)
    output:
        bed = temp(CFG["dirs"]["filter_cnv"] + "from--grch37/{seq_type}/{tumour_id}--{normal_id}--{pair_status}.{caller}.bed"),
        header = temp(CFG["dirs"]["filter_cnv"] + "from--grch37/{seq_type}/{tumour_id}--{normal_id}--{pair_status}.{caller}.bed.header")
    params:
        opts = CFG["options"]["seg2bed2seg"],
        chr_colNum = 2,
        start_colNum = 3,
        end_colNum = 4,
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.seg} 
        --output {output.bed} 
        --chromColnum {params.chr_colNum} 
        --startColnum {params.start_colNum} 
        --endColnum {params.end_colNum}
        """)


# Convert the bed file in hg19 coordinates into hg38 coordinates
rule _run_liftover:
    input:
        native = rules._combine_cnv_seg2bed.output.bed,
        chains = reference_files("genomes/grch37/chains/grch37/hg19ToHg38.over.chain")
    output:
        lifted = temp(CFG["dirs"]["filter_cnv"] + "from--grch37/{seq_type}/{tumour_id}--{normal_id}--{pair_status}.{caller}.lifted_hg19ToHg38.bed"),
        unmapped = temp(CFG["dirs"]["filter_cnv"] + "from--grch37/{seq_type}/{tumour_id}--{normal_id}--{pair_status}.{caller}.lifted_hg19ToHg38.unmapped.bed")
    params:
        mismatch = CFG["options"]["min_mismatch"]
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        liftOver -minMatch={params.mismatch}
        {input.native} {input.chains} 
        {output.lifted} {output.unmapped}
        """)

# Sort liftover output
# Here, the perl line will filter out non-standard chromosomes from the output
rule _liftover_sort:
    input:
        lifted = rules._run_liftover.output.lifted
    output:
        lifted_sorted = temp(CFG["dirs"]["filter_cnv"] + "from--grch37/{seq_type}/{tumour_id}--{normal_id}--{pair_status}.{caller}.lifted_hg19ToHg38.sorted.bed")
    shell:
        op.as_one_line("""
        sort -k1,1 -k2,2n -V {input.lifted} |
        perl -ne 'print if /^(chr)*[\dX]+\s.+/'
        > {output.lifted_sorted}
        """)


# Convert the bed file in hg19 coordinates into seg format
rule _liftover_bed2seg:
    input:
        lifted_sorted = str(rules._liftover_sort.output.lifted_sorted),
        headers = str(rules._combine_cnv_seg2bed.output.header)
    output:
        seg_lifted = temp(CFG["dirs"]["filter_cnv"] + "from--grch37/{seq_type}/{tumour_id}--{normal_id}--{pair_status}.{caller}.lifted_hg19ToHg38.seg")
    params:
        opts = CFG["options"]["seg2bed2seg"]
    conda:
        CFG["conda_envs"]["liftover-366"]
    shell:
        op.as_one_line("""
        python {params.opts} 
        --input {input.lifted_sorted}
        --column-header {input.headers}
        --output {output.seg_lifted} 
        """)

rule _combine_cnv_cleanup:
    input:
        seg = str(rules._liftover_bed2seg.output.seg_lifted)
    output:
        seg_complete = CFG["dirs"]["filter_cnv"] + "from--grch37/{seq_type}/{tumour_id}--{normal_id}--{pair_status}.{caller}.complete.lifted_hg19ToHg38.seg"
    conda:
        CFG["conda_envs"]["CNVfilteR"]
    threads:
        CFG["threads"]["combine_cnv"]
    resources:
        **CFG["resources"]["combine_cnv"]
    script:
        "src/R/cleanup.R"



# Example variant filtering rule (single-threaded; can be run on cluster head node)
rule _combine_cnv_fill_segments:
    input:
        seg = str(rules._combine_cnv_cleanup.output.seg_complete)
    output:
        seg = CFG["dirs"]["combine_cnv"] + "{caller}/{seq_type}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.filtered.filled.{caller}.seg"
    params:
        chromArm = CFG["options"]["chromArm"]
    conda:
        CFG["conda_envs"]["fill_segments"]
    threads:
        CFG["threads"]["combine_cnv"]
    resources:
        **CFG["resources"]["combine_cnv"]
    script:
        "src/python/fill.py"


# Example variant filtering rule (single-threaded; can be run on cluster head node)
rule _combine_cnv_merge_segs:
    input:
        seg_file = expand(CFG["dirs"]["combine_cnv"] + "{{caller}}/{{seq_type}}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.filtered.filled.{{caller}}.seg", zip,
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                pair_status=CFG["runs"]["pair_status"])

    output:
        seg = CFG["dirs"]["combine_cnv"] + "{caller}/{seq_type}/{caller}--merged.filtered.filled.seg"
    params:
        concatenate = CFG["options"]["concatenate_segs"],
        input_dir = CFG["dirs"]["combine_cnv"] + "{caller}/{seq_type}"
    threads:
        CFG["threads"]["combine_cnv"]
    resources:
        **CFG["resources"]["combine_cnv"]
    shell:
        op.as_one_line("""
        bash {params.concatenate} {params.input_dir} {output.seg} 
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _combine_cnv_output_seg:
    input:
        merged_seg = str(rules._combine_cnv_merge_segs.output.seg),
    output:
        merged_seg = CFG["dirs"]["outputs"] + "{caller}/{seq_type}/{caller}--merged.filtered.filled.seg"
    run:
        op.relative_symlink(input.merged_seg, output.merged_seg, in_module=True)



# Generates the target sentinels for each run, which generate the symlinks
rule _combine_cnv_all:
    input:
        expand(
            [
                str(rules._combine_cnv_output_seg.output.merged_seg)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            caller = callers)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
