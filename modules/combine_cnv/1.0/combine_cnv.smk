#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
#import glob

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
    _combine_cnv_input_matched,
    _combine_cnv_input_unmatched
    _combine_cnv_seg2bed,
    _combine_cnv_fill_segments,
    _combine_cnv_output_seg,
    _combine_cnv_all,


CFG["inputs"]["seg"] = dict(zip(CFG["inputs"]["names"], CFG["inputs"]["sample_seg"]))
callers = list(CFG["inputs"]["seg"].keys())
callers = [caller.lower() for caller in callers]

#for caller in callers:
#    CFG = config["lcr-modules"]["combine_cnv"]
#    path = CFG["inputs"]["test_key"] + caller + "-" + CFG["module_versions"][caller] + "/99-outputs/"
#    print(path)
#    existing_outputs = glob.glob(path + r'/**/*.seg', recursive=True)
#    print(existing_outputs)

#quit()

##### RULES #####

def _combine_cnv_get_bams(wildcards):
    CFG = config["lcr-modules"]["combine_cnv"]
    tbl = CFG["runs"]
    this_pair_status = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["pair_status"]
    if str(this_pair_status[0]) == "matched":
        supported_callers = ["battenberg", "sequenza"]
        this_caller = [value for value in callers if value in supported_callers]
        print(this_caller)
        print(type(this_caller))

    print(this_path)
    for i in this_caller:
        print("new: " + i)
        return(str(i))
    


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _combine_cnv_input_matched:
    input:
        seg = lambda w: config["lcr-modules"]["combine_cnv"]["inputs"]["seg"][w.caller],
        #seg = _combine_cnv_get_bams,
    output:
        seg = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{caller}/{tumour_id}--{normal_id}--{pair_status}.seg"
    wildcard_constraints:
        pair_status = "matched",
        caller = "battenberg|sequenza"
    run:
        op.absolute_symlink(input.seg, output.seg)
        op.absolute_symlink(input.vcf, output.vcf)

rule _combine_cnv_input_unmatched:
    input:
        seg = lambda w: config["lcr-modules"]["combine_cnv"]["inputs"]["seg"][w.caller]
    output:
        seg = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{caller}/{tumour_id}--{normal_id}--{pair_status}.seg"
    wildcard_constraints:
        pair_status = "unmatched",
        caller = "controlfreec|cnvkit"
    run:
        op.absolute_symlink(input.seg, output.seg)
        op.absolute_symlink(input.vcf, output.vcf)


rule _combine_cnv_seg2bed:
    input:
        seg = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{caller}/{tumour_id}--{normal_id}--{pair_status}.seg"
    output:
        bed = temp(CFG["dirs"]["filter_cnv"] + "{seq_type}--{genome_build}/{caller}/{tumour_id}--{normal_id}--{pair_status}.{caller}.bed"),
        header = temp(CFG["dirs"]["filter_cnv"] + "{seq_type}--{genome_build}/{caller}/{tumour_id}--{normal_id}--{pair_status}.{caller}.bed.header")
    params:
        opts = config["lcr-modules"]["_shared"]["lcr-modules"] + "/modules/liftover/2.0/src/convert_for_liftover.py",
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

def get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")

# Convert the bed file in hg19 coordinates into hg38 coordinates
rule _combine_cnv_run_liftover:
    input:
        native = rules._combine_cnv_seg2bed.output.bed,
        chains = get_chain
    output:
        lifted = temp(CFG["dirs"]["filter_cnv"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{caller}.{genome_build}.bed"),
        unmapped = temp(CFG["dirs"]["filter_cnv"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{caller}.{genome_build}.unmapped.bed")
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
rule _combine_cnv_liftover_sort:
    input:
        lifted = rules._combine_cnv_run_liftover.output.lifted
    output:
        lifted_sorted = temp(CFG["dirs"]["filter_cnv"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{caller}.{genome_build}.sorted.bed")
    shell:
        op.as_one_line("""
        sort -k1,1 -k2,2n -V {input.lifted} |
        perl -ne 'print if /^(chr)*[\dX]+\s.+/'
        > {output.lifted_sorted}
        """)


# Convert the bed file in hg19 coordinates into seg format
rule _combine_cnv_liftover_bed2seg:
    input:
        lifted_sorted = str(rules._combine_cnv_liftover_sort.output.lifted_sorted),
        headers = str(rules._combine_cnv_seg2bed.output.header)
    output:
        seg_lifted = temp(CFG["dirs"]["filter_cnv"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{caller}.{genome_build}.seg")
    params:
        opts = config["lcr-modules"]["_shared"]["lcr-modules"] + "/modules/liftover/2.0/src/convert_for_liftover.py"
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
        seg = str(rules._combine_cnv_liftover_bed2seg.output.seg_lifted)
    output:
        seg_complete = CFG["dirs"]["filter_cnv"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{caller}.complete.{genome_build}.seg"
    conda:
        CFG["conda_envs"]["CNVfilteR"]
    threads:
        CFG["threads"]["combine_cnv"]
    resources:
        **CFG["resources"]["combine_cnv"]
    script:
        "src/R/cleanup.R"


def get_chrArm(wildcards):
    if "38" in str({wildcards.genome_build}):
        return config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/1.0/src/chromArm.hg38.tsv"
    else:
        return config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/1.0/src/chromArm.hg19.tsv"

# Example variant filtering rule (single-threaded; can be run on cluster head node)
rule _combine_cnv_fill_segments:
    input:
        seg = str(rules._combine_cnv_cleanup.output.seg_complete),
        chrArm = get_chrArm
    output:
        seg = CFG["dirs"]["combine_cnv"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.filtered.filled.{caller}.{genome_build}.seg"
    params:
        chrArm = get_chrArm
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
        seg_file = expand(CFG["dirs"]["combine_cnv"] + "{{seq_type}}--projection/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.filtered.filled.{{caller}}.{{genome_build}}.seg", zip,
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                pair_status=CFG["runs"]["pair_status"])

    output:
        seg = CFG["dirs"]["combine_cnv"] + "{seq_type}--projection/{caller}--merged.filtered.filled.{genome_build}.seg"
    params:
        concatenate = CFG["options"]["concatenate_segs"],
        input_dir = CFG["dirs"]["combine_cnv"] + "{seq_type}--projection"
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
        merged_seg = CFG["dirs"]["outputs"] + "{seq_type}--projection/{caller}--merged.filtered.filled.{genome_build}.seg"
    run:
        op.relative_symlink(input.merged_seg, output.merged_seg, in_module=True)



# Generates the target sentinels for each run, which generate the symlinks
rule _combine_cnv_merge:
    input:
        expand(
            [
                CFG["dirs"]["outputs"] + "{seq_type}--projection/{caller}--merged.filtered.filled.{{genome_build}}.seg"
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"].drop_duplicates().tolist()*len(CFG["projections"]),
            caller = callers*len(CFG["projections"]))
    output:
        CFG["dirs"]["outputs"] + "master_merge/merged.filtered.filled.{genome_build}.seg"
    shell:
        "cat input > output"

rule _combine_cnv_all:
    input:
        expand(
            [
                str(rules._combine_cnv_merge.output)
            ],
            zip,
            genome_build=CFG["projections"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
