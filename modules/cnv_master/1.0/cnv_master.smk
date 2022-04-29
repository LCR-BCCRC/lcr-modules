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
# `CFG` is a shortcut to `config["lcr-modules"]["cnv_master"]`
CFG = op.setup_module(
    name = "cnv_master",
    version = "1.0",
    subdirectories = ["inputs", "merges", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _cnv_master_input,
    _cnv_master_all

callers = list(CFG["inputs"].keys())
callers = [caller.lower() for caller in callers]


##### RULES #####

# this function will find the existing outputs for every sample in the sample table
# it will preferrentially use the output from the tool according to user-sepcified order
def _find_best_seg(wildcards):
    CFG = config["lcr-modules"]["cnv_master"]
    possible_outputs = []
    for caller in callers:
        possible_outputs = possible_outputs + (expand(CFG["inputs"][caller],
                              tumour_id=wildcards.tumour_id,
                              normal_id=wildcards.normal_id,
                              pair_status=wildcards.pair_status,
                              tool=CFG["inputs"].keys(),
                              projection=wildcards.projection,
                              seq_type=list(CFG["runs"]["tumour_seq_type"].unique())))
    possible_outputs = [path for path in possible_outputs if os.path.exists(path)]
    assert (len(possible_outputs) >= 1), (
        f"No ouput was found for the sample {wildcards.tumour_id} in projection {wildcards.projection}. "
        f"Please ensure it exists at one of the paths specified through config."
    )
    print(f"Will use file {possible_outputs[0]} for the sample {wildcards.tumour_id}.")
    return(possible_outputs[0])

# symlink the input files into 00-inputs
rule _cnv_master_input:
    input:
        seg = _find_best_seg
    output:
        seg = CFG["dirs"]["inputs"] + "seg/{tumour_id}--{normal_id}--{pair_status}--{seq_type}--{projection}.seg"
    run:
        op.absolute_symlink(input.seg, output.seg)


def _get_all_segs(this_seq_type):
    CFG = config["lcr-modules"]["cnv_master"]
    all_segs = op.filter_samples(CFG["runs"], tumour_seq_type = this_seq_type)
    this_set = expand(str(rules._cnv_master_input.output.seg),
                            zip,
                            allow_missing = True,
                            tumour_id=all_segs["tumour_sample_id"],
                            normal_id=all_segs["normal_sample_id"],
                            pair_status=all_segs["pair_status"],
                            seq_type=all_segs["tumour_seq_type"])
    return(this_set)


# Generates the merged files for the genomes if they are present in sample table
rule _cnv_master_merge_genome_projections:
    input:
        seg_file = _get_all_segs("genome")
    output:
        merge = CFG["dirs"]["merges"] + "{seq_type}/projection--{projection}.seg",
        contents = CFG["dirs"]["merges"] + "{seq_type}/projection--{projection}.contents"
    wildcard_constraints:
        seq_type="genome"
    conda:
        CFG["conda_envs"]["R"]
    threads:
        CFG["threads"]["cnv_master"]
    resources:
        **CFG["resources"]["cnv_master"]
    group: "cnv_master"
    script:
        "src/R/merge_segs.R"


# Generates the merged files for the captures if they are present in sample table
rule _cnv_master_merge_capture_projections:
    input:
        seg_file = _get_all_segs("capture")
    output:
        merge = CFG["dirs"]["merges"] + "{seq_type}/projection--{projection}.seg",
        contents = CFG["dirs"]["merges"] + "{seq_type}/projection--{projection}.contents"
    wildcard_constraints:
        seq_type="capture"
    conda:
        CFG["conda_envs"]["R"]
    threads:
        CFG["threads"]["cnv_master"]
    resources:
        **CFG["resources"]["cnv_master"]
    group: "cnv_master"
    script:
        "src/R/merge_segs.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _cnv_master_output_genome_merges:
    input:
        genome_merge = str(rules._cnv_master_merge_genome_projections.output.merge),
        genome_content = str(rules._cnv_master_merge_genome_projections.output.contents)
    output:
        genome_merge = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}.seg",
        genome_content = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}.contents"
    wildcard_constraints:
        seq_type="genome"
    group: "cnv_master"
    run:
        op.relative_symlink(input.genome_merge, output.genome_merge, in_module = True)
        op.relative_symlink(input.genome_content, output.genome_content, in_module = True)


rule _cnv_master_output_capture_merges:
    input:
        capture_merge = str(rules._cnv_master_merge_capture_projections.output.merge),
        capture_content = str(rules._cnv_master_merge_capture_projections.output.contents)
    output:
        capture_merge = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}.seg",
        capture_content = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}.contents"
    wildcard_constraints:
        seq_type="capture"
    group: "cnv_master"
    run:
        op.relative_symlink(input.capture_merge, output.capture_merge, in_module = True)
        op.relative_symlink(input.capture_content, output.capture_content, in_module = True)

rule _cnv_master_all:
    input:
        expand(
            [
                str(rules._cnv_master_output_genome_merges.output.genome_merge),
                str(rules._cnv_master_output_genome_merges.output.genome_content),
                str(rules._cnv_master_output_capture_merges.output.capture_merge),
                str(rules._cnv_master_output_capture_merges.output.capture_content)
            ],
            projection=CFG["projections"],
            seq_type=list(CFG["runs"]["tumour_seq_type"].unique())
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
