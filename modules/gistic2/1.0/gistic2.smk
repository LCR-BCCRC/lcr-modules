#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Sierra Gillis
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
# `CFG` is a shortcut to `config["lcr-modules"]["gistic2"]`
CFG = op.setup_module(
    name = "gistic2",
    version = "1.0",
    subdirectories = ["inputs", "prepare_seg", "standard_chr", "markers", "gistic2", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _gistic2_input_seg,
    _gistic2_input_sample_sets,
    _gistic2_prepare_seg,
    _gistic2_standard_chr,
    _gistic2_make_markers,
    _gistic2_download_ref,
    _gistic2_output,
    _gistic2_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _gistic2_input_seg:
    input:
        seg = CFG["inputs"]["seg"]
    output:
        seg = CFG["dirs"]["inputs"] + "{seq_type}--projection/all--{projection}.seg"
    group: 
        "input_and_format"
    run:
        op.absolute_symlink(input.seg, output.seg)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _gistic2_input_sample_sets:
    input:
        all_sample_sets = CFG["inputs"]["all_sample_sets"]
    output:
        all_sample_sets = CFG["dirs"]["inputs"] + "sample_sets/all_sample_sets.tsv"
    group: 
        "input_and_format"
    run:
        op.absolute_symlink(input.all_sample_sets, output.all_sample_sets)

# Download refgene reference MAT file
rule _gistic2_download_ref:
    output:
        refgene_mat = CFG["dirs"]["inputs"] + "references/{projection}.refgene.mat"
    params:
        url = "http://bcgsc.ca/downloads/morinlab/gistic2_references/{projection}.refgene.mat",
        folder = CFG["dirs"]["inputs"] + "references"
    shell:
        op.as_one_line("""
        wget -P {params.folder} {params.url} 
        """)

def _get_seg_input_params(input_files = [str(rules._gistic2_input_seg.output.seg)]):
    # match "genome" and/or "capture" in the input file names, add flag to string
    genome_param = ["--genome " + v for v in input_files if "genome" in v]
    capture_param = ["--capture " + v for v in input_files if "capture" in v]

    # combine into one string
    return(" ".join(genome_param + capture_param))

# Merges capture and genome seg files if available, and subset to the case_set provided
rule _gistic2_prepare_seg:
    input:
        seg = expand(str(rules._gistic2_input_seg.output.seg),
                    allow_missing=True,
                    seq_type=CFG["samples"]["seq_type"].unique()
                    ),
        all_sample_sets = str(rules._gistic2_input_sample_sets.output.all_sample_sets),
        seg_dir = CFG["dirs"]["inputs"]
    output:
        seg = CFG["dirs"]["prepare_seg"] + "{case_set}--{projection}.seg"
    log:
        stdout = CFG["logs"]["prepare_seg"] + "{case_set}--{projection}.stdout.log",
        stderr = CFG["logs"]["prepare_seg"] + "{case_set}--{projection}.stderr.log"
    params:
        script = CFG["prepare_seg"],
        case_set = CFG["case_set"],
        seg_input_params = _get_seg_input_params
    group: 
        "input_and_format"
    shell:
        op.as_one_line("""
        Rscript {params.script} 
        {params.seg_input_params} 
        --output_dir $(dirname {output.seg})/ 
        --all_sample_sets {input.all_sample_sets} 
        --case_set {params.case_set} 
        > {log.stdout} 2> {log.stderr}
        """)

# Removes entries for non-standard chromosomes from the seg file
rule _gistic2_standard_chr:
    input:
        seg = str(rules._gistic2_prepare_seg.output.seg)
    output:
        seg = CFG["dirs"]["standard_chr"] + "{case_set}--{projection}--standard_Chr.seg"
    log:
        stderr = CFG["logs"]["standard_chr"] + "{case_set}--{projection}--standard_Chr.stderr.log"
    shell:
        op.as_one_line("""
        awk 'BEGIN{{IGNORECASE = 1}} {{FS="\t"}} $2 !~ /Un|random|alt/ {{print}}' {input.seg} > {output.seg}
        2> {log.stderr}
        """)

# Create a markers file that has every segment start and end that appears in the seg file
rule _gistic2_make_markers:
    input:
        seg = str(rules._gistic2_standard_chr.output.seg)
    output:
        temp_markers = temp(CFG["dirs"]["markers"] + "temp_markers--{case_set}--{projection}.txt"),
        markers = CFG["dirs"]["markers"] + "markers--{case_set}--{projection}.txt"
    log:
        stderr = CFG["logs"]["markers"] + "{case_set}--{projection}--markers.stderr.log"
    shell: 
        op.as_one_line("""
        sed '1d' {input.seg} | cut -f2,3 > {output.temp_markers}
            &&
        sed '1d' {input.seg} | cut -f2,4 >> {output.temp_markers}
            &&
        sort -V -k1,1 -k2,2nr {output.temp_markers} | uniq | nl > {output.markers} 
        2> {log.stderr}
        """)

# Run gistic2 for a single seq_type (capture, genome) for the confidence thresholds listed in the config
rule _gistic2_run:
    input:
        seg = str(rules._gistic2_standard_chr.output.seg),
        refgene_mat = str(rules._gistic2_download_ref.output.refgene_mat),
        markers = str(rules._gistic2_make_markers.output.markers)
    output:
        all_data_by_genes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/all_data_by_genes.txt",
        all_lesions = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/all_lesions.conf_{conf}.txt",
        all_thresholded_by_genes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/all_thresholded.by_genes.txt",
        amp_genes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/amp_genes.conf_{conf}.txt",
        amp_qplot_pdf = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/amp_qplot.pdf",
        amp_qplot_png = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/amp_qplot.png",
        broad_data_by_genes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/broad_data_by_genes.txt",
        broad_significance_results = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/broad_significance_results.txt",
        broad_values_by_arm = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/broad_values_by_arm.txt",
        del_genes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/del_genes.conf_{conf}.txt",
        del_qplot_pdf = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/del_qplot.pdf",
        del_qplot_png = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/del_qplot.png",
        focal_data_by_genes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/focal_data_by_genes.txt",
        freqarms_vs_ngenes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/freqarms_vs_ngenes.pdf",
        raw_copy_number_pdf = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/raw_copy_number.pdf",
        raw_copy_number_png = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/raw_copy_number.png",
        regions_track = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/regions_track.conf_{conf}.bed",
        sample_cutoffs = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/sample_cutoffs.txt",
        sample_seg_counts = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/sample_seg_counts.txt",
        scores = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/scores.gistic"
    log:
        stdout = CFG["logs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/gistic2.stdout.log",
        stderr = CFG["logs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/gistic2.stderr.log"
    params:
        base_dir = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/conf_{conf}/",
        opts = CFG["options"]["gistic2_run"]
    conda:
        CFG["conda_envs"]["gistic2"]
    threads:
        CFG["threads"]["gistic2_run"]
    resources:
        **CFG["resources"]["gistic2_run"]
    shell:
        op.as_one_line("""
        gistic2 -b {params.base_dir} -seg {input.seg} -mk {input.markers} -refgene {input.refgene_mat} 
        -genegistic 1 -broad 1 -savegene 1 -conf 0.{wildcards.conf} -v 30 -saveseg 0 -savedata 0
        {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _gistic2_output:
    input:
        all_data_by_genes = str(rules._gistic2_run.output.all_data_by_genes),
        all_lesions = str(rules._gistic2_run.output.all_lesions),
        all_thresholded_by_genes = str(rules._gistic2_run.output.all_thresholded_by_genes),
        amp_genes = str(rules._gistic2_run.output.amp_genes),
        amp_qplot_pdf = str(rules._gistic2_run.output.amp_qplot_pdf),
        amp_qplot_png = str(rules._gistic2_run.output.amp_qplot_png),
        broad_data_by_genes = str(rules._gistic2_run.output.broad_data_by_genes),
        broad_significance_results = str(rules._gistic2_run.output.broad_significance_results),
        broad_values_by_arm = str(rules._gistic2_run.output.broad_values_by_arm),
        del_genes = str(rules._gistic2_run.output.del_genes),
        del_qplot_pdf = str(rules._gistic2_run.output.del_qplot_pdf),
        del_qplot_png = str(rules._gistic2_run.output.del_qplot_png),
        focal_data_by_genes = str(rules._gistic2_run.output.focal_data_by_genes),
        freqarms_vs_ngenes = str(rules._gistic2_run.output.freqarms_vs_ngenes),
        raw_copy_number_pdf = str(rules._gistic2_run.output.raw_copy_number_pdf),
        raw_copy_number_png = str(rules._gistic2_run.output.raw_copy_number_png),
        regions_track = str(rules._gistic2_run.output.regions_track),
        sample_cutoffs = str(rules._gistic2_run.output.sample_cutoffs),
        sample_seg_counts = str(rules._gistic2_run.output.sample_seg_counts),
        scores = str(rules._gistic2_run.output.scores)
    output:
        all_data_by_genes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/all_data_by_genes.txt",
        all_lesions = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/all_lesions.conf_{conf}.txt",
        all_thresholded_by_genes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/all_thresholded.by_genes.txt",
        amp_genes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/amp_genes.conf_{conf}.txt",
        amp_qplot_pdf = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/amp_qplot.pdf",
        amp_qplot_png = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/amp_qplot.png",
        broad_data_by_genes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/broad_data_by_genes.txt",
        broad_significance_results = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/broad_significance_results.txt",
        broad_values_by_arm = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/broad_values_by_arm.txt",
        del_genes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/del_genes.conf_{conf}.txt",
        del_qplot_pdf = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/del_qplot.pdf",
        del_qplot_png = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/del_qplot.png",
        focal_data_by_genes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/focal_data_by_genes.txt",
        freqarms_vs_ngenes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/freqarms_vs_ngenes.pdf",
        raw_copy_number_pdf = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/raw_copy_number.pdf",
        raw_copy_number_png = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/raw_copy_number.png",
        regions_track = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/regions_track.conf_{conf}.bed",
        sample_cutoffs = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/sample_cutoffs.txt",
        sample_seg_counts = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/sample_seg_counts.txt",
        scores = CFG["dirs"]["outputs"] + "{case_set}--{projection}/conf_{conf}/scores.gistic"

    run:
        op.relative_symlink(input.all_data_by_genes, output.all_data_by_genes, in_module= True),
        op.relative_symlink(input.all_lesions, output.all_lesionss, in_module= True),
        op.relative_symlink(input.all_thresholded_by_genes, output.all_thresholded_by_genes, in_module= True),
        op.relative_symlink(input.amp_genes, output.amp_genes, in_module= True),
        op.relative_symlink(input.amp_qplot_pdf, output.amp_qplot_pdf, in_module= True),
        op.relative_symlink(input.amp_qplot_png, output.amp_qplot_png, in_module= True),
        op.relative_symlink(input.broad_data_by_genes, output.broad_data_by_genes, in_module= True),
        op.relative_symlink(input.broad_significance_results, output.broad_significance_results, in_module= True),
        op.relative_symlink(input.broad_values_by_arm, output.broad_values_by_arm, in_module= True),
        op.relative_symlink(input.del_genes, output.del_genes, in_module= True),
        op.relative_symlink(input.del_qplot_pdf, output.del_qplot_pdf, in_module= True),
        op.relative_symlink(input.del_qplot_png, output.del_qplot_png, in_module= True),
        op.relative_symlink(input.focal_data_by_genes, output.focal_data_by_genes, in_module= True),
        op.relative_symlink(input.freqarms_vs_ngenes, output.freqarms_vs_ngenes, in_module= True),
        op.relative_symlink(input.raw_copy_number_pdf, output.raw_copy_number_pdf, in_module= True),
        op.relative_symlink(input.raw_copy_number_png, output.raw_copy_number_png, in_module= True),
        op.relative_symlink(input.regions_track, output.regions_track, in_module= True),
        op.relative_symlink(input.sample_cutoffs, output.sample_cutoffs, in_module= True),
        op.relative_symlink(input.sample_seg_counts, output.sample_seg_counts, in_module= True),
        op.relative_symlink(input.scores, output.scores, in_module= True)


rule _gistic2_all:
    input:
        expand(
            [
                str(rules._gistic2_output.output.all_data_by_genes),
                str(rules._gistic2_output.output.all_lesions),
                str(rules._gistic2_output.output.all_thresholded_by_genes),
                str(rules._gistic2_output.output.amp_genes),
                str(rules._gistic2_output.output.amp_qplot_pdf),
                str(rules._gistic2_output.output.amp_qplot_png),
                str(rules._gistic2_output.output.broad_data_by_genes),
                str(rules._gistic2_output.output.broad_significance_results),
                str(rules._gistic2_output.output.broad_values_by_arm),
                str(rules._gistic2_output.output.del_genes),
                str(rules._gistic2_output.output.del_qplot_pdf),
                str(rules._gistic2_output.output.del_qplot_png),
                str(rules._gistic2_output.output.focal_data_by_genes),
                str(rules._gistic2_output.output.freqarms_vs_ngenes),
                str(rules._gistic2_output.output.raw_copy_number_pdf),
                str(rules._gistic2_output.output.raw_copy_number_png),
                str(rules._gistic2_output.output.regions_track),
                str(rules._gistic2_output.output.sample_cutoffs),
                str(rules._gistic2_output.output.sample_seg_counts),
                str(rules._gistic2_output.output.scores)
            ],              
            conf = CFG["options"]["conf_level"],
            projection = CFG["projections"],
            case_set = CFG["case_set"]
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
