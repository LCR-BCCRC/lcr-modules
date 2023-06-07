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
    subdirectories = ["inputs", "markers", "gistic2", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _gistic2_input_seg,
    _gistic2_make_markers
    _gistic2_output,
    _gistic2_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _gistic2_input_seg:
    input:
        seg = CFG["inputs"]["seg"]
    output:
        seg = CFG["dirs"]["inputs"] + "{seq_type}--projection/all--{projection}.seg"
    run:
        op.absolute_symlink(input.seg, output.seg)

# Create a markers file that has every segment start and end that appears in the seg file
rule _gistic2_make_markers:
    input:
        seg = str(rules._gistic2_input_seg.output.seg)
    output:
        temp_markers = temp(CFG["dirs"]["markers"] + "{seq_type}--projection/temp_markers--{projection}.txt"),
        markers = CFG["dirs"]["markers"] + "{seq_type}--projection/markers--{projection}.txt"
    log:
        stderr = CFG["logs"]["markers"] + "{seq_type}--projection/all--{projection}/markers.stderr.log"
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
        seg = str(rules._gistic2_input_seg.output.seg),
        refgene_mat = "/home/sgillis/cancer_docker_singularity/gistic2/reference/hg38.UCSC.add_miR.160920.refgene.mat",
        #reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
        markers = str(rules._gistic2_make_markers.output.markers)
    output:
        all_data_by_genes = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/all_data_by_genes.txt",
        all_lesions = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/all_lesions.conf_{conf}.txt",
        all_thresholded_by_genes = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/all_thresholded.by_genes.txt",
        amp_genes = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/amp_genes.conf_{conf}.txt",
        amp_qplot_pdf = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/amp_qplot.pdf",
        amp_qplot_png = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/amp_qplot.png",
        broad_data_by_genes = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/broad_data_by_genes.txt",
        broad_significance_results = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/broad_significance_results.txt",
        broad_values_by_arm = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/broad_values_by_arm.txt",
        del_genes = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/del_genes.conf_{conf}.txt",
        del_qplot_pdf = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/del_qplot.pdf",
        del_qplot_png = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/del_qplot.png",
        focal_data_by_genes = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/focal_data_by_genes.txt",
        freqarms_vs_ngenes = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/freqarms_vs_ngenes.pdf",
        raw_copy_number_pdf = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/raw_copy_number.pdf",
        raw_copy_number_png = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/raw_copy_number.png",
        regions_track = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/regions_track.conf_{conf}.bed",
        sample_cutoffs = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/sample_cutoffs.txt",
        sample_seg_counts = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/sample_seg_counts.txt",
        scores = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/scores.gistic"
    log:
        stdout = CFG["logs"]["gistic2"] + "{seq_type}--projection/all--{projection}/gistic2.stdout.log",
        stderr = CFG["logs"]["gistic2"] + "{seq_type}--projection/all--{projection}/gistic2.stderr.log"
    params:
        base_dir = CFG["dirs"]["gistic2"] + "{seq_type}--projection/all--{projection}/",
        opts = CFG["options"]["gistic2_run"]
    conda:
        CFG["conda_envs"]["gistic2"]
    threads:
        CFG["threads"]["gistic2_run"]
    resources:
        **CFG["resources"]["gistic2_run"]

    shell:
        op.as_one_line("""
        gistic2 -b {params.basedir} -seg {input.seg} -mk {input.markers} -refgene {input.refgene_mat} 
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
        all_data_by_genes = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/all_data_by_genes.txt",
        all_lesions = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/all_lesions.conf_{conf}.txt",
        all_thresholded_by_genes = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/all_thresholded.by_genes.txt",
        amp_genes = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/amp_genes.conf_{conf}.txt",
        amp_qplot_pdf = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/amp_qplot.pdf",
        amp_qplot_png = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/amp_qplot.png",
        broad_data_by_genes = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/broad_data_by_genes.txt",
        broad_significance_results = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/broad_significance_results.txt",
        broad_values_by_arm = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/broad_values_by_arm.txt",
        del_genes = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/del_genes.conf_{conf}.txt",
        del_qplot_pdf = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/del_qplot.pdf",
        del_qplot_png = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/del_qplot.png",
        focal_data_by_genes = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/focal_data_by_genes.txt",
        freqarms_vs_ngenes = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/freqarms_vs_ngenes.pdf",
        raw_copy_number_pdf = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/raw_copy_number.pdf",
        raw_copy_number_png = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/raw_copy_number.png",
        regions_track = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/regions_track.conf_{conf}.bed",
        sample_cutoffs = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/sample_cutoffs.txt",
        sample_seg_counts = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/sample_seg_counts.txt",
        scores = CFG["dirs"]["outputs"] + "{seq_type}--projection/all--{projection}/scores.gistic"

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
            expand(
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
                str(rules._gistic2_output.output.scores),
                zip,
                projection=CFG["projections"],
                seq_type=list(CFG["runs"]["tumour_seq_type"].unique()),
                allow_missing = True # Allows snakemake to expand on these wildcards and fill conf below
            ),
            conf = CFG["options"]["gistic2_run"]["conf_level"] # Specified as a list
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
