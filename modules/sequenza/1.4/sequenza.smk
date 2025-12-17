#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     Bruno Grande


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import pandas as pd

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
# `CFG` is a shortcut to `config["lcr-modules"]["sequenza"]`
CFG = op.setup_module(
    name = "sequenza",
    version = "1.4",
    subdirectories = ["inputs", "seqz", "sequenza", "igv_seg", "convert_coordinates", "fill_regions", "normalize", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _sequenza_input_bam,
    _sequenza_input_chroms,
    _sequenza_input_dbsnp_pos,
    _sequenza_output_seg,
    _sequenza_output_projection,
    _sequenza_output_sub,
    _sequenza_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _sequenza_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Pulls in list of chromosomes for the genome builds
rule _sequenza_input_chroms:
    input:
        txt = ancient(reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt"))
    output:
        txt = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.absolute_symlink(input.txt, output.txt)



rule _sequenza_input_dbsnp_pos:
    input:
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        pos = CFG["dirs"]["inputs"] + "dbsnp/{genome_build}/dbsnp.common_all-151.pos"
    log:
        stderr = CFG["logs"]["inputs"] + "{genome_build}/sequenza_input_dbsnp_pos.stderr.log"
    resources:
        **CFG["resources"]["vcf_sort"]
    shell:
        op.as_one_line("""
        gzip -dc {input.vcf}
            |
        awk 'BEGIN {{FS="\t"}} $0 !~ /^#/ {{print $1 ":" $2}}' 2>> {log.stderr}
            |
        LC_ALL=C sort -S {resources.mem_mb}M > {output.pos} 2>> {log.stderr}
        """)


rule _sequenza_bam2seqz:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        gc_wiggle = reference_files("genomes/{genome_build}/annotations/gc_wiggle.window_50.wig.gz"),
        genome = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        seqz = temp(CFG["dirs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/chromosomes/{chrom}.binned.seqz.gz")
    log:
        stderr = CFG["logs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_bam2seqz.{chrom}.stderr.log"
    params:
        bam2seqz_opts = CFG["options"]["bam2seqz"],
        seqz_binning_opts = CFG["options"]["seqz_binning"]
    conda:
        CFG["conda_envs"]["sequenza-utils"]
    threads:
        CFG["threads"]["bam2seqz"]
    resources:
        **CFG["resources"]["bam2seqz"]
    shell:
        op.as_one_line("""
        sequenza-utils bam2seqz {params.bam2seqz_opts} -gc {input.gc_wiggle} --fasta {input.genome}
        --normal {input.normal_bam} --tumor {input.tumour_bam} --chromosome {wildcards.chrom} 2>> {log.stderr}
            |
        sequenza-utils seqz_binning {params.seqz_binning_opts} --seqz - 2>> {log.stderr}
            |
        gzip > {output} 2>> {log.stderr}
        """)


def _sequenza_request_chrom_seqz_files(wildcards):
    CFG = config["lcr-modules"]["sequenza"]
    with open(checkpoints._sequenza_input_chroms.get(**wildcards).output.txt) as f:
        mains_chroms = f.read().rstrip("\n").split("\n")
    seqz_files = expand(
        CFG["dirs"]["seqz"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/chromosomes/{chrom}.binned.seqz.gz",
        chrom=mains_chroms
    )
    return seqz_files


rule _sequenza_merge_seqz:
    input:
        seqz = _sequenza_request_chrom_seqz_files,
        merge_seqz = CFG["inputs"]["merge_seqz"],
        gc = reference_files("genomes/{genome_build}/annotations/gc_wiggle.window_50.wig.gz")
    output:
        seqz = CFG["dirs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.binned.unfiltered.seqz.gz"
    log:
        stderr = CFG["logs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_merge_seqz.stderr.log"
    threads:
        CFG["threads"]["merge_seqz"]
    resources:
        **CFG["resources"]["merge_seqz"]
    shell:
        op.as_one_line("""
        bash {input.merge_seqz} {input.seqz} 2>> {log.stderr}
            |
        gzip > {output.seqz} 2>> {log.stderr}
        """)


rule _sequenza_filter_seqz:
    input:
        seqz = str(rules._sequenza_merge_seqz.output.seqz),
        filter_seqz = CFG["inputs"]["filter_seqz"],
        dbsnp_pos = str(rules._sequenza_input_dbsnp_pos.output.pos),
        blacklist = reference_files("genomes/{genome_build}/encode/encode-blacklist.{genome_build}.bed")
    output:
        seqz = CFG["dirs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.binned.filtered.seqz.gz"
    log:
        stderr = CFG["logs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_filter_seqz.stderr.log"
    conda:
        CFG["conda_envs"]["bedtools"]
    threads:
        CFG["threads"]["filter_seqz"]
    resources:
        **CFG["resources"]["filter_seqz"]
    shell:
        op.as_one_line("""
        SEQZ_BLACKLIST_BED_FILES='{input.blacklist}'
        {input.filter_seqz} {input.seqz} {input.dbsnp_pos} 2>> {log.stderr}
            |
        gzip > {output.seqz} 2>> {log.stderr}
        """)


rule _sequenza_run:
    input:
        seqz = CFG["dirs"]["seqz"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merged.binned.{filter_status}.seqz.gz",
        run_sequenza = CFG["inputs"]["run_sequenza"],
        assembly = reference_files("genomes/{genome_build}/version.txt"),
        chroms = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt"),
        x_chrom = reference_files("genomes/{genome_build}/genome_fasta/chromosome_x.txt")
    output:
        segments = CFG["dirs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{filter_status}/sequenza_segments.txt"
    log:
        stdout = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_run.{filter_status}.stdout.log",
        stderr = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_run.{filter_status}.stderr.log"
    conda:
        CFG["conda_envs"]["r-sequenza"]
    threads:
        CFG["threads"]["sequenza"]
    resources:
        **CFG["resources"]["sequenza"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {input.run_sequenza} {input.seqz} {input.assembly} {input.chroms} {input.x_chrom}
        $(dirname {output.segments}) {threads} > {log.stdout} 2> {log.stderr}
        """)


# Fill subclones-like output with empty regions for compatibility with downstream tools
rule _sequenza_fill_txt:
    input:
        sub = str(rules._sequenza_run.output.segments)
    output:
        sub = CFG["dirs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{filter_status}/sequenza_segments_filled.txt"
    log:
        stderr = CFG["logs"]["sequenza"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sequenza_fill_txt.{filter_status}.stderr.log"
    threads: 1
    group: "sequenza_post_process"
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"],
        script = "fill_segments.sh",
        arm_file = lambda w: "src/chromArm.hg38.bed" if "38" in str({w.genome_build}) else "src/chromArm.grch37.bed",
        blacklist_file = lambda w: "src/blacklisted.hg38.bed" if "38" in str({w.genome_build}) else "src/blacklisted.grch37.bed"
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.path}{params.script}
        {params.path}{params.arm_file}
        {input.sub}
        {params.path}{params.blacklist_file}
        {output.sub}
        {wildcards.tumour_id}
        sequenza
        2>> {log.stderr}
        """)


rule _sequenza_cnv2igv:
    input:
        segments = str(rules._sequenza_run.output.segments),
        cnv2igv =  ancient(CFG["inputs"]["cnv2igv"])
    output:
        igv = CFG["dirs"]["igv_seg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{filter_status}/sequenza_segments.igv.seg"
    log:
        stderr = CFG["logs"]["igv_seg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cnv2igv.{filter_status}.stderr.log"
    params:
        opts = CFG["options"]["preserve"]
    conda:
        CFG["conda_envs"]["cnv2igv"]
    group: "sequenza_post_process"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id} on $(hostname) at $(date)" > {log.stderr};
        python {input.cnv2igv} --mode sequenza {params.opts} --sample {wildcards.tumour_id}
        {input.segments} > {output.igv} 2>> {log.stderr}
        """)


def _sequenza_get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")


# Convert the coordinates of seg file to a different genome build
rule _sequenza_convert_coordinates:
    input:
        sequenza_native = str(rules._sequenza_cnv2igv.output.igv).replace("{filter_status}", "filtered"),
        sequenza_chain = _sequenza_get_chain
    output:
        sequenza_lifted = CFG["dirs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.stderr.log"
    threads: 1
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        liftover_minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    group: "sequenza_post_process"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.liftover_script}
        SEG
        {input.sequenza_native}
        {output.sequenza_lifted}
        {input.sequenza_chain}
        YES
        {params.liftover_minmatch}
        2>> {log.stderr}
        """)

# ensure to request the correct files for each projection and drop wildcards that won't be used downstream
def _sequenza_prepare_projection(wildcards):
    CFG = config["lcr-modules"]["sequenza"]
    tbl = CFG["runs"]
    this_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"].tolist()

    prefixed_projections = CFG["options"]["prefixed_projections"]
    non_prefixed_projections = CFG["options"]["non_prefixed_projections"]

    if any(substring in this_genome_build[0] for substring in prefixed_projections):
        hg38_projection = str(rules._sequenza_cnv2igv.output.igv).replace("{genome_build}", this_genome_build[0]).replace("{filter_status}", "filtered")
        grch37_projection = str(rules._sequenza_convert_coordinates.output.sequenza_lifted).replace("{genome_build}", this_genome_build[0])
        # handle the hg19 (prefixed) separately
        if "38" in str(this_genome_build[0]):
            grch37_projection = grch37_projection.replace("{chain}", "hg38ToHg19")
        else:
            grch37_projection = grch37_projection.replace("{chain}", "hg19ToHg38")

    elif any(substring in this_genome_build[0] for substring in non_prefixed_projections):
        grch37_projection = str(rules._sequenza_cnv2igv.output.igv).replace("{genome_build}", this_genome_build[0]).replace("{filter_status}", "filtered")
        hg38_projection = str(rules._sequenza_convert_coordinates.output.sequenza_lifted).replace("{genome_build}", this_genome_build[0])
        # handle the grch38 (non-prefixed) separately
        if "38" in str(this_genome_build[0]):
            hg38_projection = hg38_projection.replace("{chain}", "hg38ToHg19")
        else:
            hg38_projection = hg38_projection.replace("{chain}", "hg19ToHg38")
    else:
        raise AttributeError(f"The specified genome build {this_genome_build[0]} is not specified in the config under options to indicate its chr prefixing.")

    return{
        "grch37_projection": grch37_projection,
        "hg38_projection": hg38_projection
    }


# Fill segments of both native and filled file
rule _sequenza_fill_segments:
    input:
        unpack(_sequenza_prepare_projection)
    output:
        grch37_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"),
        hg38_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg")
    log:
        stderr = CFG["logs"]["fill_regions"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}_fill_segments.stderr.log"
    threads: 1
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"]
    conda:
        CFG["conda_envs"]["bedtools"]
    group: "sequenza_post_process"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        echo "Filling grch37 projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.grch37.bed
        {input.grch37_projection}
        {params.path}src/blacklisted.grch37.bed
        {output.grch37_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        echo "Filling hg38 projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.hg38.bed
        {input.hg38_projection}
        {params.path}src/blacklisted.hg38.bed
        {output.hg38_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        """)

def _sequenza_determine_projection(wildcards):
    CFG = config["lcr-modules"]["sequenza"]
    if any(substring in wildcards.projection for substring in ["hg19", "grch37", "hs37d5"]):
        this_file = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"
    elif any(substring in wildcards.projection for substring in ["hg38", "grch38"]):
        this_file = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg"
    return (this_file)


# Normalize chr prefix of the output file
rule _sequenza_normalize_projection:
    input:
        filled = _sequenza_determine_projection,
        chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
    output:
        projection = CFG["dirs"]["normalize"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    resources:
        **CFG["resources"]["post_sequenza"]
    threads: 1
    group: "sequenza_post_process"
    run:
        # read the main chromosomes file of the projection
        chromosomes = pd.read_csv(input.chrom_file, sep = "\t", names=["chromosome"], header=None)
        # handle chr prefix
        if "chr" in chromosomes["chromosome"][0]:
            seg_open = pd.read_csv(input.filled, sep = "\t")
            chrom = list(seg_open['chrom'])
            # avoid cases of chrchr1 if the prefix already there
            for i in range(len(chrom)):
                if 'chr' not in str(chrom[i]):
                    chrom[i]='chr'+str(chrom[i])
            seg_open.loc[:, 'chrom']=chrom
            seg_open.to_csv(output.projection, sep="\t", index=False, na_rep='NA')
        else:
            # remove chr prefix
            seg_open = pd.read_csv(input.filled, sep = "\t")
            seg_open["chrom"] = seg_open["chrom"].astype(str).str.replace('chr', '')
            seg_open.to_csv(output.projection, sep="\t", index=False, na_rep='NA')


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sequenza_output_projection:
    input:
        projection = str(rules._sequenza_normalize_projection.output.projection)
    output:
        projection = CFG["dirs"]["outputs"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    wildcard_constraints:
        projection = "|".join(CFG["requested_projections"]),
        pair_status = "|".join(set(CFG["runs"]["pair_status"].tolist()))
    group: "sequenza_post_process"
    run:
        op.relative_symlink(input.projection, output.projection, in_module = True)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sequenza_output_seg:
    input:
        seg = str(rules._sequenza_cnv2igv.output.igv).replace("{filter_status}", "filtered")
    output:
        seg = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.seg"
    wildcard_constraints:
        projection = "|".join(CFG["requested_projections"]),
        pair_status = "|".join(set(CFG["runs"]["pair_status"].tolist()))
    group: "sequenza_post_process"
    run:
        op.relative_symlink(input.seg, output.seg, in_module=True)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sequenza_output_sub:
    input:
        sub = str(rules._sequenza_fill_txt.output.sub).replace("{filter_status}", "filtered")
    output:
        sub = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.txt"
    group: "sequenza_post_process"
    run:
        op.relative_symlink(input.sub, output.sub, in_module=True)

# Generates the target sentinels for each run, which generate the symlinks
rule _sequenza_all:
    input:
        expand(
            [
                str(rules._sequenza_output_seg.output.seg),
                str(rules._sequenza_output_sub.output.sub),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]),
        expand(
            expand(
            [
                str(rules._sequenza_output_projection.output.projection)
            ],
            zip,  # Run expand() with zip(), not product()
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            seq_type=CFG["runs"]["tumour_seq_type"],
            pair_status=CFG["runs"]["pair_status"],
            allow_missing=True),
            tool="sequenza",
            projection=CFG["requested_projections"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
