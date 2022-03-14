#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A

##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import glob
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

# Setup module and store module-specific configuration in `CFG_battenberg`
# `CFG_battenberg` is a shortcut to `config["lcr-modules"]["battenberg"]`
CFG_battenberg = op.setup_module(
    name = "battenberg",
    version = "1.2",
    subdirectories = ["inputs", "infer_sex","battenberg", "normalize", "outputs"],
)

#set variable for prepending to PATH based on config
BATTENBERG_SCRIPT_PATH = CFG_battenberg['inputs']['src_dir']
#this is used in place of the shell.prefix() because that was not working consistently. This is not ideal. 

#this preserves the variable when using lambda functions
#_battenberg_CFG_battenberg = CFG_battenberg

# Define rules to be run locally when using a compute cluster
localrules:
    _battenberg_all

BATTENBERG_VERSION_MAP = {
    "hg19": "grch37",
    "grch37": "grch37",
    "hs37d5": "grch37",
    "hg38": "hg38",
    "grch38": "hg38",
    "grch38-legacy": "hg38"

}

##### RULES #####

# Downloads the reference files into the module results directory (under '00-inputs/') from https://www.bcgsc.ca/downloads/morinlab/reference/ . 
rule _battenberg_get_reference:
    output:
        battenberg_impute =  directory(CFG_battenberg["dirs"]["inputs"] + "reference/{genome_build}/battenberg_impute_v3"),
        impute_info = CFG_battenberg["dirs"]["inputs"] + "reference/{genome_build}/impute_info.txt",
        probloci =  CFG_battenberg["dirs"]["inputs"] + "reference/{genome_build}/probloci.txt.gz",
        battenberg_wgs_replic_correction = directory(CFG_battenberg["dirs"]["inputs"] + "reference/{genome_build}/battenberg_wgs_replic_correction_1000g_v3"),
        battenberg_gc_correction = directory(CFG_battenberg["dirs"]["inputs"] + "reference/{genome_build}/battenberg_wgs_gc_correction_1000g_v3"),  
        genomesloci = directory(CFG_battenberg["dirs"]["inputs"] + "reference/{genome_build}/battenberg_1000genomesloci2012_v3")
    params:
        url = "https://www.bcgsc.ca/downloads/morinlab/reference",
        alt_build = lambda w: BATTENBERG_VERSION_MAP[w.genome_build],
        folder = CFG_battenberg["dirs"]["inputs"] + "reference/{genome_build}",
        build = "{genome_build}",
        battenberg_path = CFG_battenberg['inputs']['src_dir']
    resources:
        **CFG_battenberg["resources"]["reference"]
    threads:
        CFG_battenberg["threads"]["reference"]
    shell:
        op.as_one_line("""
        wget -qO-  {params.url}/battenberg_impute_{params.alt_build}.tar.gz  |
        tar -xvz > {output.battenberg_impute} -C {params.folder}
        &&
        wget -qO- {params.url}/battenberg_{params.alt_build}_gc_correction.tar.gz  |
        tar -xvz > {output.battenberg_gc_correction} -C {params.folder}
        &&
        wget -qO- {params.url}/battenberg_1000genomesloci_{params.alt_build}.tar.gz | 
        tar -xvz > {output.genomesloci} -C {params.folder}
        &&
        wget -O {output.impute_info} {params.url}/impute_info_{params.alt_build}.txt
        &&
        python {params.battenberg_path}/reference_correction.py {params.build} $(dirname $(readlink -f {output.impute_info}))
        &&
        wget -qO-  {params.url}/battenberg_{params.alt_build}_replic_correction.tar.gz |
        tar -xvz > {output.battenberg_wgs_replic_correction} -C {params.folder}
        &&
        wget -O {output.probloci} {params.url}/probloci_{params.alt_build}.txt.gz
        
        """)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _battenberg_input_bam:
    input:
        bam = CFG_battenberg["inputs"]["sample_bam"]
    output:
        bam = CFG_battenberg["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG_battenberg["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG_battenberg["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    group: "setup_run"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bam + ".bai", output.bai)
        op.absolute_symlink(input.bam + ".bai", output.crai)

# Installs the Battenberg R dependencies and associated software (impute2, alleleCounter)
# Currently I think this rule has to be run twice for it to work properly because the conda environment is created here. 
# I am open to suggestions for how to get around this.
rule _install_battenberg:
    output:
        complete = CFG_battenberg["dirs"]["inputs"] + "battenberg_dependencies_installed.success"
    conda:
        CFG_battenberg["conda_envs"]["battenberg"]
    log:
        input = CFG_battenberg["logs"]["inputs"] + "input.log"
    shell:
        """
        R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")' >> {log.input} && ##move some of this to config?
        R -q -e 'devtools::install_github("morinlab/battenberg")' >> {log.input} &&              ##move some of this to config?
        touch {output.complete}"""

# this process is very fast on bam files and painfully slow on cram files. 
# The result of calc_sex_status.sh is stored in a file to avoid having to rerun it unnecessarily
rule _infer_patient_sex:
    input: 
        normal_bam = CFG_battenberg["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: sex_result = CFG_battenberg["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex"
    resources:
        **CFG_battenberg["resources"]["infer_sex"]
    log:
        stderr = CFG_battenberg["logs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}_infer_sex_stderr.log"
    conda:
        CFG_battenberg["conda_envs"]["samtools"]
    group: "setup_run"
    threads: 8
    shell:
        op.as_one_line(""" 
        PATH={BATTENBERG_SCRIPT_PATH}:$PATH;
        echo "running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr} ;
        calc_sex_status.sh {input.normal_bam} {input.fasta} {wildcards.normal_id} > {output.sex_result} 2>> {log.stderr} &&
        echo "DONE running {rule} for {wildcards.normal_id} on $(hostname) at $(date)" >> {log.stderr} 
        """)


# This rule runs the entire Battenberg pipeline. Eventually we may want to set this rule up to allow re-starting
# of partially completed jobs (e.g. if they run out of RAM and are killed by the cluster, they can automatically retry)
rule _run_battenberg:
    input:
        tumour_bam = CFG_battenberg["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG_battenberg["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        installed = CFG_battenberg["dirs"]["inputs"] + "battenberg_dependencies_installed.success",
        sex_result = CFG_battenberg["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        impute_info = str(rules._battenberg_get_reference.output.impute_info)

    output:
        refit=CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_refit_suggestion.txt",
        sub=CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.txt",
        ac=temp(CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleCounts.tab"),
        mb=temp(CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantBAF.tab"),
        mlrg=temp(CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR_gcCorrected.tab"),
        mlr=temp(CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR.tab"),
        nlr=temp(CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalLogR.tab"),
        nb=temp(CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalBAF.tab"),
        cp=CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cellularity_ploidy.txt"
    log:
        stdout = CFG_battenberg["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stdout.log",
        stderr = CFG_battenberg["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stderr.log"
    params:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        script = CFG_battenberg["inputs"]["battenberg_script"],
        out_dir = CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        ref = CFG_battenberg["dirs"]["inputs"] + "reference/{genome_build}"
    conda:
        CFG_battenberg["conda_envs"]["battenberg"]
    resources:
        **CFG_battenberg["resources"]["battenberg"]
    threads:
        CFG_battenberg["threads"]["battenberg"]
    shell:
       op.as_one_line("""
        if [[ $(head -c 4 {params.fasta}) == ">chr" ]]; then chr_prefixed='true'; else chr_prefixed='false'; fi;
        echo "$chr_prefixed"
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stdout};
        sex=$(cut -f 4 {input.sex_result}| tail -n 1); 
        echo "setting sex as $sex";
        Rscript {params.script} -t {wildcards.tumour_id} 
        -n {wildcards.normal_id} --tb $(readlink -f {input.tumour_bam}) --nb $(readlink -f {input.normal_bam}) -f {input.fasta} --reference $(readlink -f {params.ref})
        -o {params.out_dir} --chr_prefixed_genome $chr_prefixed --sex $sex --cpu {threads} >> {log.stdout} 2>> {log.stderr} &&  
        echo "DONE {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" >> {log.stdout}; 
        """)


# Convert the subclones.txt (best fit) to igv-friendly SEG files. 
rule _battenberg_to_igv_seg:
    input:
        sub = rules._run_battenberg.output.sub,
        cnv2igv = CFG_battenberg["inputs"]["cnv2igv"]
    output:
        seg = CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.igv.seg"
    log:
        stderr = CFG_battenberg["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_seg2igv.stderr.log"
    threads: 1
    group: "post_process"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        python {input.cnv2igv} --mode battenberg --sample {wildcards.tumour_id} 
        {input.sub} > {output.seg} 2>> {log.stderr}
        """)


# Fill subclones.txt with empty regions for compatibility with downstream tools
rule _battenberg_fill_subclones:
    input:
        sub = str(rules._run_battenberg.output.sub)
    output:
        sub = CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.filled.txt"
    log:
        stderr = CFG_battenberg["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_fill_subclones.stderr.log"
    threads: 1
    group: "post_process"
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/1.0/",
        script = "fill_segments.sh",
        arm_file = lambda w: "src/chromArm.hg38.bed" if "38" in str({w.genome_build}) else "src/chromArm.grch37.bed",
        blacklist_file = lambda w: "src/blacklisted.hg38.bed" if "38" in str({w.genome_build}) else "src/blacklisted.grch37.bed"
    conda:
        CFG_battenberg["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.path}{params.script}
        {params.path}{params.arm_file}
        {input.sub}
        {params.path}{params.blacklist_file}
        {output.sub}
        {wildcards.tumour_id}
        subclones
        2>> {log.stderr}
        """)


#due to the large number of files (several per chromosome) that are not explicit outputs, do some glob-based cleaning in the output directory
rule _battenberg_cleanup:
    input:
        rules._battenberg_to_igv_seg.output.seg
    output:
        complete = CFG_battenberg["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cleanup_complete.txt"
    group: "post_process"
    shell:
        op.as_one_line("""
        d=$(dirname {output});
        rm -f $d/*impute_input* &&
        rm -f $d/*alleleFrequencies* &&
        rm -f $d/*aplotype* &&
        rm -f $d/*BAFsegmented* && 
        touch {output.complete}
        """)


configfile: config["lcr-modules"]["_shared"]["lcr-modules"] + "modules/liftover/" + CFG_battenberg["options"]["liftover_module_version"] + "/config/default.yaml"
config["lcr-modules"]["liftover"]["runs"] = CFG_battenberg["runs"]
config["lcr-modules"]["liftover"]["dirs"]["_parent"] = (config["lcr-modules"]["_shared"]["root_output_dir"] +
    "battenberg" +
    "-" +
    config["lcr-modules"]["battenberg"]["version"] +
    "_liftover-" + CFG_battenberg["options"]["liftover_module_version"])
config["lcr-modules"]["liftover"]["tool"] = "battenberg"
config["lcr-modules"]["liftover"]["input_type"] = "seg"
config["lcr-modules"]["liftover"]["inputs"]["sample_file"] = str(rules._battenberg_to_igv_seg.output.seg)

include: "../../liftover/" + CFG_battenberg["options"]["liftover_module_version"] + "/liftover.smk"

def _battenberg_determine_projection(wildcards):
    CFG_battenberg = config["lcr-modules"]["battenberg"]
    tbl = CFG_battenberg["runs"]
    this_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"]
    
    prefixed_projections = CFG_battenberg["options"]["prefixed_projections"]
    non_prefixed_projections = CFG_battenberg["options"]["non_prefixed_projections"]

    if any(substring in this_genome_build[0] for substring in prefixed_projections):
        print(f"This genome build {this_genome_build[0]} for sample {wildcards.tumour_id} is specified to be chr prefixed.")
        hg38_projection = str(rules._liftover_native_output.output).replace("{genome_build}", this_genome_build[0])
        grch37_projection = str(rules._liftover_output.output).replace("{genome_build}", this_genome_build[0]).replace("{chain}", "hg38ToHg19")

    elif any(substring in this_genome_build[0] for substring in non_prefixed_projections):
        print(f"This genome build {this_genome_build[0]} for sample {wildcards.tumour_id} is specified to be non-prefixed.")
        grch37_projection = str(rules._liftover_native_output.output).replace("{genome_build}", this_genome_build[0])
        hg38_projection = str(rules._liftover_output.output).replace("{genome_build}", this_genome_build[0]).replace("{chain}", "hg19ToHg38")
    
    else:
        raise AttributeError(f"The specified genome build {this_genome_build[0]} is not specified in the config under options to indicate its chr prefixing.")

    return{
        "grch37_projection": grch37_projection,
        "hg38_projection": hg38_projection
    }


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _battenberg_normalize_projection:
    input:
        unpack(_battenberg_determine_projection)
    output:
        grch37_projection = CFG_battenberg["dirs"]["normalize"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg",
        hg38_projection = CFG_battenberg["dirs"]["normalize"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg"
    threads: 1
    group: "post_process"
    run:
        # grch37 is always non-prefixed
        seg_open = pd.read_csv(input.grch37_projection, sep = "\t")
        seg_open["chrom"] = seg_open["chrom"].astype(str).str.replace('chr', '')
        seg_open.to_csv(output.grch37_projection, sep="\t", index=False)
        # hg38 is always prefixed, but add it if currently missing
        seg_open = pd.read_csv(input.hg38_projection, sep = "\t")
        chrom = list(seg_open['chrom'])
        for i in range(len(chrom)):
            if 'chr' not in str(chrom[i]):
                chrom[i]='chr'+str(chrom[i])
        seg_open.loc[:, 'chrom']=chrom
        seg_open.to_csv(output.hg38_projection, sep="\t", index=False)

rule _battenberg_output_projection:
    input:
        grch37_projection = str(rules._battenberg_normalize_projection.output.grch37_projection),
        hg38_projection = str(rules._battenberg_normalize_projection.output.hg38_projection)
    output:
        grch37_projection = CFG_battenberg["dirs"]["outputs"] + CFG_battenberg["output"]["seg"]["grch37_projection"],
        hg38_projection = CFG_battenberg["dirs"]["outputs"] + CFG_battenberg["output"]["seg"]["hg38_projection"]
    threads: 1
    group: "post_process"
    run:
        op.relative_symlink(input.grch37_projection, output.grch37_projection, in_module = True)
        op.relative_symlink(input.hg38_projection, output.hg38_projection, in_module = True)

# Symlinks the final output files into the module results directory (under '99-outputs/')
# All plots generated by Battenberg are symlinked using a glob for convenience

rule _battenberg_output_seg:
    input:
        seg = rules._battenberg_to_igv_seg.output.seg,
        sub = rules._battenberg_fill_subclones.output.sub,
        cp = rules._run_battenberg.output.cp
    output:
        seg = CFG_battenberg["dirs"]["outputs"] + CFG_battenberg["output"]["seg"]["original"],
        sub = CFG_battenberg["dirs"]["outputs"] + CFG_battenberg["output"]["txt"]["subclones"],
        cp = CFG_battenberg["dirs"]["outputs"] + CFG_battenberg["output"]["txt"]["cell_ploidy"]
    params: 
        batt_dir = CFG_battenberg["dirs"]["battenberg"] + "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        png_dir = CFG_battenberg["dirs"]["outputs"] + "png/{seq_type}--{genome_build}"
    group: "post_process"
    run:
        plots = glob.glob(params.batt_dir + "/*.png")
        for png in plots:
            bn = os.path.basename(png)
            op.relative_symlink(png, params.png_dir + "/" + bn,in_module=True)
        op.relative_symlink(input.seg, output.seg,in_module=True)
        op.relative_symlink(input.sub, output.sub,in_module=True)
        op.relative_symlink(input.cp, output.cp,in_module=True)

# Generates the target sentinels for each run, which generate the symlinks
rule _battenberg_all:
    input:
        expand(
            [ 
                rules._run_battenberg.output.sub,
                rules._battenberg_output_seg.output.seg,
                rules._battenberg_cleanup.output.complete
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG_battenberg["runs"]["tumour_seq_type"],
            genome_build=CFG_battenberg["runs"]["tumour_genome_build"],
            tumour_id=CFG_battenberg["runs"]["tumour_sample_id"],
            normal_id=CFG_battenberg["runs"]["normal_sample_id"],
            pair_status=CFG_battenberg["runs"]["pair_status"]),
        expand(
            [
                str(rules._battenberg_output_projection.output.grch37_projection),
                str(rules._battenberg_output_projection.output.hg38_projection)
            ],
            zip,  # Run expand() with zip(), not product()
            tumour_id=CFG_battenberg["runs"]["tumour_sample_id"],
            normal_id=CFG_battenberg["runs"]["normal_sample_id"],
            #genome_build = CFG_battenberg["runs"]["tumour_genome_build"],
            seq_type=CFG_battenberg["runs"]["tumour_seq_type"],
            pair_status=CFG_battenberg["runs"]["pair_status"],
            #repeat the tool name N times in expand so each pair in run is used
            tool=[config["lcr-modules"]["liftover"]["tool"]] * len(CFG_battenberg["runs"]["tumour_sample_id"])
            #chain=["hg38ToHg19" if "38" in str(x) else "hg19ToHg38" for x in CFG_battenberg["runs"]["tumour_genome_build"]]
            )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG_battenberg` variable
CFG=CFG_battenberg
op.cleanup_module(CFG_battenberg)