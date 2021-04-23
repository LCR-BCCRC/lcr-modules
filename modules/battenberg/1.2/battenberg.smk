#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A

##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import glob

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
# `CFG` is a shortcut to `config["lcr-modules"]["battenberg"]`
CFG = op.setup_module(
    name = "battenberg",
    version = "1.1",
    subdirectories = ["inputs", "infer_sex","battenberg", "outputs"],
)

#set variable for prepending to PATH based on config
SCRIPT_PATH = CFG['inputs']['src_dir']
#this is used in place of the shell.prefix() because that was not working consistently. This is not ideal. 

#this preserves the variable when using lambda functions
_battenberg_CFG = CFG

# Define rules to be run locally when using a compute cluster
localrules:
    _battenberg_all

VERSION_MAP = {
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
        battenberg_impute =  directory(CFG["dirs"]["inputs"] + "reference/{genome_build}/battenberg_impute_v3"),
        impute_info = CFG["dirs"]["inputs"] + "reference/{genome_build}/impute_info.txt",
        probloci =  CFG["dirs"]["inputs"] + "reference/{genome_build}/probloci.txt.gz",
        battenberg_wgs_replic_correction = directory(CFG["dirs"]["inputs"] + "reference/{genome_build}/battenberg_wgs_replic_correction_1000g_v3"),
        battenberg_gc_correction = directory(CFG["dirs"]["inputs"] + "reference/{genome_build}/battenberg_wgs_gc_correction_1000g_v3"),  
        genomesloci = directory(CFG["dirs"]["inputs"] + "reference/{genome_build}/battenberg_1000genomesloci2012_v3")
    params:
        url = "https://www.bcgsc.ca/downloads/morinlab/reference",
        alt_build = lambda w: VERSION_MAP[w.genome_build],
        folder = CFG["dirs"]["inputs"] + "reference/{genome_build}",
        build = "{genome_build}",
        PATH = CFG['inputs']['src_dir']
    resources:
        **CFG["resources"]["reference"]
    threads:
        CFG["threads"]["reference"]
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
        python {params.PATH}/reference_correction.py {params.build}
        &&
        wget -qO-  {params.url}/battenberg_{params.alt_build}_replic_correction.tar.gz |
        tar -xvz > {output.battenberg_wgs_replic_correction} -C {params.folder}
        &&
        wget -O {output.probloci} {params.url}/probloci_{params.alt_build}.txt.gz
        
        """)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _battenberg_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
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
        complete = "config/envs/battenberg_dependencies_installed.success"
    conda:
        CFG["conda_envs"]["battenberg"]
    log:
        input = CFG["logs"]["inputs"] + "input.log"
    shell:
        """
        R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")' >> {log.input} && ##move some of this to config?
        R -q -e 'devtools::install_github("morinlab/battenberg")' >> {log.input} &&              ##move some of this to config?
        touch {output.complete}"""

# this process is very fast on bam files and painfully slow on cram files. 
# The result of calc_sex_status.sh is stored in a file to avoid having to rerun it unnecessarily
rule _infer_patient_sex:
    input: 
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output: sex_result = CFG["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex"
    resources:
        **CFG["resources"]["infer_sex"]
    log:
        stderr = CFG["logs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}_infer_sex_stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    group: "setup_run"
    threads: 8
    shell:
        op.as_one_line(""" 
        echo "{params.checker}";
        
        """)


# This rule runs the entire Battenberg pipeline. Eventually we may want to set this rule up to allow re-starting
# of partially completed jobs (e.g. if they run out of RAM and are killed by the cluster, they can automatically retry)
rule _run_battenberg:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        installed = "config/envs/battenberg_dependencies_installed.success",
        sex_result = CFG["dirs"]["infer_sex"] + "{seq_type}--{genome_build}/{normal_id}.sex",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        impute_info = CFG["dirs"]["inputs"] + "reference/{genome_build}/impute_info.txt"

    output:
        refit=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_refit_suggestion.txt",
        sub=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.txt",
        ac=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_alleleCounts.tab"),
        mb=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantBAF.tab"),
        mlrg=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR_gcCorrected.tab"),
        mlr=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_mutantLogR.tab"),
        nlr=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalLogR.tab"),
        nb=temp(CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_normalBAF.tab"),
        cp=CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cellularity_ploidy.txt"
    log:
        stdout = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stdout.log",
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_battenberg.stderr.log"
    params:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        script = CFG["inputs"]["battenberg_script"],
        out_dir = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        ref = CFG["dirs"]["inputs"] + "reference/{genome_build}"
    conda:
        CFG["conda_envs"]["battenberg"]
    resources:
        **CFG["resources"]["battenberg"]
    threads:
        CFG["threads"]["battenberg"]
    shell:
       op.as_one_line("""
        if [[ $(head -c 4 {params.fasta}) == ">chr" ]]; then chr_prefixed='true'; else chr_prefixed='false'; fi;
        echo "$chr_prefixed"
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stdout};
        sex=$(cut -f 4 {input.sex_result}| tail -n 1); 
        echo "setting sex as $sex";
        Rscript {params.script} -t {wildcards.tumour_id} 
        -n {wildcards.normal_id} --tb {input.tumour_bam} --nb {input.normal_bam} -f {input.fasta} --reference {params.ref}
        -o {params.out_dir} --chr_prefixed_genome $chr_prefixed --sex $sex --cpu {threads} >> {log.stdout} 2>> {log.stderr} &&  
        echo "DONE {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" >> {log.stdout}; 
        """)


# Convert the subclones.txt (best fit) to igv-friendly SEG files. 
rule _battenberg_to_igv_seg:
    input:
        sub = rules._run_battenberg.output.sub,
        cnv2igv = CFG["inputs"]["cnv2igv"]
    output:
        seg = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_subclones.igv.seg"
    log:
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_seg2igv.stderr.log"
    threads: 1
    group: "post_process"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        python {input.cnv2igv} --mode battenberg --sample {wildcards.tumour_id} 
        {input.sub} > {output.seg} 2>> {log.stderr}
        """)


#due to the large number of files (several per chromosome) that are not explicit outputs, do some glob-based cleaning in the output directory
rule _battenberg_cleanup:
    input:
        rules._battenberg_to_igv_seg.output.seg
    output:
        complete = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}/{tumour_id}_cleanup_complete.txt"
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

# Symlinks the final output files into the module results directory (under '99-outputs/')
# All plots generated by Battenberg are symlinked using a glob for convenience

rule _battenberg_output_seg:
    input:
        seg = rules._battenberg_to_igv_seg.output.seg,
        sub = rules._run_battenberg.output.sub,
        cp = rules._run_battenberg.output.cp
    output:
        seg = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_subclones.igv.seg",
        sub = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_subclones.txt",
        cp = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_cellularity_ploidy.txt"
    params: 
        batt_dir = CFG["dirs"]["battenberg"] + "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}",
        png_dir = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}"
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
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])





##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
