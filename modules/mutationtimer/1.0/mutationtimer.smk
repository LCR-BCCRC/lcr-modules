#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Sierra Gillis
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
# `CFG` is a shortcut to `config["lcr-modules"]["mutationtimer"]`
CFG = op.setup_module(
    name = "mutationtimer",
    version = "1.0",
    subdirectories = ["inputs", "convert2bed", "liftover", "mutationtimer", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _mutationtimer_input_bb,
    _mutationtimer_convert2bed,
    _mutationtimer_liftover,
    _mutationtimer_output_tsvs,
    _mutationtimer_all,


##### RULES #####

# Get GAMBLR config into run dir
rule _mutationtimer_gamblr_config:
    params:
        config_url = CFG["inputs"]["gamblr_config_url"],
    output:
        config = "config.yml"
    shell:
        op.as_one_line("""
        wget -qO {output.config} {params.config_url}
        """)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mutationtimer_input_bb:
    input:
        bb = CFG["inputs"]["battenberg_txt"]
    output:
        bb = CFG["dirs"]["inputs"] + "battenberg--{genome_build}/{tumour_id}--{normal_id}_subclones.txt"
    run:
        op.absolute_symlink(input.bb, output.bb)

rule _mutationtimer_input_cellularity:
    input:
        cellularity = CFG["inputs"]["cellularity_txt"]
    output:
        cellularity = CFG["dirs"]["inputs"] + "battenberg--{genome_build}/{tumour_id}--{normal_id}_cellularity_ploidy.txt"
    run:
        op.absolute_symlink(input.cellularity, output.cellularity)


# Subsets subclones.txt input to only necessary cols and outputs as bed for liftover
rule _mutationtimer_convert2bed:
    input:
        bb = str(rules._mutationtimer_input_bb.output.bb)
    output:
        bed = CFG["dirs"]["convert2bed"] + "battenberg--{genome_build}/{tumour_id}--{normal_id}_subclones.bed"
    log:
        stderr = CFG["logs"]["mutationtimer"] + "{tumour_id}--{normal_id}/{genome_build}/convert2bed.stderr.log"
    threads: 1
    shell:
        op.as_one_line("""
        awk -F"\t" -v OFS="\t" '{{print $1,$2,$3,$8,$9,$10,$11,$12,$13}}' {input.bb} >> {output.bed} 2> {log.stderr}
        """)

# Get chain for liftover based on genome build
def get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")

rule _mutationtimer_liftover:
    input:
        bed = str(rules._mutationtimer_convert2bed.output.bed),
        chain = get_chain
    output:
        lifted = CFG["dirs"]["liftover"] + "from--{genome_build}/{tumour_id}--{normal_id}_lifted_{chain}.bed",
        unampped = CFG["dirs"]["liftover"] + "from--{genome_build}/{tumour_id}--{normal_id}_lifted_{chain}.unmapped.bed"
    log:
        stderr = CFG["logs"]["mutationtimer"] + "{tumour_id}--{normal_id}/{genome_build}/liftover_{chain}.stderr.log"
    threads: 1
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    wildcard_constraints:
        chain = "hg38ToHg19|hg19ToHg38"
    shell:
        op.as_one_line("""
        bash {params.liftover_script}  BED {input.bed}
        {output.lifted} {input.chain}
        YES {params.minmatch}
        2> {log.stderr}
        """)

# Ensures the correct subclones file, naive or lifted, is used as input
def _prepare_mt_inputs(wildcards):
    if "19" in wildcards.projection:
        genome_list = ["hg19", "grch37", "hs37d5"]
    else:
        genome_list = ["hg38"]

    CFG = config["lcr-modules"]["mutationtimer"]
    tbl = CFG["runs"]
    tumor_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.normal_sample_id == wildcards.normal_id)]["tumour_genome_build"].tolist()

    # build and projection "match"
    if str(tumor_genome_build[0]) in genome_list:
        cellularity = str(rules._mutationtimer_input_cellularity.output.cellularity).replace("{genome_build}", tumor_genome_build[0])
        bb = str(rules._mutationtimer_convert2bed.output.bed).replace("{genome_build}", tumor_genome_build[0])
    # build and projection differ, and projection is hg19, means build was hg38
    elif "19" in wildcards.projection:
        cellularity = str(rules._mutationtimer_input_cellularity.output.cellularity).replace("{genome_build}", tumor_genome_build[0])
        bb = str(rules._mutationtimer_liftover.output.lifted).replace("{genome_build}", tumor_genome_build[0]).replace("{chain}", "hg38ToHg19")
    # build and projection differ, and projection is not hg19, means build was hg19
    else:
        cellularity = str(rules._mutationtimer_input_cellularity.output.cellularity).replace("{genome_build}", tumor_genome_build[0])
        bb = str(rules._mutationtimer_liftover.output.lifted).replace("{genome_build}", tumor_genome_build[0]).replace("{chain}", "hg19ToHg38")

    return{
        "bb": bb,
        "cellularity": cellularity
    }


rule  _mutationtimer_run:
    input:
        unpack(_prepare_mt_inputs),
        gamblr = ancient(rules._mutationtimer_gamblr_config.output.config)
    output:
        timed_ssm = CFG["dirs"]["mutationtimer"] + "{tumour_id}--{normal_id}/{tumour_id}_timed_ssm.{projection}.tsv",
        timed_cna = CFG["dirs"]["mutationtimer"] + "{tumour_id}--{normal_id}/{tumour_id}_timed_cna.{projection}.tsv"
    log:
        stderr = CFG["logs"]["mutationtimer"] + "{tumour_id}--{normal_id}/{projection}/mutationtimer.stderr.log"
    params:
        script = CFG["options"]["mutationtimer_script"],
        n_bootstrap = CFG["options"]["n_bootstrap"]
    conda:
        CFG["conda_envs"]["mutationtimer"]
    threads:
        CFG["threads"]["mutationtimer"]
    resources:
        **CFG["resources"]["mutationtimer"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {params.script}
        {input.bb}
        {input.cellularity}
        {params.n_bootstrap}
        {output.timed_ssm}
        {output.timed_cna}
        {wildcards.tumour_id}
        {wildcards.projection}
        {log}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mutationtimer_output_tsvs:
    input:
        timed_ssm = str(rules._mutationtimer_run.output.timed_ssm),
        timed_cna = str(rules._mutationtimer_run.output.timed_cna)
    output:
        timed_ssm = CFG["dirs"]["outputs"] + "timed_ssm/{projection}/{tumour_id}--{normal_id}_timed_ssm.{projection}.tsv",
        timed_cna = CFG["dirs"]["outputs"] + "timed_cna/{projection}/{tumour_id}--{normal_id}_timed_cna.{projection}.tsv"
    run:
        op.relative_symlink(input.timed_ssm, output.timed_ssm, in_module= True)
        op.relative_symlink(input.timed_cna, output.timed_cna, in_module= True)



# Generates the target sentinels for each run, which generate the symlinks
rule _mutationtimer_all:
    input:
        expand([str(rules._mutationtimer_liftover.output.lifted)],
            zip,
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            chain=["hg38ToHg19" if "38" in str(x) else "hg19ToHg38" for x in CFG["runs"]["tumour_genome_build"]]
        ),
        expand(
            expand(
            [
                str(rules._mutationtimer_output_tsvs.output.timed_ssm),
                str(rules._mutationtimer_output_tsvs.output.timed_cna)
            ],
            zip,
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            allow_missing=True),
            projection=["hg19", "hg38"]
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
