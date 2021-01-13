#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package oncopipe and inspect packages. 
import oncopipe as op
import inspect

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

# Use modified setup_module function to skip the global `CFG` check. This module assigns the module config to a variable called `CFG_SLMS3` so that the submodules can still use `CFG`. 
def setup_module(name, version, subdirectories):

    # Get namespace where module is being set up
    module_frame = inspect.currentframe().f_back
    module_globals = module_frame.f_globals
    config = module_globals["config"]

    # Ensure minimum version of Snakemake
    snakemake.utils.min_version("5.4.0")

    # Ensure that the lcr-modules _shared config is loaded
    assert "lcr-modules" in config and "_shared" in config["lcr-modules"], (
        "Shared lcr-modules configuration is not loaded. "
        "See README.md in lcr-modules for more information."
    )

    # Ensure that this module's config is loaded
    assert name in config["lcr-modules"], (
        f"The configuration for the {name!r} module is not loaded. "
        "It should be loaded before the module Snakefile (.smk) is "
        "included. See README.md for more information."
    )

    # Get configuration for the given module and create samples shorthand
    mconfig = copy.deepcopy(config["lcr-modules"]["_shared"])
    snakemake.utils.update_config(mconfig, config["lcr-modules"][name])
    msamples = mconfig["samples"].copy()

    # Check whether there are "None" strings
    op.check_for_none_strings(mconfig, name)

    # Check whether there are "__UPDATE__" strings
    op.check_for_update_strings(mconfig, name)

    # Drop samples whose seq_types do not appear in pairing_config
    assert "pairing_config" in mconfig, "`pairing_config` missing from module config."
    sample_seq_types = msamples["seq_type"].unique()
    pairing_config = mconfig["pairing_config"]
    supported_seq_types = [
        k for k, v in pairing_config.items() if "run_paired_tumours" in v
    ]
    unsupported_seq_types = set(sample_seq_types) - set(supported_seq_types)
    if len(unsupported_seq_types) > 0:
        logger.warning(
            f"Some samples have seq_types {unsupported_seq_types} that are "
            f"not configured in the pairing config for the {name} module. "
            "They will be excluded from the analysis."
        )
    msamples = msamples[msamples["seq_type"].isin(supported_seq_types)].copy()
    mconfig["samples"] = msamples

    # Set module name and version
    mconfig["name"] = name
    mconfig["version"] = version

    # Ensure that common module sub-fields are present
    subfields = ["inputs", "dirs", "conda_envs", "options", "threads", "mem_mb"]
    for subfield in subfields:
        if subfield not in mconfig:
            mconfig[subfield] = dict()

    # Check reference
    assert (
        "genome_build" in msamples
    ), "Add a `genome_build` column to your samples data frame."

    # Update placeholders in any string in the module-specific config
    def update_placeholders(obj, **placeholders):
        if isinstance(obj, str):
            result = obj
            for placeholder, value in placeholders.items():
                result = result.replace("{" + placeholder + "}", value)
        else:
            result = obj
        return result

    # Find repository and module directories
    repodir = os.path.normpath(mconfig["lcr-modules"])
    modsdir = os.path.join(repodir, "modules", name, version)
    scriptsdir = os.path.normpath(mconfig["lcr-scripts"])

    placeholders = {
        "REPODIR": repodir,
        "MODSDIR": modsdir,
        "SCRIPTSDIR": scriptsdir,
    }
    mconfig = op.walk_through_dict(mconfig, update_placeholders, **placeholders)

    # Validate samples data frame
    schemas_dir = os.path.join(modsdir, "schemas")
    schemas = os.listdir(schemas_dir)
    for schema in schemas:
        snakemake.utils.validate(msamples, schema=os.path.join(schemas_dir, schema))

    # Configure output directory if not specified and create it
    if mconfig["dirs"].get("_parent") is None:
        root_output_dir = mconfig.get("root_output_dir", "results")
        output_dir = os.path.join(root_output_dir, f"{name}-{version}")
        mconfig["dirs"]["_parent"] = output_dir
    mconfig["dirs"]["_parent"] = mconfig["dirs"]["_parent"].rstrip("/") + "/"
    os.makedirs(mconfig["dirs"]["_parent"], exist_ok=True)

    # Update paths to conda environments to be relative to the module directory
    for env_name, env_val in mconfig["conda_envs"].items():
        if env_val is not None:
            mconfig["conda_envs"][env_name] = os.path.realpath(env_val)

    # Setup output sub-directories
    scratch_subdirs = mconfig.get("scratch_subdirectories", [])
    mconfig = op.setup_subdirs(mconfig, subdirectories, scratch_subdirs)

    # Setup log sub-directories
    mconfig["logs"] = dict()
    parent_dir = mconfig["dirs"]["_parent"]
    launched_fmt = op._session.launched_fmt
    logs_parent_dir = os.path.join(parent_dir, "logs", launched_fmt)
    logs_parent_dir = logs_parent_dir.rstrip("/") + "/"
    mconfig["logs"]["_parent"] = logs_parent_dir
    os.makedirs(logs_parent_dir, exist_ok=True)
    for subdir, value in mconfig["dirs"].items():
        if subdir == "_parent":
            continue
        mconfig["logs"][subdir] = value.replace(parent_dir, logs_parent_dir)

    # Generate runs
    assert "pairing_config" in mconfig, "Module config must have 'pairing_config'."
    runs = op.generate_runs(
        msamples, mconfig["pairing_config"], mconfig.get("unmatched_normal_ids")
    )

    # Split runs based on pair_status
    mconfig["runs"] = runs
    mconfig["paired_runs"] = runs[runs.pair_status != "no_normal"]
    mconfig["unpaired_runs"] = runs[runs.pair_status == "no_normal"]

    # Return module-specific configuration
    config["lcr-modules"][name] = mconfig
    return mconfig


# Setup module and store module-specific configuration in `CFG_SLMS3`
# `CFG_SLMS3` is a shortcut to `config["lcr-modules"]["slms_3"]`
# This is not assigned to `CFG` as in other modules because the sub-modules here will use `CFG`. 
CFG_SLMS3 = setup_module(
    name = "slms_3",
    version = "1.0",
    subdirectories = ["inputs", "strelka_gnomad", "lofreq_gnomad", "strelka_lofreq_union", "mutect2_depth_filt", "outputs"],
)

# Ensure each of the submodule versions meets requirements. 

assert float(CFG_SLMS3["module_versions"]["manta"]) >= 2.0, (
    f"The current manta module version is {CFG_SLMS3['module_versions']['manta']}. "
    "SLMS-3 requires manta module version 2.0 or higher. "
)

assert float(CFG_SLMS3["module_versions"]["starfish"]) >= 2.0, (
    f"The current starfish module version is {CFG_SLMS3['module_versions']['starfish']}. "
    "SLMS-3 requires starfish module version 2.0 or higher. "
)

assert float(CFG_SLMS3["module_versions"]["mutect2"]) >= 2.0, (
    f"The current mutect2 module version is {CFG_SLMS3['module_versions']['mutect2']}. "
    "SLMS-3 requires mutect2 module version 2.0 or higher. "
)

assert float(CFG_SLMS3["module_versions"]["strelka"]) >= 1.1, (
    f"The current strelka module version is {CFG_SLMS3['module_versions']['strelka']}. "
    "SLMS-3 requires strelka module version 1.1 or higher. "
)

# Define rules to be run locally when using a compute cluster
localrules:
    _slms_3_input_sage_vcf,
    _slms_3_input_strelka_vcf,
    _slms_3_input_mutect_vcf, 
    _slms_3_input_lofreq_vcf, 
    _slms_3_mutect2_samples_table, 
    _slms_3_output_vcf, 
    _slms_3_all,

##### MODULE CONFIGFILES #####

# Load default module configs and update with values from slms_3 config

for tool in CFG_SLMS3["module_versions"].keys():
    lcr_modules = config["lcr-modules"]["_shared"]["lcr-modules"] + "modules/"
    configfile: lcr_modules + "/" +  tool + "/" + CFG_SLMS3["module_versions"][tool] + "/config/default.yaml"
    snakemake.utils.update_config(
        config["lcr-modules"][tool], 
        config["lcr-modules"]["slms_3"][tool]
    )
    config["lcr-modules"][tool]["inputs"]["sample_bam"] = CFG_SLMS3["inputs"]["sample_bam"]
    config["lcr-modules"][tool]["inputs"]["sample_bai"] = CFG_SLMS3["inputs"]["sample_bai"]    


##### FIRST PASS VARIANT CALLING MODULE SNAKEFILES #####

lcr_modules = config["lcr-modules"]["_shared"]["lcr-modules"] + "modules/"

# Load Manta first
include: lcr_modules + "/manta/" + CFG_SLMS3["module_versions"]["manta"] + "/manta.smk"

# Update Strelka config to use Manta Candidate Small Indels    
config["lcr-modules"]["strelka"]["inputs"]["candidate_small_indels"] = expand(str(rules._manta_output_vcf.output.vcf), vcf_name = "candidateSmallIndels", allow_missing=True)

# Load Strelka, SAGE, and Lofreq
include: lcr_modules + "/strelka/" + CFG_SLMS3["module_versions"]["strelka"] + "/strelka.smk"
include: lcr_modules + "/sage/" + CFG_SLMS3["module_versions"]["sage"] + "/sage.smk"
include: lcr_modules + "/lofreq/" + CFG_SLMS3["module_versions"]["lofreq"] + "/lofreq.smk"
 
##### FIRST PASS VARIANT CALLING RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')

rule _slms_3_input_strelka_vcf: 
    input:
        vcf = str(rules._strelka_output_filtered_vcf.output.vcf), 
        tbi = str(rules._strelka_output_filtered_vcf.output.vcf_tbi)
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.tbi, output.tbi)

rule _slms_3_input_sage_vcf: 
    input:
        vcf = str(rules._sage_output_vcf.output.combined), 
        tbi = str(rules._sage_output_vcf.output.combined_tbi)
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.tbi, output.tbi)

rule _slms_3_input_lofreq_vcf: 
    input:
        vcf = str(rules._lofreq_output_vcf.output.vcf_all), 
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lofreq.snvs.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lofreq.snvs.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.vcf + ".tbi", output.tbi)

# Annotate Strelka VCF and remove common GnomAD variants

rule _slms_3_annotate_strelka_gnomad:
    input:
        vcf = str(rules._slms_3_input_strelka_vcf.output.vcf),
        tbi = str(rules._slms_3_input_strelka_vcf.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = CFG_SLMS3["dirs"]["strelka_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.combined.gnomad.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["strelka_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka.combined.gnomad.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["strelka_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_gnomad.stderr.log"
    conda:
        CFG_SLMS3["conda_envs"]["bcftools"]
    threads:
        CFG_SLMS3["threads"]["strelka_gnomad"]
    resources:
        **CFG_SLMS3["resources"]["strelka_gnomad"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads} 
        -a {input.gnomad} -c INFO/AF {input.vcf} | 
        awk 'BEGIN {{FS=OFS="\\t"}} {{ if ($1 !~ /^#/ && $8 !~ ";AF=") $8=$8";AF=0"; print $0; }}' | 
        bcftools view -i 'INFO/AF < 0.0001 && INFO/SomaticEVS >= 10 && FMT/DP[1] >= 10' -Oz -o {output.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


# Annotate LoFreq VCF and remove common GnomAD variants

rule _slms_3_annotate_lofreq_gnomad:
    input: 
        vcf = str(rules._slms_3_input_lofreq_vcf.output.vcf), 
        tbi = str(rules._slms_3_input_lofreq_vcf.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output: 
        vcf = CFG_SLMS3["dirs"]["lofreq_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.snvs.gnomad.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["lofreq_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.snvs.gnomad.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["lofreq_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_gnomad.stderr.log"
    conda: 
        CFG_SLMS3["conda_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["lofreq_gnomad"]
    threads: 
        CFG_SLMS3["threads"]["lofreq_gnomad"]
    shell: 
        op.as_one_line("""
        bcftools annotate --threads {threads} 
        -a {input.gnomad} -c INFO/AF {input.vcf} | 
        awk 'BEGIN {{FS=OFS="\\t"}} {{ if ($1 !~ /^#/ && $8 !~ ";AF=") $8=$8";AF=0"; print $0; }}' | 
        bcftools view -i 'INFO/AF < 0.0001' -Oz -o {output.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

# Create a union of Strelka and Lofreq as candidates for Mutect2

rule _slms_3_strelka_lofreq_union: 
    input: 
        strelka = str(rules._slms_3_annotate_strelka_gnomad.output.vcf),
        strelka_tbi = str(rules._slms_3_annotate_strelka_gnomad.output.tbi), 
        lofreq = str(rules._slms_3_annotate_lofreq_gnomad.output.vcf), 
        lofreq_tbi = str(rules._slms_3_annotate_lofreq_gnomad.output.tbi)
    output: 
        vcf = CFG_SLMS3["dirs"]["strelka_lofreq_union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka_lofreq_union_gnomad.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["strelka_lofreq_union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.strelka_lofreq_union_gnomad.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["strelka_lofreq_union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_lofreq_union_gnomad.stderr.log"
    conda: 
        CFG_SLMS3["conda_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["strelka_lofreq_union"]
    threads: 
        CFG_SLMS3["threads"]["strelka_lofreq_union"]
    shell: 
        op.as_one_line("""
        bcftools merge --force-samples -m all --threads {threads} 
        -Oz -o {output.vcf} {input.strelka} {input.lofreq} 2> {log.stderr} && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

##### SECOND PASS VARIANT CALLING MODULE SNAKEFILES #####

# Upate the Mutect2 config with the _strelka_lofreq_vcf

config["lcr-modules"]["mutect2"]["inputs"]["candidate_positions"] = str(rules._slms_3_strelka_lofreq_union.output.vcf)

include: lcr_modules + "/mutect2/" + CFG_SLMS3["module_versions"]["mutect2"] + "/mutect2.smk"

##### SECOND PASS VARIANT CALLING RULES #####

rule _slms_3_input_mutect_vcf: 
    input:
        vcf = str(rules._mutect2_output_vcf.output.vcf), 
        tbi = str(rules._mutect2_output_vcf.output.tbi)
    output:
        vcf = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.mutect2.combined.vcf", 
        tbi = CFG_SLMS3["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.mutect2.combined.vcf.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.tbi, output.tbi)


# Depth/VAF filter the Mutect2 output. This must be completed here and not in the mutect2 module so that the TUMOR column can be properly indexed. 
rule _slms_3_mutect2_samples_table: 
    output: 
        table = CFG_SLMS3["dirs"]["mutect2_depth_filt"] + "samples.txt"
    shell: 
        """echo "TUMOR" > {output.table}"""

rule _slms_3_mutect2_depth_filt: 
    input: 
        vcf = str(rules._slms_3_input_mutect_vcf.output.vcf), 
        table = str(rules._slms_3_mutect2_samples_table.output.table)
    output: 
        vcf = CFG_SLMS3["dirs"]["mutect2_depth_filt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.depthfilt.mutect2.combined.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["mutect2_depth_filt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.depthfilt.mutect2.combined.vcf.gz.tbi"
    log:
        stderr = CFG_SLMS3["logs"]["mutect2_depth_filt"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_depth_filt.stderr.log"
    conda: 
        CFG_SLMS3["conda_envs"]["bcftools"]
    resources: 
        **CFG_SLMS3["resources"]["mutect2_depth_filt"]
    threads: 
        CFG_SLMS3["threads"]["mutect2_depth_filt"]
    shell: 
        op.as_one_line("""
        tsamp=$(zgrep "##tumor_sample=" {input.vcf} | sed 's|##tumor_sample=||g');
        nsamp=$(zgrep "##normal_sample=" {input.vcf} | sed 's|##normal_sample=||g');
        bcftools view {input.vcf} | 
        sed "s|$tsamp|TUMOR|g" | sed "s|$nsamp|NORMAL|g" | 
        bcftools view -i 'FMT/DP[@{input.table}] >= 10 && FMT/AD[@{input.table}:1] >= 4 && FMT/AF[@{input.table}:0] >= 0.1' 
        -Oz -o {output.vcf} 2> {log.stderr} && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

# Update Starfish config to use outputs generated here

snakemake.utils.update_config(config["lcr-modules"]["starfish"], {
    "dirs": {
        "_parent": CFG_SLMS3["dirs"]["_parent"] + "starfish-" + CFG_SLMS3["module_versions"]["starfish"]
    },
    "inputs": {
        "vcf": {
            "mutect2": str(rules._slms_3_mutect2_depth_filt.output.vcf),
            "sage": str(rules._slms_3_input_sage_vcf.output.vcf), 
            "lofreq": str(rules._slms_3_annotate_lofreq_gnomad.output.vcf), 
            "strelka": str(rules._slms_3_annotate_strelka_gnomad.output.vcf)
        }
    }
})

include: lcr_modules + "/starfish/" + CFG_SLMS3["module_versions"]["starfish"] + "/starfish.smk"

# Symlinks the final output files into the module results directory (under '99-outputs/')

def _slms_3_get_starfish_output(wildcards): 
    CFG = config["lcr-modules"]["starfish"]
    # Get the path to the output directory from the checkpoint _starfish_rename_output
    # **wildcards unpacks the wildcards object from the rule where this input function is called
    checkpoint_output = os.path.dirname(checkpoints._starfish_rename_output.get(**wildcards).output.renamed)
    # Obtain the vcf_name wildcards using the glob_wildcards function
    vcf = expand(
        CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz",
        **wildcards,
        vcf_name = glob_wildcards(os.path.join(checkpoint_output, "{vcf_name, '3+'}.vcf.gz")).vcf_name
        )
    return vcf

rule _slms_3_output_vcf:
    input:
        _slms_3_get_starfish_output, 
    output:
        vcf = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.final.vcf.gz", 
        tbi = CFG_SLMS3["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.final.vcf.gz.tbi"
    run:
        op.relative_symlink(input, output.vcf, in_module = True)
        op.relative_symlink(input + ".tbi", output.tbi, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _slms_3_all:
    input:
        # rules._strelka_all.input, 
        # rules._lofreq_all.input, 
        # rules._mutect2_all.input, 
        # rules._sage_all.input, 
        # rules._starfish_all.input, 
        expand(
            [
                str(rules._slms_3_output_vcf.output.vcf),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG_SLMS3["runs"]["tumour_seq_type"],
            genome_build=CFG_SLMS3["runs"]["tumour_genome_build"],
            tumour_id=CFG_SLMS3["runs"]["tumour_sample_id"],
            normal_id=CFG_SLMS3["runs"]["normal_sample_id"],
            pair_status=CFG_SLMS3["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
# op.cleanup_module(CFG_SLMS3)
