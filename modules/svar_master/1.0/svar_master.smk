#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
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
# `CFG` is a shortcut to `config["lcr-modules"]["svar_master"]`
CFG_SV = op.setup_module(
    name = "svar_master",
    version = "1.0",
    subdirectories = ["inputs", "intersect_svs", "annotate_svs", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _svar_master_input_gridss,
    _svar_master_input_manta,
    _svar_master_output_native, 
    _svar_master_output_lifted,
    _svar_master_all


##### MODULE CONFIGFILES #####

# Load default module configs and update with values from svar_master config

min_module_versions = {
    "gridss": "2.0", 
    "manta": "2.0", 
    "hmftools": "1.1", 
    "liftover": "2.0"
}

for tool in CFG_SV["module_versions"].keys():
    # Ensure specified module version meets minimum requirements
    assert float(CFG_SV["module_versions"][tool]) >= float(min_module_versions[tool]), (
        f"The current {tool} module version is {CFG_SV['module_versions'][tool]}. "
        f"SVAR_MASTER requires {tool} module version {min_module_versions[tool]} or higher. "
    )
    lcr_modules = config["lcr-modules"]["_shared"]["lcr-modules"] + "modules/"
    configfile: lcr_modules + "/" +  tool + "/" + CFG_SV["module_versions"][tool] + "/config/default.yaml"
    snakemake.utils.update_config(
        config["lcr-modules"][tool], 
        config["lcr-modules"]["svar_master"][tool]
    )
    if not tool == "liftover": 
        config["lcr-modules"][tool]["inputs"]["sample_bam"] = CFG_SV["inputs"]["sample_bam"]
        config["lcr-modules"][tool]["inputs"]["sample_bai"] = CFG_SV["inputs"]["sample_bai"]    


##### MODULE SNAKEFILES #####

# Run Manta and GRIDSS 
include: "../../manta/" + CFG_SV["module_versions"]["manta"] + "/manta.smk"
include: "../../gridss/" + CFG_SV["module_versions"]["gridss"] + "/gridss.smk"

# Update hmftools config to use GRIDSS outputs
config['lcr-modules']['hmftools']['inputs']['slms3_vcf'] =  CFG_SV["inputs"]["slms_3_dir"] + "/99-outputs/vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.final.vcf.gz"
config['lcr-modules']['hmftools']['inputs']["gridss_somatic"] = str(rules._gridss_output_somatic_vcf.output.somatic)
config['lcr-modules']['hmftools']['inputs']["gridss_somatic_tbi"] = str(rules._gridss_output_somatic_vcf.output.somatic_tbi)
config['lcr-modules']['hmftools']['inputs']["gridss_somatic_filtered"] = str(rules._gridss_output_somatic_vcf.output.filtered)
config['lcr-modules']['hmftools']['inputs']["gridss_somatic_filtered_tbi"] = str(rules._gridss_output_somatic_vcf.output.filtered_tbi)

# Run hmftools 
include: "../../hmftools/" + CFG_SV["module_versions"]["hmftools"] + "/hmftools.smk"

##### Intersect Manta and GRIDSS #####

rule _svar_master_input_gridss:
    input:
        vcf = str(rules._gridss_output_somatic_vcf.output.filtered), 
        tbi = str(rules._gridss_output_somatic_vcf.output.filtered_tbi)
    output:
        vcf = CFG_SV["dirs"]["inputs"] + "sv_vcfs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss.somatic_filtered.vcf.gz", 
        tbi = CFG_SV["dirs"]["inputs"] + "sv_vcfs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss.somatic_filtered.vcf.gz.tbi"
    group: 
        "input_and_step_1"
    run:
        op.absolute_symlink(input.vcf, output.vcf)
        op.absolute_symlink(input.tbi, output.tbi)

def get_manta_vcf(wildcards): 
    if wildcards.pair_status in ["matched", "unmatched"]: 
        vcf_name = "somaticSV"
    elif wildcards.pair_status == "no_normal": 
        vcf_name = "tumorSV"
    vcf = expand(
        str(rules._manta_output_vcf.output.vcf), vcf_name = vcf_name, allow_missing=True
    )
    return vcf

rule _svar_master_input_manta:
    input:
        get_manta_vcf
    output:
        vcf = CFG_SV["dirs"]["inputs"] + "sv_vcfs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.manta.vcf", 
    group: 
        "input_and_step_1"
    run:
        op.absolute_symlink(input, output.vcf)

rule _svar_master_intersect: 
    input: 
        gridss = str(rules._svar_master_input_gridss.output.vcf), 
        manta = str(rules._svar_master_input_manta.output.vcf)
    output: 
        bedpe = CFG_SV["dirs"]["intersect_svs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_manta.bedpe" 
    log: 
        CFG_SV["logs"]["intersect_svs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.log", 
        CFG_SV["logs"]["intersect_svs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.RData"
    params: 
        regions_bed = op.switch_on_wildcard("genome_build", CFG_SV["switches"]["intersect"]["regions_bed"]), 
        bedpe_cols = CFG_SV["options"]["intersect"]["bedpe_cols"], 
        minvaf = CFG_SV["options"]["intersect"]["minvaf"], 
        mindp = CFG_SV["options"]["intersect"]["mindp"], 
        maxgap = CFG_SV["options"]["intersect"]["maxgap"]
    conda: 
        CFG_SV["conda_envs"]["filter_svs"]
    threads: CFG_SV["threads"]["intersect"]
    resources: 
        **CFG_SV["resources"]["intersect"]
    script: 
        CFG_SV["options"]["intersect"]["combine_svs"]

rule _svar_master_annotate: 
    input: 
        bedpe = str(rules._svar_master_intersect.output.bedpe)
    output: 
        a_bed = temp(CFG_SV["dirs"]["annotate_svs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/A_combined.bed"), 
        b_bed = temp(CFG_SV["dirs"]["annotate_svs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/B_combined.bed")
    log: 
        CFG_SV["logs"]["annotate_svs"] + "log/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/bedtools_closest.log"
    params: 
        annotations = op.switch_on_wildcard("genome_build", CFG_SV["switches"]["annotate"]["annotation_bed"])
    conda: 
        CFG_SV["conda_envs"]["bedtools"]
    threads: 
        CFG_SV["threads"]["annotate"]
    resources: 
        **CFG_SV["resources"]["annotate"]
    shell: 
        op.as_one_line("""
        tail -n+2 {input.bedpe} | 
        cut -f 1-3,7-15 | sort -k1,1 -k2,2n | 
        bedtools closest -D ref -a stdin -b {params.annotations} > {output.a_bed}; 
        tail -n+2 {input.bedpe} | 
        cut -f 4-15 | sort -k1,1 -k2,2n | 
        bedtools closest -D ref -a stdin -b {params.annotations} > {output.b_bed};
        """) 

rule _svar_master_annotate_combine: 
    input: 
        a_combined = str(rules._svar_master_annotate.output.a_bed), 
        b_combined = str(rules._svar_master_annotate.output.b_bed) 
    output: 
        bedpe = CFG_SV["dirs"]["annotate_svs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotated.bedpe"
    log: 
        CFG_SV["dirs"]["annotate_svs"] + "log/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/bedtools_closest.log"
    conda: 
        CFG_SV["conda_envs"]["filter_svs"]
    threads: 
        CFG_SV["threads"]["combine"]
    resources: 
        **CFG_SV["resources"]["combine"]
    script: 
        CFG_SV["options"]["combine_annotated"]["script"]


config["lcr-modules"]["liftover"]["dirs"]["_parent"] = CFG_SV["dirs"]["_parent"] + "liftover-" + CFG_SV["module_versions"]["liftover"]
config["lcr-modules"]["liftover"]["inputs"]["sample_file"] = str(rules._svar_master_annotate_combine.output.bedpe)

include: "../../liftover/" + CFG_SV["module_versions"]["liftover"] + "/liftover.smk"

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _svar_master_output_native:
    input:
        native = str(rules._liftover_input_file.output.another_tsv) 
    output:
        native = CFG_SV["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}.native.bedpe"
    run:
        op.relative_symlink(input.native, output.native, in_module = True)

rule _svar_master_output_lifted:
    input:
        lifted = str(rules._liftover_output.output)
    output:
        lifted = CFG_SV["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}.lifted_{chain}.bedpe"
    run:
        op.relative_symlink(input.lifted, output.lifted, in_module = True)



# Generates the target sentinels for each run, which generate the symlinks
rule _svar_master_all:
    input:
        rules._gridss_all.input, 
        rules._manta_all.input, 
        rules._hmftools_all.input, 
        rules._liftover_all.input, 
        expand(
            [
                str(rules._svar_master_output_native.output.native),
                str(rules._svar_master_output_lifted.output.lifted)
            ],
            zip,  # Run expand() with zip(), not product()
            tumour_id=CFG_SV["runs"]["tumour_sample_id"],
            normal_id=CFG_SV["runs"]["normal_sample_id"],
            genome_build = CFG_SV["runs"]["tumour_genome_build"],
            seq_type=CFG_SV["runs"]["tumour_seq_type"],
            pair_status=CFG_SV["runs"]["pair_status"],
            #repeat the tool name N times in expand so each pair in run is used
            tool=[CFG_SV["liftover"]["tool"]] * len(CFG_SV["runs"]["tumour_sample_id"]),
            chain=["hg38ToHg19" if "38" in str(x) else "hg19ToHg38" for x in CFG_SV["runs"]["tumour_genome_build"]]
            )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
CFG=CFG_SV
op.cleanup_module(CFG)
