#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Giuliano Banco
# Module Author:    Giuliano Banco
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import os
import glob 
import string
import copy
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

# Setup module and store module-specific configuration in `CFG`
# `CFG_NANOVOTE` is a shortcut to `config["lcr-modules"]["nanovote"]`
CFG_NANOVOTE = op.setup_module(
    name = "nanovote",
    version = "1.0",
    subdirectories = ["inputs", "starfish_input", "starfish", "union", "starfish_output", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _nanovote_input_clairs_vcf,
    _nanovote_input_clairs_to_vcf,
    _nanovote_input_deepsomatic_vcf,
    _nanovote_starfish_input_vcf,
    _nanovote_starfish_rename_output,
    _nanovote_starfish_output_vcf,
    _nanovote_starfish_dispatch,
    _nanovote_rename_samples_all,
    _nanovote_output_vcf,
    _nanovote_all,


# Ensure each of the submodule versions meets requirements. 

assert float(CFG_NANOVOTE["module_versions"]["clairs"]) >= 1.1, (
    f"The current ClairS module version is {CFG_NANOVOTE['module_versions']['clairs']}. "
    "NanoVote requires ClairS module version 1.1 or higher. "
)

assert float(CFG_NANOVOTE["module_versions"]["clairs_to"]) >= 1.1, (
    f"The current ClairS-TO module version is {CFG_NANOVOTE['module_versions']['clairs_to']}. "
    "NanoVote requires ClairS-TO module version 1.1 or higher. "
)

# assert float(CFG_NANOVOTE["module_versions"]["starfish"]) >= 2.0, (
#     f"The current starfish module version is {CFG_NANOVOTE['module_versions']['starfish']}. "
#     "NanoVote requires starfish module version 2.0 or higher. "
# )

assert float(CFG_NANOVOTE["module_versions"]["deepsomatic"]) >= 1.0, (
    f"The current deepsomatic module version is {CFG_NANOVOTE['module_versions']['deepsomatic']}. "
    "NanoVote requires deepsomatic module version 1.0 or higher. "
)

# Load default module configs and update with values from nanovote config
# Create NanoVote-level normal_name from tumour_chemistry
normal_map = CFG_NANOVOTE["options"]["normal_name"]

CFG_NANOVOTE["runs"]["normal_name"] = (
    CFG_NANOVOTE["runs"]["tumour_chemistry"].map(normal_map)
)

if CFG_NANOVOTE["runs"]["normal_name"].isna().any():
    missing = CFG_NANOVOTE["runs"].loc[
        CFG_NANOVOTE["runs"]["normal_name"].isna(),
        "tumour_chemistry"
    ].unique()
    raise ValueError(
        f"Some chemistry values have no normal_name mapping: {missing}. "
        f"Available mappings are: {normal_map}"
    )


# Load default module configs and update with values from NanoVote config
for tool in ["clairs", "clairs_to", "deepsomatic"]:
    lcr_modules = config["lcr-modules"]["_shared"]["lcr-modules"] + "modules/"
    configfile: lcr_modules + "/" + tool + "/" + CFG_NANOVOTE["module_versions"][tool] + "/config/default.yaml"

    # Apply module-specific overrides from nanovote config, if present
    snakemake.utils.update_config(
        config["lcr-modules"][tool],
        config["lcr-modules"]["nanovote"].get(tool, {})
    )

    # Shared tumour BAM inputs
    config["lcr-modules"][tool]["inputs"]["sample_bam"] = CFG_NANOVOTE["inputs"]["sample_bam"]
    config["lcr-modules"][tool]["inputs"]["sample_bai"] = CFG_NANOVOTE["inputs"]["sample_bai"]

    # Only paired/unmatched-normal callers need normal BAMs and normal_name mapping
    if tool in ["clairs", "deepsomatic"]:
        config["lcr-modules"][tool]["options"]["normal_name"] = CFG_NANOVOTE["options"]["normal_name"]

        config["lcr-modules"][tool]["inputs"]["normal_bam"] = CFG_NANOVOTE["inputs"]["normal_bam"]
        config["lcr-modules"][tool]["inputs"]["normal_bai"] = CFG_NANOVOTE["inputs"]["normal_bai"]


# Determine how DeepSomatic is being run within NanoVote
NANOVOTE_DEEPSOMATIC_CALLING_MODE = (
    config["lcr-modules"]["deepsomatic"]["options"]["calling_mode"]
)

NANOVOTE_DEEPSOMATIC_CALLING_MODE_DIRS = {
    "unmatched": "unmatched",
    "tumor_only": "tumor-only",
}

if NANOVOTE_DEEPSOMATIC_CALLING_MODE not in NANOVOTE_DEEPSOMATIC_CALLING_MODE_DIRS:
    raise ValueError(
        "Invalid DeepSomatic calling_mode in NanoVote: "
        f"{NANOVOTE_DEEPSOMATIC_CALLING_MODE}. "
        "Expected 'unmatched' or 'tumor_only'."
    )

NANOVOTE_DEEPSOMATIC_CALLING_MODE_DIR = (
    NANOVOTE_DEEPSOMATIC_CALLING_MODE_DIRS[
        NANOVOTE_DEEPSOMATIC_CALLING_MODE
    ]
)

NANOVOTE_RUN_ID = (
    "{seq_type}--{genome_build}/"
    "{tumour_id}--{normal_name}--{chemistry}"
)

NANOVOTE_DEEPSOMATIC_RUN_ID = (
    "{seq_type}--{genome_build}/"
    "{tumour_id}--{normal_name}--{chemistry}--"
    + NANOVOTE_DEEPSOMATIC_CALLING_MODE_DIR
)


##### RULES #####

##### FIRST PASS VARIANT CALLING MODULE SNAKEFILES #####

# Load ClairS, ClairS-TO, and DeepSomatic
include: "../../clairs/" + CFG_NANOVOTE["module_versions"]["clairs"] + "/clairs.smk"
include: "../../clairs_to/" + CFG_NANOVOTE["module_versions"]["clairs_to"] + "/clairs_to.smk"
include: "../../deepsomatic/" + CFG_NANOVOTE["module_versions"]["deepsomatic"] + "/deepsomatic.smk"

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _nanovote_input_clairs_vcf: 
    input:
        vcf = str(rules._clairs_gnomad_annotation.output.vcf), 
        tbi = str(rules._clairs_gnomad_annotation.output.tbi),
        cleanup = str(rules._clairs_clean.output.cleanup)
    output:
        vcf = CFG_NANOVOTE["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.vcf.gz", 
        tbi = CFG_NANOVOTE["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _nanovote_input_clairs_to_vcf: 
    input:
        vcf = str(rules._clairs_to_gnomad_annotation.output.vcf), 
        tbi = str(rules._clairs_to_gnomad_annotation.output.tbi),
        cleanup = str(rules._clairs_to_clean.output.cleanup_complete)
    output:
        vcf = CFG_NANOVOTE["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--tumor_only/clairs_to.vcf.gz", 
        tbi = CFG_NANOVOTE["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--tumor_only/clairs_to.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _nanovote_input_deepsomatic_vcf:
    input:
        vcf = str(rules._deepsomatic_gnomad_annotation.output.vcf),
        tbi = str(rules._deepsomatic_gnomad_annotation.output.tbi)
    output:
        vcf = (CFG_NANOVOTE["dirs"]["inputs"] + "vcf/" + NANOVOTE_DEEPSOMATIC_RUN_ID + "/deepsomatic.vcf.gz"),
        tbi = (CFG_NANOVOTE["dirs"]["inputs"] + "vcf/" + NANOVOTE_DEEPSOMATIC_RUN_ID + "/deepsomatic.vcf.gz.tbi")
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module = True)
        op.relative_symlink(input.tbi, output.tbi, in_module = True)


##### NANOVOTE STARFISH/VOTING RULES #####

nanovote_callers = ["clairs", "clairs_to", "deepsomatic"]

nanovote_vcfs = {
    "clairs": str(rules._nanovote_input_clairs_vcf.output.vcf),
    "clairs_to": str(rules._nanovote_input_clairs_to_vcf.output.vcf),
    "deepsomatic": str(rules._nanovote_input_deepsomatic_vcf.output.vcf),
}

assert len(nanovote_callers) < 6, (
    "Starfish can only handle a maximum of 6 input VCFs."
)


def _nanovote_fill_wildcards(pattern, wildcards):
    return pattern.format(
        seq_type = wildcards.seq_type,
        genome_build = wildcards.genome_build,
        tumour_id = wildcards.tumour_id,
        normal_name = wildcards.normal_name,
        chemistry = wildcards.chemistry,
    )


# Symlink NanoVote caller VCFs into a common Starfish-style input location
rule _nanovote_starfish_input_vcf:
    input:
        vcf = lambda w: _nanovote_fill_wildcards(nanovote_vcfs[w.caller], w),
        tbi = lambda w: (_nanovote_fill_wildcards(nanovote_vcfs[w.caller], w) + ".tbi")
    output:
        vcf = (CFG_NANOVOTE["dirs"]["starfish_input"] + "vcf/" + NANOVOTE_RUN_ID + ".{caller}.vcf.gz"),
        tbi = (CFG_NANOVOTE["dirs"]["starfish_input"] + "vcf/" + NANOVOTE_RUN_ID + ".{caller}.vcf.gz.tbi")
    run:
        op.absolute_symlink(input.vcf, output.vcf)
        op.absolute_symlink(input.tbi, output.tbi)


# Run Starfish on ClairS, ClairS-TO, and DeepSomatic
rule _nanovote_starfish_run:
    input:
        vcfs = expand(
            CFG_NANOVOTE["dirs"]["starfish_input"] + "vcf/" + NANOVOTE_RUN_ID + ".{caller}.vcf.gz",
            caller = nanovote_callers,
            allow_missing = True
        ),
        tbis = expand(
            CFG_NANOVOTE["dirs"]["starfish_input"] + "vcf/" + NANOVOTE_RUN_ID + ".{caller}.vcf.gz.tbi",
            caller = nanovote_callers,
            allow_missing = True
        ),
        reference = ancient(reference_files("genomes/{genome_build}/sdf")),
        starfish_script = CFG_NANOVOTE["inputs"]["starfish_script"]
    output:
        complete = touch(CFG_NANOVOTE["dirs"]["starfish"] + NANOVOTE_RUN_ID + "/starfish.complete"),
        two_plus_vcf = (CFG_NANOVOTE["dirs"]["starfish"] + NANOVOTE_RUN_ID + "/2+.vcf.gz"),
        two_plus_tbi = (CFG_NANOVOTE["dirs"]["starfish"] + NANOVOTE_RUN_ID + "/2+.vcf.gz.tbi")
    log:
        stdout = CFG_NANOVOTE["logs"]["starfish"] + NANOVOTE_RUN_ID + "/starfish_run.stdout.log",
        stderr = CFG_NANOVOTE["logs"]["starfish"] + NANOVOTE_RUN_ID + "/starfish_run.stderr.log"
    params:
        opts = CFG_NANOVOTE["options"]["starfish_run"],
        caller_names = " ".join(nanovote_callers)
    conda:
        CFG_NANOVOTE["conda_envs"]["starfish"]
    container:
        CFG_NANOVOTE["container_envs"]["starfish"]
    threads:
        CFG_NANOVOTE["threads"]["starfish_run"]
    resources:
        **CFG_NANOVOTE["resources"]["starfish_run"]
    shell:
        op.as_one_line("""
        if [[ -e $(dirname {output.complete})/temp ]]; then 
            rm -rf $(dirname {output.complete})/temp; 
        fi
        &&
        {input.starfish_script} --sdf {input.reference} -O $(dirname {output.complete})
        --names {params.caller_names} --threads {threads}
        {params.opts} --venn --verbose
        -V {input.vcfs} > {log.stdout} 2> {log.stderr}
        """)


checkpoint _nanovote_starfish_rename_output:
    input:
        complete = CFG_NANOVOTE["dirs"]["starfish"] + NANOVOTE_RUN_ID + "/starfish.complete"
    output:
        renamed = touch(CFG_NANOVOTE["dirs"]["starfish"] + NANOVOTE_RUN_ID + "/starfish.renamed")
    run:
        names = list(string.ascii_uppercase[0:(len(nanovote_callers))])
        name_dict = {names[i]: nanovote_callers[i] for i in range(len(nanovote_callers))}
        outdir = os.path.dirname(input.complete)
        for src in glob.glob(outdir + "/[A-Z]*.vcf.gz*"):
            bname = os.path.basename(src)
            lhs, rhs = bname.split(".", 1)
            lhs = "_".join(lhs)
            for key in list(name_dict.keys()):
                lhs = lhs.replace(key, name_dict[key])
            dest = os.path.join(outdir, ".".join([lhs, rhs]))
            os.rename(src, dest)


rule _nanovote_starfish_output_vcf:
    input:
        vcf = (
            CFG_NANOVOTE["dirs"]["starfish"]
            + NANOVOTE_RUN_ID
            + "/{vcf_name}.vcf.gz"
        ),
        tbi = (
            CFG_NANOVOTE["dirs"]["starfish"]
            + NANOVOTE_RUN_ID
            + "/{vcf_name}.vcf.gz.tbi"
        )
    output:
        vcf = (
            CFG_NANOVOTE["dirs"]["starfish_output"]
            + "vcf/"
            + NANOVOTE_RUN_ID
            + ".{vcf_name}.vcf.gz"
        ),
        tbi = (
            CFG_NANOVOTE["dirs"]["starfish_output"]
            + "vcf/"
            + NANOVOTE_RUN_ID
            + ".{vcf_name}.vcf.gz.tbi"
        )
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module=True), 
        op.relative_symlink(input.tbi, output.tbi, in_module=True)


def _nanovote_starfish_get_output_target(wildcards):
    checkpoint_output = os.path.dirname(
        checkpoints._nanovote_starfish_rename_output.get(
            **wildcards
        ).output.renamed
    )
    vcf_names = glob_wildcards(
        os.path.join(checkpoint_output, "{vcf_name}.vcf.gz")
    ).vcf_name
    targets = expand(
        [
            (
                CFG_NANOVOTE["dirs"]["starfish_output"]
                + "vcf/"
                + NANOVOTE_RUN_ID
                + ".{vcf_name}.vcf.gz"
            ),
            (
                CFG_NANOVOTE["dirs"]["starfish_output"]
                + "vcf/"
                + NANOVOTE_RUN_ID
                + ".{vcf_name}.vcf.gz.tbi"
            ),
        ],
        **wildcards,
        vcf_name=vcf_names
    )
    return targets


rule _nanovote_starfish_dispatch:
    input:
        _nanovote_starfish_get_output_target
    output:
        touch(CFG_NANOVOTE["dirs"]["outputs"] + "dispatched/" + NANOVOTE_RUN_ID + ".dispatched")


# Create a final union VCF file summarizing which variants were called by which variant callers
# First rename the sample columns in each VCF as TUMOR_{caller} or NORMAL_{caller}
rule _nanovote_rename_samples_all:
    input:
        vcf = lambda w: _nanovote_fill_wildcards(
            nanovote_vcfs[w.caller],
            w
        )
    output:
        vcf = temp(
            CFG_NANOVOTE["dirs"]["union"]
            + "{seq_type}--{genome_build}/"
            + "{tumour_id}--{normal_name}--{chemistry}/"
            + "{caller}.tmp.vcf.gz"
        ),
        tbi = temp(
            CFG_NANOVOTE["dirs"]["union"]
            + "{seq_type}--{genome_build}/"
            + "{tumour_id}--{normal_name}--{chemistry}/"
            + "{caller}.tmp.vcf.gz.tbi"
        ),
        samples = temp(
            CFG_NANOVOTE["dirs"]["union"]
            + "{seq_type}--{genome_build}/"
            + "{tumour_id}--{normal_name}--{chemistry}/"
            + "{caller}.samples.txt"
        )
    log:
        stderr = (
            CFG_NANOVOTE["logs"]["union"]
            + "{seq_type}--{genome_build}/"
            + "{tumour_id}--{normal_name}--{chemistry}/"
            + "rename_samples_{caller}.stderr.log"
        )
    conda:
        CFG_NANOVOTE["conda_envs"]["bcftools"]
    container:
        CFG_NANOVOTE["container_envs"]["bcftools"]
    resources:
        **CFG_NANOVOTE["resources"]["rename_all"]
    threads:
        CFG_NANOVOTE["threads"]["rename_all"]
    shell:
        op.as_one_line("""
        printf
            "TUMOR\\tTUMOR_{wildcards.caller}\\n"
            > {output.samples}
        &&
        bcftools reheader
            -s {output.samples}
            -o {output.vcf}.reheadered.vcf.gz
            {input.vcf}
            2> {log.stderr}
        &&
        if bcftools view -h {output.vcf}.reheadered.vcf.gz |
            awk '/^##FORMAT=<ID=NAF,/ {{found=1}} END {{exit !found}}';
        then
            bcftools annotate
                -x FORMAT/NAF
                -Oz
                -o {output.vcf}
                {output.vcf}.reheadered.vcf.gz
                2>> {log.stderr};
        else
            mv
                {output.vcf}.reheadered.vcf.gz
                {output.vcf};
        fi
        &&
        tabix
            -f
            -p vcf
            {output.vcf}
            2>> {log.stderr}
        &&
        rm -f {output.vcf}.reheadered.vcf.gz
        """)


# Merge VCFs
rule _nanovote_union_vcf: 
    input: 
        vcf = expand(
            rules._nanovote_rename_samples_all.output.vcf, 
            caller = nanovote_callers,
            allow_missing = True
        ), 
        tbi = expand(
            rules._nanovote_rename_samples_all.output.tbi, 
            caller = nanovote_callers, 
            allow_missing = True
        )
    output: 
        vcf = CFG_NANOVOTE["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/union.vcf.gz", 
        tbi = CFG_NANOVOTE["dirs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/union.vcf.gz.tbi"
    log:
        stderr = CFG_NANOVOTE["logs"]["union"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/union.stderr.log"
    conda:
        CFG_NANOVOTE["conda_envs"]["bcftools"]
    container:
        CFG_NANOVOTE["container_envs"]["bcftools"]
    resources: 
        **CFG_NANOVOTE["resources"]["union_vcf"]
    threads: 
        CFG_NANOVOTE["threads"]["union_vcf"]
    shell: 
        op.as_one_line("""
        bcftools merge --threads {threads} -m both -Oz -o {output.vcf} {input.vcf} 2> {log.stderr}
        && 
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)

# Retain only 2+ calls that are also present in ClairS-TO
rule _nanovote_require_clairs_to:
    input:
        two_plus_vcf = str(rules._nanovote_starfish_run.output.two_plus_vcf),
        two_plus_tbi = str(rules._nanovote_starfish_run.output.two_plus_tbi),
        clairs_to_vcf = str(rules._nanovote_input_clairs_to_vcf.output.vcf),
        clairs_to_tbi = str(rules._nanovote_input_clairs_to_vcf.output.tbi)
    output:
        vcf = (
            CFG_NANOVOTE["dirs"]["starfish_output"]
            + NANOVOTE_RUN_ID
            + ".clairs_to_required.vcf.gz"
        ),
        tbi = (
            CFG_NANOVOTE["dirs"]["starfish_output"]
            + NANOVOTE_RUN_ID
            + ".clairs_to_required.vcf.gz.tbi"
        )
    log:
        stderr = (
            CFG_NANOVOTE["logs"]["starfish_output"]
            + NANOVOTE_RUN_ID
            + ".clairs_to_required.stderr.log"
        )
    conda:
        CFG_NANOVOTE["conda_envs"]["bcftools"]
    container:
        CFG_NANOVOTE["container_envs"]["bcftools"]
    resources:
        **CFG_NANOVOTE["resources"]["union_vcf"]
    threads:
        CFG_NANOVOTE["threads"]["union_vcf"]
    shell:
        op.as_one_line("""
        bcftools isec
            -n=2
            -w1
            {input.two_plus_vcf}
            {input.clairs_to_vcf}
            -Oz
            -o {output.vcf}
            2> {log.stderr}
        &&
        tabix
            -f
            -p vcf
            {output.vcf}
            2>> {log.stderr}
        """)


rule _nanovote_output_vcf:
    input:
        isec_vcf = str(rules._nanovote_require_clairs_to.output.vcf),
        isec_tbi = str(rules._nanovote_require_clairs_to.output.tbi),
        union_vcf = str(rules._nanovote_union_vcf.output.vcf),
        union_tbi = str(rules._nanovote_union_vcf.output.tbi),
        dispatched = str(rules._nanovote_starfish_dispatch.output),
        two_plus_vcf = str(rules._nanovote_starfish_run.output.two_plus_vcf),
        two_plus_tbi = str(rules._nanovote_starfish_run.output.two_plus_tbi)
    output:
        isec_vcf = (CFG_NANOVOTE["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}.nanovote.final.vcf.gz"),
        isec_tbi = (CFG_NANOVOTE["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}.nanovote.final.vcf.gz.tbi"),
        union_vcf = (CFG_NANOVOTE["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}.nanovote.union.vcf.gz"),
        union_tbi = (CFG_NANOVOTE["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}.nanovote.union.vcf.gz.tbi"),
        two_plus_vcf = (CFG_NANOVOTE["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}.nanovote.2+.vcf.gz"),
        two_plus_tbi = (CFG_NANOVOTE["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}.nanovote.2+.vcf.gz.tbi")
    params:
        rm_files = CFG_NANOVOTE["options"]["cleanup_vcfs"],
        files_to_rm = [
            str(rules._nanovote_input_clairs_vcf.output.vcf),
            str(rules._nanovote_input_clairs_vcf.output.tbi),
            str(rules._nanovote_input_clairs_to_vcf.output.vcf),
            str(rules._nanovote_input_clairs_to_vcf.output.tbi),
            str(rules._nanovote_input_deepsomatic_vcf.output.vcf),
            str(rules._nanovote_input_deepsomatic_vcf.output.tbi),
        ]
    run:
        op.relative_symlink(str(input.isec_vcf), str(output.isec_vcf), in_module=True)
        op.relative_symlink(str(input.isec_tbi), str(output.isec_tbi), in_module=True)
        op.relative_symlink(str(input.union_vcf), str(output.union_vcf), in_module=True)
        op.relative_symlink(str(input.union_tbi), str(output.union_tbi), in_module=True)
        op.relative_symlink(str(input.two_plus_vcf), str(output.two_plus_vcf), in_module=True)
        op.relative_symlink(str(input.two_plus_tbi), str(output.two_plus_tbi), in_module=True)
        if params.rm_files: 
            for file in params.files_to_rm: 
                if os.path.exists(file): 
                    os.remove(file)


# Generates the target sentinels for each run, which generate the symlinks
rule _nanovote_all:
    input:
        rules._clairs_to_all.input,
        rules._clairs_all.input,
        rules._deepsomatic_all.input,
        expand(
            [
                str(rules._nanovote_output_vcf.output.isec_vcf),
                str(rules._nanovote_output_vcf.output.isec_tbi),
                str(rules._nanovote_output_vcf.output.union_vcf),
                str(rules._nanovote_output_vcf.output.union_tbi),
                str(rules._nanovote_output_vcf.output.two_plus_vcf),
                str(rules._nanovote_output_vcf.output.two_plus_tbi)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type = CFG_NANOVOTE["runs"]["tumour_seq_type"],
            genome_build = CFG_NANOVOTE["runs"]["tumour_genome_build"],
            tumour_id = CFG_NANOVOTE["runs"]["tumour_sample_id"],
            platform = CFG_NANOVOTE["runs"]["tumour_platform"],
            chemistry = CFG_NANOVOTE["runs"]["tumour_chemistry"],
            normal_name = CFG_NANOVOTE["runs"]["normal_name"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
CFG=CFG_NANOVOTE
op.cleanup_module(CFG)
