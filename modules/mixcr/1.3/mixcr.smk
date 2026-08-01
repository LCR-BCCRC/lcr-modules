#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Anita dos Santos
# Module Author:    Manuela Cruz
# Contributors:     Laura Hilton
# 1.4:              adapted from 1.2 for MiXCR 4.x (analyze presets, JDK17, license, native SHM).


##### SETUP #####


import oncopipe as op
import os

min_oncopipe_version="1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

CFG = op.setup_module(
    name = "mixcr",
    version = "1.3",
    subdirectories = ["inputs", "mixcr", "outputs"],
)

# MiXCR 4.x requires a MiLaboratories license (activated per job via MI_LICENSE_FILE). Require the
# path up front so a missing license fails here with a clear message, not deep inside a job.
_license = CFG["inputs"].get("mixcr_license", "")
assert _license and _license != "__UPDATE__", (
    "MiXCR 4.x needs a license: set inputs.mixcr_license to the path of your MiLaboratories "
    "license file in the project config (gitignored / never committed). "
    "Obtain one at https://licensing.milaboratories.com/."
)
assert os.path.exists(_license), (
    f"inputs.mixcr_license = '{_license}' but that file does not exist; place your MiLaboratories "
    "license there before running."
)

RECEPTORS = CFG["receptors"]

ig_type = ["IGH","IGK","IGL"]
tr_type = ["TRA", "TRB", "TRD", "TRG"]

RANGE = ig_type + tr_type

if isinstance(RECEPTORS, str):
    RECEPTORS = RECEPTORS.split(" ")
if isinstance(RECEPTORS, list) and len(RECEPTORS) == 1:
    if RECEPTORS[0] == "ALL":
        RECEPTORS = ig_type + tr_type
    if RECEPTORS[0] == "BCR":
        RECEPTORS = ig_type
    if RECEPTORS[0] == "TCR":
        RECEPTORS = tr_type

assert all(receptor in RANGE for receptor in RECEPTORS), (
    "Config 'receptors' value is necessary for specifying target receptors. "
    "Ensure desired receptors are included in config and uppercase. "
    "Choose from: 'ALL', 'BCR', 'TCR' or list of IGH, IGK, IGL, TRA, TRB, TRD, TRG. "
)

# Define rules to be run locally when using a compute cluster
localrules:
    _install_mixcr,
    _mixcr_input_fastq,
    _mixcr_output_txt,
    _mixcr_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mixcr_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"],
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R2.fastq.gz",
    run:
        op.absolute_symlink(input.fastq_1, output.fastq_1)
        op.absolute_symlink(input.fastq_2, output.fastq_2)

# Installs MiXCR from github if not already present
rule _install_mixcr:
    params:
        mixcr = CFG["inputs"]["mixcr_exec"],
        version = CFG["options"]["mixcr_version"]
    output:
        complete = CFG["inputs"]["mixcr_exec"] + "/mixcr_dependencies_installed.success"
    shell:
        '''
        download_url="https://github.com/milaboratory/mixcr/releases/download/v{params.version}/mixcr-{params.version}.zip";
        mkdir -p {params.mixcr};

        if [ ! -f {params.mixcr}/mixcr ]; then
            wget -cO - $download_url > {params.mixcr}/mixcr.zip && unzip {params.mixcr}/mixcr.zip -d {params.mixcr}/ && rm {params.mixcr}/mixcr.zip;
            mv {params.mixcr}/mixcr*/* {params.mixcr}/ && rm -r {params.mixcr}/mixcr*/;
        fi

        touch  {output.complete};
        '''

# MiXCR 4.x is split into analyze (heavy, ~hours) and export (cheap, seconds) so an
# export tweak or transient failure never re-triggers the expensive analyze. The
# rna-seq preset ends in assembleContigs, so the checkpoint clns is <prefix>.contigs.clns.
rule _mixcr_analyze:
    input:
        fastq_1 = str(rules._mixcr_input_fastq.output.fastq_1),
        fastq_2 = str(rules._mixcr_input_fastq.output.fastq_2),
        fastq_1_real = CFG["inputs"]["sample_fastq_1"], # prevent premature deletion of temp fastqs
        fastq_2_real = CFG["inputs"]["sample_fastq_2"],
        installed = str(rules._install_mixcr.output.complete)
    output:
        clns = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.contigs.clns"
    log:
        stdout = CFG["logs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr_analyze.stdout.log",
        stderr = CFG["logs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr_analyze.stderr.log"
    resources:
        **CFG["resources"]["mixcr_run"]
    params:
        preset = op.switch_on_wildcard("seq_type", CFG["options"]["analyze"]),
        prefix = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}",
        mixcr = CFG["inputs"]["mixcr_exec"] + "/mixcr",
        license = CFG["inputs"]["mixcr_license"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda: CFG["conda_envs"]["java"]
    container: CFG["container_envs"]["java"]
    threads:
        CFG["threads"]["mixcr_run"]
    shell:
        op.as_one_line("""
        export JAVA_OPTS="-Xmx{params.jvmheap}m";
        export MI_LICENSE_FILE="{params.license}";
        {params.mixcr} analyze {params.preset} -t {threads} -f
        {input.fastq_1} {input.fastq_2} {params.prefix} > {log.stdout} 2> {log.stderr};
        """)

# Cheap: clns -> ALL + per-chain exportClones (--dont-split-files, .tsv->.txt) + report.
rule _mixcr_export:
    input:
        clns = str(rules._mixcr_analyze.output.clns)
    output:
        txt = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clonotypes.ALL.txt",
        report = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.report",
        results = expand(CFG["dirs"]["mixcr"] + "{{seq_type}}/{{sample_id}}/mixcr.{{sample_id}}.clonotypes.{chain}.txt", chain = RECEPTORS)
    log:
        stdout = CFG["logs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr_export.stdout.log",
        stderr = CFG["logs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr_export.stderr.log"
    resources:
        mem_mb = 8000
    params:
        export = CFG["options"]["export_clones"],
        prefix = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}",
        mixcr = CFG["inputs"]["mixcr_exec"] + "/mixcr",
        license = CFG["inputs"]["mixcr_license"],
        chains = " ".join(RECEPTORS)
    conda: CFG["conda_envs"]["java"]
    container: CFG["container_envs"]["java"]
    shell:
        op.as_one_line("""
        export JAVA_OPTS="-Xmx6500m";
        export MI_LICENSE_FILE="{params.license}";
        {params.mixcr} exportClones {params.export} --dont-split-files -f {input.clns} {output.txt}.tsv > {log.stdout} 2> {log.stderr};
        [ -s {output.txt}.tsv ] || : > {output.txt}.tsv;
        mv {output.txt}.tsv {output.txt};
        for chain in {params.chains}; do
        {params.mixcr} exportClones --chains $chain {params.export} --dont-split-files -f {input.clns}
        {params.prefix}.clonotypes.$chain.txt.tsv >> {log.stdout} 2>> {log.stderr};
        [ -s {params.prefix}.clonotypes.$chain.txt.tsv ] || head -1 {output.txt} > {params.prefix}.clonotypes.$chain.txt.tsv;
        mv {params.prefix}.clonotypes.$chain.txt.tsv {params.prefix}.clonotypes.$chain.txt;
        done;
        {params.mixcr} exportReports {input.clns} {output.report}.txt >> {log.stdout} 2>> {log.stderr};
        [ -s {output.report}.txt ] || : > {output.report}.txt;
        mv {output.report}.txt {output.report};
        """)

rule _mixcr_output_txt:
    input:
        results = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clonotypes.{chain}.txt"
    output:
        results = CFG["dirs"]["outputs"] + "txt/{seq_type}/mixcr.{sample_id}.clonotypes.{chain}.txt"
    wildcard_constraints:
        chain = '[A-Z]+'
    run:
        op.relative_symlink(input.results, output.results, in_module=True)

# Generates the target sentinels for each run, which generate the symlinks

rule _mixcr_all:
    input:
        expand(
            expand(
                [
                    str(rules._install_mixcr.output.complete),
                    str(rules._mixcr_output_txt.output.results),
                ],
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["samples"]["seq_type"],
                sample_id=CFG["samples"]["sample_id"],
                allow_missing=True),
            chain=RECEPTORS)


##### CLEANUP #####


op.cleanup_module(CFG)
