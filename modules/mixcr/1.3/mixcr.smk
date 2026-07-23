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

assert type(CFG["igblastn"])==bool, (
    "Ensure 'igblastn' is set to either True or False in config. "
    "True: also runs IgBLAST reannotation (% identity to IMGT) alongside MiXCR's native V identity."
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
    _mixcr_to_fasta,
    _igblastn_run,
    _update_mixcr_results,
    _symlink_mixcr_update,
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

# Installs the latest MiXCR release (4.x) from github if not already present
rule _install_mixcr:
    params:
        mixcr = CFG["inputs"]["mixcr_exec"]
    output:
        complete = CFG["inputs"]["mixcr_exec"] + "/mixcr_dependencies_installed.success"
    shell:
        '''
        download_url=$(curl --silent "https://api.github.com/repos/milaboratory/mixcr/releases/latest" | grep '"browser_download_url":' | sed -E 's/.*\"([^\"]+)\".*/\\1/');
        mkdir -p {params.mixcr};

        if [ ! -f {params.mixcr}/mixcr ]; then
            wget -cO - $download_url > {params.mixcr}/mixcr.zip && unzip {params.mixcr}/mixcr.zip -d {params.mixcr}/ && rm {params.mixcr}/mixcr.zip;
            mv {params.mixcr}/mixcr*/* {params.mixcr}/ && rm -r {params.mixcr}/mixcr*/;
        fi

        touch  {output.complete};
        '''

# MiXCR 4.x: analyze (preset) -> ALL + per-chain exportClones -> report.
# License activated per job via MI_LICENSE_FILE.
rule _mixcr_run:
    input:
        fastq_1 = str(rules._mixcr_input_fastq.output.fastq_1),
        fastq_2 = str(rules._mixcr_input_fastq.output.fastq_2),
        fastq_1_real = CFG["inputs"]["sample_fastq_1"], # prevent premature deletion of temp fastqs
        fastq_2_real = CFG["inputs"]["sample_fastq_2"],
        installed = str(rules._install_mixcr.output.complete)
    output:
        txt = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clonotypes.ALL.txt",
        report = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.report",
        results = expand(CFG["dirs"]["mixcr"] + "{{seq_type}}/{{sample_id}}/mixcr.{{sample_id}}.clonotypes.{chain}.txt", chain = RECEPTORS)
    log:
        stdout = CFG["logs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr_run.stdout.log",
        stderr = CFG["logs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr_run.stderr.log"
    resources:
        **CFG["resources"]["mixcr_run"]
    params:
        preset = op.switch_on_wildcard("seq_type", CFG["options"]["analyze"]),
        export = CFG["options"]["export_clones"],
        prefix = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}",
        # rna-seq preset ends in assembleContigs -> final clns is <prefix>.contigs.clns
        clns = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.contigs.clns",
        mixcr = CFG["inputs"]["mixcr_exec"] + "/mixcr",
        license = CFG["inputs"]["mixcr_license"],
        chains = " ".join(RECEPTORS),
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
        {params.mixcr} exportClones {params.export} -f {params.clns} {output.txt} >> {log.stdout} 2>> {log.stderr};
        for chain in {params.chains}; do
        {params.mixcr} exportClones --chains $chain {params.export} -f {params.clns}
        {params.prefix}.clonotypes.$chain.txt >> {log.stdout} 2>> {log.stderr};
        done;
        {params.mixcr} exportReports {params.clns} {output.report} >> {log.stdout} 2>> {log.stderr};
        touch "{output.txt}";
        """)

if CFG["igblastn"]:

    # Optional: reannotate each clone with IgBLAST IMGT % identity, next to MiXCR's native identity.

    rule _mixcr_to_fasta:
        input:
            mixcr_finished = str(rules._mixcr_run.output.txt),
            mixcr_chains = rules._mixcr_run.output.results,
            mixcr_results = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clonotypes.{chain}.txt",
            script = CFG["igblast_scripts"]["mixcr2fasta"]
        output:
            fasta = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clonotypes.{chain}.VDJseq.fasta",
            seq_info = temp(CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clonotypes.{chain}.regions.txt")
        shell:
            "{input.script} -i {input.mixcr_results} -o {output.fasta} -s {output.seq_info}"

    receptor_dict = {"IGH": "ig", "IGK": "ig", "IGL": "ig", "TRA": "tr", "TRB": "tr", "TRD":"tr", "TRG": "tr"}
    run_dict = {"IGH":"Ig", "IGK":"Ig", "IGL":"Ig", "TRA":"TCR","TRB":"TCR","TRD":"TCR","TRG":"TCR"}

    rule _igblastn_run:
        input:
            fasta = str(rules._mixcr_to_fasta.output.fasta),
            ig_db = reference_files("genomes/no_build/igblast/database/imgt_database.success")
        output:
            db = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clonotypes.{chain}.igblastn.fmt7"
        params:
            receptor_type = lambda wildcards: run_dict[wildcards.chain],
            aux = reference_files("downloads/igblast/optional_file/human_gl.aux"),
            gdv = lambda wildcards: (reference_files("genomes/no_build/igblast/database/imgt_human_" + receptor_dict[wildcards.chain] + "_v.ndb")).replace(".ndb",""),
            gdj = lambda wildcards: (reference_files("genomes/no_build/igblast/database/imgt_human_" + receptor_dict[wildcards.chain] + "_j.ndb")).replace(".ndb",""),
            gdd = lambda wildcards: (reference_files("genomes/no_build/igblast/database/imgt_human_" + receptor_dict[wildcards.chain] + "_d.ndb")).replace(".ndb",""),
            run_flags = CFG["options"]["igblast_run"]["run_flags"],
            form = "7 std btop " + CFG["options"]["igblast_run"]["form"]
        conda:
            CFG["conda_envs"]["igblast"]
        container:
            CFG["container_envs"]["igblast"]
        shell:
            op.as_one_line("""
            igblastn -query {input.fasta} -out {output.db}
            -ig_seqtype {params.receptor_type} -organism human
            -auxiliary_data {params.aux}
            -germline_db_V {params.gdv} -germline_db_J {params.gdj} -germline_db_D {params.gdd}
            {params.run_flags} -outfmt '{params.form}' -domain_system imgt
            """)

    rule _update_mixcr_results:
        input:
            db = str(rules._igblastn_run.output.db),
            og_mixcr = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clonotypes.{chain}.txt",
            seq_info = str(rules._mixcr_to_fasta.output.seq_info),
            script = CFG["igblast_scripts"]["igblastn2mixcr"]
        output:
            txt = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clonotypes.{chain}.igblast.txt"
        shell:
            "{input.script} -d {input.db} -m {input.og_mixcr} -o {output.txt} -s {input.seq_info}"

    rule _symlink_mixcr_update:
        input:
            txt = str(rules._update_mixcr_results.output.txt)
        output:
            txt = CFG["dirs"]["outputs"] + "txt/{seq_type}/mixcr.{sample_id}.clonotypes.{chain}.igblast.txt"
        wildcard_constraints:
            chain = "|".join(RECEPTORS)
        run:
            op.relative_symlink(input.txt, output.txt, in_module=True)


# Symlinks the final output files into the module results directory (under '99-outputs/')
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

if CFG["igblastn"]:
    rule _mixcr_all:
        input:
            expand(
                expand(
                    [
                        str(rules._install_mixcr.output.complete),
                        str(rules._symlink_mixcr_update.output.txt),
                        str(rules._mixcr_output_txt.output.results),
                    ],
                    zip,
                    seq_type=CFG["samples"]["seq_type"],
                    sample_id=CFG["samples"]["sample_id"],
                    allow_missing=True),
                chain=RECEPTORS)
elif not CFG["igblastn"]:
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
