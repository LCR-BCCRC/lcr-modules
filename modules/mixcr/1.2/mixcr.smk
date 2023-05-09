#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Anita dos Santos
# Module Author:    Manuela Cruz
# Contributors:     Laura Hilton


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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mixcr"]`
CFG = op.setup_module(
    name = "mixcr",
    version = "1.2",
    subdirectories = ["inputs", "mixcr", "outputs"],
)

assert type(CFG["igblastn"])==bool, (
    "Ensure 'igblastn' is set to either True or False in config. "
    "True: Runs igblastn alignment to retrieve % identity to IMGT sequences. This requires the '--contig-assembly' parameter when running MiXCR and may increase running time. "
    )

if CFG["igblastn"]:
    parameters_ok = False
    for seq_type in list(CFG["samples"]["seq_type"]):
        seq_type_parameters = CFG["options"]["mixcr_run"][seq_type]
        if "--contig-assembly" in seq_type_parameters: 
            parameters_ok = True
        assert parameters_ok, (
            "Config 'igblastn' value is set to True, but ['mixcr_run']['options'] does not include '--contig-assembly'. "
            "Include '--contig-assembly' in CFG['mixcr_run']['options'] to run igblastn.\n "
            "\n***** NOTE: using '--contig-assembly' may increase run duration. *****\n"
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

config_run_parameters = CFG["options"]["mixcr_run"]

for seq_type, run_parameters in config_run_parameters.items():
    if seq_type in list(CFG["samples"]["seq_type"]):
        override = False
        if "--receptor-type" in run_parameters:
            if not "--receptor-type xcr" in run_parameters:
                param_list = run_parameters.split("--")
                receptor_param = list(filter(lambda x: x.startswith("receptor-type"), param_list))
                rec_type = receptor_param[0].split(' ')[1]
                if len(RECEPTORS) == 1:
                    desired_run_1 = "--receptor-type " + RECEPTORS[0].lower()
                    if RECEPTORS[0] in ig_type:
                        desired_run_2 = "--receptor-type bcr"
                    elif RECEPTORS[0] in tr_type:
                        desired_run_2 = "--receptor-type tcr"
                    if not desired_run_1 in run_parameters and not desired_run_2 in run_parameters:
                        override = True
                        desired_run = desired_run_2
                if len(RECEPTORS) > 1:
                    if all(receptor in ig_type for receptor in RECEPTORS):
                        desired_run = "--receptor-type bcr"
                        if rec_type != "bcr":
                            override = True
                    if all(receptor in tr_type for receptor in RECEPTORS):
                        desired_run = "--receptor-type tcr"
                        if rec_type != "tcr":
                            override = True
                    if any(receptor in ig_type for receptor in RECEPTORS) and any(receptor in tr_type for receptor in RECEPTORS):
                        desired_run = "--receptor-type xcr"
                        override = True

        if override:
            print(f"----- Desired receptors: {RECEPTORS} \n----- Replacing receptor type specified in config for {seq_type} run: '--receptor-type {rec_type}' to '{desired_run}'")
    
            override_parameters = run_parameters.replace("--receptor-type " + rec_type, desired_run)
            CFG["options"]["mixcr_run"][seq_type] = override_parameters

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

# Installs latest MiXCR release from github if the mixcr folder is not present yet
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

# Run MiXCR rule
rule _mixcr_run:
    input:
        fastq_1 = str(rules._mixcr_input_fastq.output.fastq_1),
        fastq_2 = str(rules._mixcr_input_fastq.output.fastq_2),
        fastq_1_real = CFG["inputs"]["sample_fastq_1"], # Prevent premature deletion of fastqs marked as temp
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
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["mixcr_run"]),
        prefix = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}", 
        mixcr = CFG["inputs"]["mixcr_exec"] + "/mixcr", 
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8) 
    conda: CFG["conda_envs"]["java"]
    threads:
        CFG["threads"]["mixcr_run"]
    shell:
        op.as_one_line("""
        {params.mixcr} analyze shotgun -Xmx{params.jvmheap}m
        -s hsa -t {threads} {params.opts} 
        {input.fastq_1} {input.fastq_2} 
        {params.prefix} > {log.stdout} 2> {log.stderr};
        touch "{output.txt}";
        """)
        
if CFG["igblastn"]:

    # Run IgBLAST if specified in config
        
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


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
