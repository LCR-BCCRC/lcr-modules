#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Jasper Wong
# Module Author:    Jasper Wong
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
# `CFG` is a shortcut to `config["lcr-modules"]["spechla"]`
CFG = op.setup_module(
    name = "spechla",
    version = "1.0",
    subdirectories = ["inputs", "hla_reads", "spechla", "loh", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _spechla_input_bam,
    _spechla_setup,
    _spechla_index,
    _spechla_output_txt,
    _spechla_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _spechla_input_bam:
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


rule _spechla_setup:
    output:
        txt = touch(CFG["dirs"]["inputs"] + "SpecHLA_" + str(CFG["options"]["spechla_version"]) + "/clone.done")
    log:
        stdout = CFG["logs"]["inputs"] + "git_clone.stdout.log"
    params:
        scriptdir = CFG["dirs"]["inputs"] + "SpecHLA_" + str(CFG["options"]["spechla_version"]),
        url = "https://github.com/deepomicslab/SpecHLA/archive/refs/tags/v" + str(CFG["options"]["spechla_version"]) + ".tar.gz"
    shell:
        op.as_one_line("""
        mkdir -p {params.scriptdir} &&
        wget -qO- {params.url} | tar xzf - -C {params.scriptdir} --strip-components 1 > {log.stdout} 2>&1 && 
        chmod +x -R {params.scriptdir}/bin/*
        """)

rule _spechla_index:
    input:
        txt = str(rules._spechla_setup.output.txt)
    output:
        txt = touch(CFG["dirs"]["inputs"] + "SpecHLA_" + str(CFG["options"]["spechla_version"]) + "/index.done")
    log:
        stdout = CFG["logs"]["inputs"] + "index.stdout.log"
    params:
        scriptdir = CFG["dirs"]["inputs"] + "SpecHLA_" + str(CFG["options"]["spechla_version"])
    threads:
        CFG["threads"]["setup"]
    resources:
        **CFG["resources"]["setup"]  
    shell:
        op.as_one_line("""
        unset LD_LIBRARY_PATH && unset LIBRARY_PATH &&
        bash {params.scriptdir}/index.sh > {log.stdout} 2>&1
        """)

# Extract HLA reads from bam and cram
def _which_genome(wildcards):
    if "38" in str({wildcards.genome_build}):
        return "hg38"
    else:
        return "hg19"

rule _spechla_extract_reads:
    input:
        txt = str(rules._spechla_setup.output.txt),
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    output:
        fastq_1 = temp(CFG["dirs"]["hla_reads"] + "{seq_type}--{genome_build}/{sample_id}_extract_1.fq.gz"),
        fastq_2 = temp(CFG["dirs"]["hla_reads"] + "{seq_type}--{genome_build}/{sample_id}_extract_2.fq.gz"),
        fastq_unpaired = temp(CFG["dirs"]["hla_reads"] + "{seq_type}--{genome_build}/{sample_id}_extract.unpaired.fq.gz"),
    log:
        stderr = CFG["logs"]["hla_reads"] + "{seq_type}--{genome_build}/{sample_id}/extract_reads.stderr.log"
    params:
        scriptdir = CFG["dirs"]["inputs"] + "SpecHLA_" + str(CFG["options"]["spechla_version"]),
        genome = _which_genome,
        outdir = CFG["dirs"]["hla_reads"] + "{seq_type}--{genome_build}/"
    conda: CFG["conda_envs"]["spechla"]
    group: "extract_and_run"
    threads:
        CFG["threads"]["extract_reads"]
    resources:
        **CFG["resources"]["extract_reads"]  
    shell:
        op.as_one_line("""
        bash {params.scriptdir}/script/ExtractHLAread.sh -s {wildcards.sample_id} -b {input.bam} -r {params.genome} -o {params.outdir} > {log.stderr} 2>&1
        """)

#Choose full-length or exon typing [0|1]. 0 indicates full-length, 1 means exon,
#default is 0. With Exome or RNA data, must select 1 (i.e., exon typing).
def _which_seqtype(wildcards):
    if any(seq_type in str({wildcards.seq_type}) for seq_type in ['exome', 'capture', 'target', 'mrna']):
        return 1
    else:
        return 0

rule _spechla_hla_typing:
    input:
        fq_1 = CFG["dirs"]["hla_reads"] + "{seq_type}--{genome_build}/{sample_id}_extract_1.fq.gz",
        fq_2 = CFG["dirs"]["hla_reads"] + "{seq_type}--{genome_build}/{sample_id}_extract_2.fq.gz",
        fastq_unpaired = CFG["dirs"]["hla_reads"] + "{seq_type}--{genome_build}/{sample_id}_extract.unpaired.fq.gz"
    output:
        hla_results = CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}/hla.result.txt",
        realigned_bam = temp(CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.realign.sort.bam"),
        realigned_bam_depth = temp(CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.realign.sort.bam.depth"),
    log:
        stderr = CFG["logs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}/extract_reads.stderr.log"
    params:
        scriptdir = CFG["dirs"]["inputs"] + "SpecHLA_" + str(CFG["options"]["spechla_version"]),
        outdir = CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/",
        which_seqtype = _which_seqtype,
        var_quality = CFG["options"]["spechla"]["var_qual"],
        var_depth = CFG["options"]["spechla"]["var_depth"],
        min_depth_masking = CFG["options"]["spechla"]["min_depth_masking"],
        long_indels = CFG["options"]["spechla"]["long_indels"],
        opts = CFG["options"]["spechla"]["opts"]
    conda: CFG["conda_envs"]["spechla"]
    group: "extract_and_run"
    threads:
        CFG["threads"]["spechla"]
    resources:
        **CFG["resources"]["spechla"]  
    shell:
        op.as_one_line("""
        bash {params.scriptdir}/script/whole/SpecHLA.sh 
        -n {wildcards.sample_id} 
        -1 {input.fq_1} 
        -2 {input.fq_2} 
        -o {params.outdir} 
        -u {params.which_seqtype} 
        -j {threads} 
        -q {params.var_quality} 
        -s {params.var_depth} 
        -k {params.min_depth_masking} 
        -v {params.long_indels} {params.opts} > {log.stderr} 2>&1
        """)
        
# Remove unwanted outputs
rule _spechla_cleanup:
    input:
        realigned_bam = CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.realign.sort.bam",
        realigned_bam_depth = CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.realign.sort.bam.depth"
    output:
        done = touch(CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}/.removed.bloat")

# Examine LOH
rule _spechla_generate_freqlist:
    input:
        hla_results = str(rules._spechla_hla_typing.output.hla_results)
    output:
        freq_list = CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}/freq.list"
    params:
        outdir = CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}",
    shell:
        op.as_one_line("""
            touch {output.freq_list} && 
            ls {params.outdir}/*_freq.txt > {output.freq_list} || true
        """)

# Fetch purity and ploidy from a table
# Table should be formatted as "sample\tploidy\tpurity"
def _get_sample_ploidy(wildcards):
    CFG = config["lcr-modules"]["spechla"]
    purity_ploidy_file = CFG["options"]["purity_ploidy_file"] 
    with open(purity_ploidy_file, 'r') as f:
        for line in f.readlines():
            line = line.rstrip('\n').rstrip('\r')
            if str({wildcards.sample_id}) in line:
                cols = line.split('\t')
                ploidy = cols[2]
                valid_ploidy = ploidy.isnumeric()
                if valid_ploidy:
                    return ploidy
                else:
                    return 2
            else:
                return 2

def _get_sample_purity(wildcards):
    CFG = config["lcr-modules"]["spechla"]
    purity_ploidy_file = CFG["options"]["purity_ploidy_file"] 
    with open(purity_ploidy_file, 'r') as f:
        for line in f.readlines():
            line = line.rstrip('\n').rstrip('\r')
            if str({wildcards.sample_id}) in line:
                cols = line.split('\t')
                purity = cols[3]
                valid_purity = purity.isnumeric()
                if valid_purity:
                    return purity
                else:
                    return 1
            else:
                return 1

rule _spechla_loh:
    input:
        hla_results = str(rules._spechla_hla_typing.output.hla_results),
        cleanup = str(rules._spechla_cleanup.output.done),
        freq_list = str(rules._spechla_generate_freqlist.output.freq_list)
    output:
        loh_info = CFG["dirs"]["loh"] + "{seq_type}--{genome_build}/{sample_id}/merge.hla.copy.txt"
    log:
        stderr = CFG["logs"]["loh"] + "{seq_type}--{genome_build}/{sample_id}/loh.stderr.log"
    params:
        scriptdir = CFG["dirs"]["inputs"] + "SpecHLA_" + str(CFG["options"]["spechla_version"]),
        sample_id = "{sample_id}",
        purity = _get_sample_purity,
        ploidy = _get_sample_ploidy,
        outdir = CFG["dirs"]["loh"] + "{seq_type}--{genome_build}/{sample_id}/"
    conda: CFG["conda_envs"]["spechla"]
    group: "loh_symlink"
    threads:
        CFG["threads"]["spechla"]
    resources:
        **CFG["resources"]["spechla"] 
    shell:
        op.as_one_line("""
        perl {params.scriptdir}/script/cal.hla.copy.pl 
        -purity {params.purity}
        -ploidy {params.ploidy}
        -S {params.sample_id} 
        -F {input.freq_list}
        -T {input.hla_results}
        -O {params.outdir} > {log.stderr} 2>&1
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _spechla_output_txt:
    input:
        hla_results = CFG["dirs"]["spechla"] + "{seq_type}--{genome_build}/{sample_id}/hla.result.txt",
        loh_info = CFG["dirs"]["loh"] + "{seq_type}--{genome_build}/{sample_id}/merge.hla.copy.txt"
    output:
        hla_results = CFG["dirs"]["outputs"] + "hla_results/{seq_type}--{genome_build}/{sample_id}.hla.result.txt",
        loh_info = CFG["dirs"]["outputs"] + "loh_info/{seq_type}--{genome_build}/{sample_id}.merge.hla.copy.txt"
    group: "loh_symlink"
    run:
        op.relative_symlink(input.hla_results, output.hla_results, in_module= True)
        op.relative_symlink(input.loh_info, output.loh_info, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _spechla_all:
    input:
        expand(
            [
                str(rules._spechla_output_txt.output.hla_results),
                str(rules._spechla_output_txt.output.loh_info),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
