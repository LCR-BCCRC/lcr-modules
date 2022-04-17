# Original Author:  Nicole Thomas
# Module Author:    Nicole Thomas
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

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
# `CFG` is a shortcut to `config["lcr-modules"]["promethion"]`
CFG = op.setup_module(
    name = "nanopolish",
    version = "1.0",
    subdirectories = ["inputs", "meth_calls", "meth_freq","merged_nanopolish_calls","merged_nanopolish_freq","outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _promethion_input,
    _install_nanopolish,
    _promethion_input_chrs,
    _nanopolish_output,
    _nanopolish_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _promethion_input:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"],
        fastq = CFG["inputs"]["sample_fastq"],
        index = CFG["inputs"]["sample_index"],
        fai = CFG["inputs"]["sample_fai"],
        gzi = CFG["inputs"]["sample_gzi"],
        readdb = CFG["inputs"]["sample_readdb"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        fastq = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.fastq",
        index = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.fastq.index",
        fai = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.fastq.index.fai",
        gzi = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.fastq.index.gzi",
        readdb = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.fastq.index.readdb",
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.fastq, output.fastq)
        op.absolute_symlink(input.index, output.index)
        op.absolute_symlink(input.fai, output.fai)
        op.absolute_symlink(input.gzi, output.gzi)
        op.absolute_symlink(input.readdb, output.readdb)

# Installs latest nanopolish release from github if the nanopolish folder is not present yet
rule _install_nanopolish:
    params:
        nanopolish = CFG["inputs"]["nanopolish_exec"]
    output: 
        complete = touch(CFG["inputs"]["nanopolish_exec"] + "/nanopolish_dependencies_installed.success")
    shell:
        op.as_one_line('''
        mkdir -p {params.nanopolish};
        git clone git@github.com:jts/nanopolish.git {params.nanopolish};
        cd {params.nanopolish};
        make;
        cd -;
        
        {output.complete};
        ''')


# Symlink chromosomes used for parallelization
checkpoint _promethion_input_chrs:
    input:
        chrs = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        chrs = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.absolute_symlink(input.chrs, output.chrs)



rule _nanopolish_meth_calls:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fastq = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.fastq",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        installed = str(rules._install_nanopolish.output.complete)
    params:
         nano = CFG["inputs"]["nanopolish_exec"] + "/nanopolish"
    threads:
        CFG["threads"]["nanopolish"] 
    log:
        stderr = CFG["logs"]["meth_calls"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.meth_calls.stderr.log"     
    resources: 
        mem_mb = CFG["mem_mb"]["meth_calls"]          
    output:
        calls = CFG["dirs"]["meth_calls"] + "{seq_type}--{genome_build}/{sample_id}/chromosomes/{chrom}.calls.tsv.gz"    

    shell:
        op.as_one_line('''
            {params.nano} call-methylation -t {threads} 
            -r {input.fastq} -b {input.bam} -g {input.fasta} 
            -q cpg -w {wildcards.chrom} | gzip > {output.calls} 
            2> {log.stderr}
            ''')


rule _nanopolish_meth_freq:
        input: 
            calls = str(rules._nanopolish_meth_calls.output.calls)
        output: 
            freq = CFG["dirs"]["meth_freq"] + "{seq_type}--{genome_build}/{sample_id}/chromosomes/{chrom}.frequency.tsv.gz"
        params:
            nano = CFG["inputs"]["nanopolish_exec"] + "/scripts/calculate_methylation_frequency.py" 
        log:
            stderr = CFG["logs"]["meth_freq"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.meth_freq.stderr.log"
        resources: 
            mem_mb = CFG["mem_mb"]["meth_freq"] 
        shell:
            "{params.nano} -s {input.calls} | gzip > {output.freq} 2> {log.stderr}"


def _nanopolish_get_chr_meth_calls(wildcards):
    CFG = config["lcr-modules"]["nanopolish"]
    chrs = checkpoints._promethion_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    calls = expand(
        CFG["dirs"]["meth_calls"] + "{{seq_type}}--{{genome_build}}/{{sample_id}}/chromosomes/{chrom}.calls.tsv.gz",
        chrom = chrs
    )
    return(calls)


def _nanopolish_get_chr_meth_freq(wildcards):
    CFG = config["lcr-modules"]["nanopolish"]
    chrs = checkpoints._promethion_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    freqs = expand(
        CFG["dirs"]["meth_freq"] + "{{seq_type}}--{{genome_build}}/{{sample_id}}/chromosomes/{chrom}.frequency.tsv.gz",
        chrom = chrs
    )
    return(freqs)


rule _merge_nanopolish_calls:
    input:
        calls = _nanopolish_get_chr_meth_calls
    output:
        calls = CFG["dirs"]["merged_nanopolish_calls"] + "{seq_type}--{genome_build}/{sample_id}.calls.tsv.gz"
    resources: 
        mem_mb = CFG["mem_mb"]["merged_nanopolish_calls"]               
    shell:
        op.as_one_line("""
        zcat {input.calls} | grep -v "chromosome" |
        cat <(zcat {input.calls} | head -n1) - | 
        gzip >  {output.calls} 
        """)


rule _merge_nanopolish_freq:
    input:
        freq = _nanopolish_get_chr_meth_freq
    output:
        freq = CFG["dirs"]["merged_nanopolish_freq"] + "{seq_type}--{genome_build}/{sample_id}.frequency.tsv.gz"
    resources: 
        mem_mb = CFG["mem_mb"]["merged_nanopolish_freq"]   
    shell:
        op.as_one_line("""
        zcat {input.freq} | grep -v "chromosome" |
        cat <(zcat {input.freq} | head -n1) - | 
        gzip >  {output.freq}
        """)


 

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _nanopolish_output:
    input:
        calls = str(rules._merge_nanopolish_calls.output.calls),
        freq = str(rules._merge_nanopolish_freq.output.freq)
    output:
        calls = CFG["dirs"]["outputs"] + "nanopolish_meth_calls/{seq_type}--{genome_build}/{sample_id}.calls.tsv.gz",
        freq = CFG["dirs"]["outputs"] + "nanopolish_meth_freq/{seq_type}--{genome_build}/{sample_id}.frequency.tsv.gz"
    run:
        op.relative_symlink(input.calls, output.calls, in_module= True),
        op.relative_symlink(input.freq, output.freq, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _nanopolish_all:
    input:
        expand(
            [
                str(rules._nanopolish_output.output.calls),
                str(rules._nanopolish_output.output.freq)
                
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)                  





