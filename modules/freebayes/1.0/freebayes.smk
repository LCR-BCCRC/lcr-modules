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
# `CFG` is a shortcut to `config["lcr-modules"]["freebayes"]`
CFG = op.setup_module(
    name = "freebayes",
    version = "1.0",
    subdirectories = ["inputs", "freebayes", "crossmap", "normalize", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _freebayes_input_bam,
    _freebayes_normalize_prefix,
    _freebayes_output_vcf,
    _freebayes_all,

##### CROSSMAP FUNCTIONS #####

# To normalize prefixes, we need to determine
# 1) Which iteration of the genome build something is (i.e. GRch37 or GRCh38)
# 2) If a genome build is chr-prefixed or not
# To determine this, lets load the reference config and re-parse it
FREEBAYES_REFERENCE_CONFIG = CFG["inputs"]["reference_config"]  # Reference config path
configfile: FREEBAYES_REFERENCE_CONFIG
# Store all the attributes we will need
FREEBAYES_GENOME_VERSION = {}  # Will be a simple {"GRCh38-SFU": "grch38"} etc.
FREEBAYES_GENOME_PREFIX = {}  # Will be a simple hash of {"GRCh38-SFU": True} if chr-prefixed
FREEBAYES_VERSION_MAP = {}

for genome_build, attributes in config['genome_builds'].items():
    try:
        genome_version = attributes["version"]
    except KeyError as e:
        # This wasn't included in the reference entry for this genome build
        # This should never happen, as the reference workflow checks for this,
        # but ¯\_(ツ)_/¯
        raise AttributeError(f"Unable to determine the \"version\" of genome {genome_version} in reference config {FREEBAYES_REFERENCE_CONFIG}") from e
    try:
        genome_provider = attributes["provider"]
    except KeyError as e:
        raise AttributeError(f"Unable to determine the \"provider\" of genome {genome_version} in reference config {FREEBAYES_REFERENCE_CONFIG}") from e

    FREEBAYES_GENOME_VERSION[genome_build] = genome_version  # What is the parent genome build?
    FREEBAYES_GENOME_PREFIX[genome_build] = True if genome_provider == "ucsc" else False  # Is this chr-prefixed?
    FREEBAYES_VERSION_MAP[genome_build] = genome_version.replace("grch", "GRCh")  # Genome build for freebayes



##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _freebayes_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"], 
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam", 
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai", 
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    group: 
        "input_and_run"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _freebayes_run:
    input:
        bam = str(rules._freebayes_input_bam.output.bam),
        bai = str(rules._freebayes_input_bam.output.bai), 
        crai = str(rules._freebayes_input_bam.output.crai),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["freebayes"] + "{seq_type}--{genome_build}/{sample_id}/freebayes.vcf"
    log:
        stderr = CFG["logs"]["freebayes"] + "{seq_type}--{genome_build}/{sample_id}/freebayes.stderr.log"
    params:
        opts = CFG["options"]["freebayes"]
    conda:
        CFG["conda_envs"]["freebayes"]
    group: 
        "input_and_run"
    threads:
        CFG["threads"]["freebayes"]
    resources:
        **CFG["resources"]["freebayes"]    # All resources necessary can be included and referenced from the config files.
    shell:
        op.as_one_line("""
            freebayes-parallel 
                <(fasta_generate_regions.py 
                {input.fasta}.fai 100000) 
                {threads} 
                -f {input.fasta} 
                {input.bam} 
                > {output.vcf} 
                2>> {log.stderr} 
        """)

def get_original_genome(wildcards):
    # Determine the original (i.e. input) reference genome for this sample
    # Since this module projects to various output genome builds, we need to parse the sample table for the starting build
    # To determine what we need to do
    sample_table = config['lcr-modules']["freebayes"]["samples"]
    sample_entry = op.filter_samples(sample_table, sample_id = wildcards.sample_id, seq_type = wildcards.seq_type)
    if len(sample_entry) == 0:
        raise AttributeError("Unable to locate a a sample with sample_id:{wildcards.sample_id}, normal_id:{wildcards.normal_id}, seq_type:{wildcards.seq_type} in the \'sample\' table")
    original_genome_build = sample_entry.genome_build.tolist()
    return original_genome_build

def get_chain(genome_build):
    # NOTE: This only currently supports hg38 and hg19. If you are using other genome builds, this will need to be handled
    genome_version = FREEBAYES_GENOME_VERSION[genome_build]
    if genome_version == "grch38":
        return reference_files("genomes/" + genome_build + "/chains/grch38/hg38ToHg19.over.chain")
    elif genome_version == "grch37":
        return reference_files("genomes/" + genome_build +"/chains/grch37/hg19ToHg38.over.chain")
    else:
        raise AttributeError(f"No supported CrossMap chain for {genome_version} within this module")


def crossmap_input(wildcards):
    CFG = config["lcr-modules"]["freebayes"]
    original_genome_build = get_original_genome(wildcards)[0]
    return {
        "vcf":
            expand(
                rules._freebayes_run.output.vcf,
                **wildcards,
                genome_build = original_genome_build
            ),
        "chain": get_chain(original_genome_build)
    }

rule _freebayes_crossmap:
    input:
        unpack(crossmap_input), 
        fasta = reference_files("genomes/{target_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["crossmap"] + "{seq_type}--{target_build}/{sample_id}/freebayes.vcf"
    log:
        stdout = CFG["logs"]["crossmap"] + "{seq_type}--{target_build}/{sample_id}/freebayes.vcf.crossmap.stdout.log",
        stderr = CFG["logs"]["crossmap"] + "{seq_type}--{target_build}/{sample_id}/freebayes.vcf.crossmap.stderr.log"
    conda:
        CFG["conda_envs"]["crossmap"]
    threads:
        CFG["threads"]["crossmap"]
    resources:
        **CFG["resources"]["crossmap"]
    wildcard_constraints:
        target_build = "hg38|hg19"  # Crossmap only converts to chr-prefixed outputs, so these are what will be generated
    shell:
        op.as_one_line("""
        CrossMap vcf 
        {input.chain} 
        {input.vcf} 
        {input.fasta} 
        {output.vcf}
        > {log.stdout} 2> {log.stderr}
        """)

def get_normalize_input(wildcards, todo_only = False):
    CFG = config["lcr-modules"]["freebayes"]
    new_genome_build  = wildcards.target_build
    # Since snakemake only knows what the TARGET genome build is, we need to find the source
    original_genome_builds = get_original_genome(wildcards)
    original_genome_version = [FREEBAYES_GENOME_VERSION[x] for x in original_genome_builds]
    new_genome_version = FREEBAYES_GENOME_VERSION[new_genome_build]

    # Do we need to run CrossMap on this? Check the genome version
    if new_genome_build in original_genome_builds:
        # Source matches. CrossMap not necessary. No vcf normalization needed.
        vcf = expand(
            str(rules._freebayes_run.output.vcf),
            **wildcards,
            genome_build = new_genome_build
        )
        todo = "symlink"
    elif new_genome_version in original_genome_version:
        # Source matches. CrossMap not necessary.
        # Get the list index of the matching item
        index = original_genome_version.index(new_genome_version)
        # Get the input vcf file
        vcf = expand(
            str(rules._freebayes_run.output.vcf),
            **wildcards,
            genome_build = original_genome_builds[index]
        )
        # Check if the chr prefix status is the same e.g. for grch37 vs. hs37d5.
        # If yes, symlink. If no, normalize.
        if FREEBAYES_GENOME_PREFIX[new_genome_build] == FREEBAYES_GENOME_PREFIX[original_genome_builds[index]]:
            todo = "symlink"
        else:
            todo = "normalize"
    else:
        # Source doesn't match. CrossMap and normalization necessary.
        target_build = "hg38" if new_genome_version == "grch38" else "hg19"
        vcf = expand(
            rules._freebayes_crossmap.output.vcf,
            target_build = target_build,
            allow_missing = True
        )
        todo = "normalize"

    if todo_only:
        return todo
    else:
        return {"vcf": vcf}

# Add or remove chr prefix as necessary
rule _freebayes_normalize_prefix:
    input:
        unpack(lambda w: get_normalize_input(w, todo_only = False))
    output:
        vcf = CFG["dirs"]["normalize"] + "{seq_type}--{target_build}/{sample_id}/freebayes.vcf"
    params:
        todo = lambda w: get_normalize_input(w, todo_only = True),
        dest_chr = lambda w: FREEBAYES_GENOME_PREFIX[w.target_build]
    group: "normalize_and_bgzip"
    wildcard_constraints:
        target_build = "|".join(CFG["options"]["target_builds"])
    run:
        logger.info("Handling " + str(wildcards.sample_id) + " for " + str(wildcards.target_build) + ":\nparams.dest_chr: " + str(params.dest_chr) + "\ntodo: " + str(params.todo))
        if params.todo == "symlink":
            op.relative_symlink(input.vcf, output.vcf, in_module = True)
        else:
            logger.info("Normalizing " + input.vcf[0] + " with params.dest_chr=" + str(params.dest_chr))
            with open(input.vcf[0]) as f, open(output.vcf, "w") as o:
                for line in f: 
                    line = line.rstrip("\n").rstrip("\r")
                    # To handle CrossMap weirdness, remove all chr-prefixes and add them back later
                    if line.startswith("chr"): # only modify the text if it starts with the prefix
                        line = line.replace("chr", "", 1)
    
                    if line.startswith("##contig=<ID="): # Remove chr prefixes from header contig lines
                        line = line.replace("chr", "", 1)
                            
                    if params.dest_chr and not line.startswith("#"): # Will evaluate to True if the destination genome is chr-prefixed
                        # Add chr prefix
                        line = "chr" + str(line)

                    if params.dest_chr and line.startswith("##contig=<ID="): # Fix chr prefixing of header contig lines
                        line = line.replace("##contig=<ID=", "##contig=<ID=chr", 1)
                    o.write(line)
                    o.write(os.linesep)

rule _freebayes_gzip_vcf: 
    input: 
        unpack(lambda w: get_normalize_input(w, todo_only = False)), # Ensure temp vcf files don't get deleted
        invcf = str(rules._freebayes_normalize_prefix.output.vcf)
    output: 
        vcf = CFG["dirs"]["normalize"] + "{seq_type}--{target_build}/{sample_id}/freebayes.vcf.gz", 
        tbi = CFG["dirs"]["normalize"] + "{seq_type}--{target_build}/{sample_id}/freebayes.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["crossmap"]
    threads:
        CFG["threads"]["gzip"]
    resources:
        **CFG["resources"]["gzip"] 
    group: "normalize_and_bgzip"
    shell: 
        op.as_one_line("""
            bcftools sort -Oz -o {output.vcf} {input.invcf} && 
            tabix -p vcf {output.vcf}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _freebayes_output_vcf:
    input:
        vcf = str(rules._freebayes_gzip_vcf.output.vcf), 
        tbi = str(rules._freebayes_gzip_vcf.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--projection/{sample_id}.freebayes.{target_build}.vcf.gz", 
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--projection/{sample_id}.freebayes.{target_build}.vcf.gz.tbi"
    group: "normalize_and_bgzip"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True)
        op.relative_symlink(input.tbi, output.tbi, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _freebayes_all:
    input:
        expand(
            expand(
                rules._freebayes_output_vcf.output,
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["samples"]["seq_type"],
                sample_id=CFG["samples"]["sample_id"], 
                allow_missing = True
            ), 
            target_build = CFG["options"]["target_builds"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
