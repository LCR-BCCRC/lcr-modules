#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Jasper
# Module Author:    Jasper
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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["controlfreec"]`
CFG = op.setup_module(
    name = "controlfreec",
    version = "1.2",
    subdirectories = ["inputs", "mpileup", "run", "calc_sig_and_plot",
    "freec2bed", "freec2circos", "cnv2igv", "convert_coordinates",
    "fill_regions", "normalize", "outputs"]
)



# Define rules to be run locally when using a compute cluster
localrules:
    _controlfreec_get_map_refs_hg19,
    _controlfreec_get_map_refs_hg38,
    _controlfreec_symlink_map_refs_unmasked,
    _controlfreec_symlink_map,
    _controlfreec_input_chroms,
    _controlfreec_check_chrFiles,
    _controlfreec_input_bam,
    _controlfreec_symlink_run_result,
    _controlfreec_plot,
    _controlfreec_output,
    _controlfreec_all


##### RULES #####

#### Rules for mappability reference
# when NOT using masked reference genome
rule _controlfreec_get_map_refs_hg19:
    output:
        tar = temp(CFG["dirs"]["inputs"] + "references/mappability/unmasked/raw/out100m2_hg19.tar.gz"),
        gem = CFG["dirs"]["inputs"] + "references/mappability/unmasked/out100m2_hg19.gem"
    params:
        outdir = CFG["dirs"]["inputs"] + "references/mappability/unmasked/"
    conda:
        CFG["conda_envs"]["wget"]
    shell:
        """
        wget -O {output.tar} http://xfer.curie.fr/get/7hZIk1C63h0/hg19_len100bp.tar.gz && \
        tar -xvf {output.tar} -C {params.outdir} --wildcards --no-anchored 'out100m2*gem'
        """
rule _controlfreec_get_map_refs_hg38:
    output:
        tar = temp(CFG["dirs"]["inputs"] + "references/mappability/unmasked/raw/out100m2_hg38.zip"),
        gem = CFG["dirs"]["inputs"] + "references/mappability/unmasked/raw/out100m2_hg38.gem"
    params:
        outdir = CFG["dirs"]["inputs"] + "references/mappability/unmasked/"
    conda:
        CFG["conda_envs"]["wget"]
    shell:
        """
        wget -O {output.tar} http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip && \
        unzip {output.tar} -d {params.outdir}
        """

def _get_unmasked_map_ref(wildcards):
    if "38" in str({wildcards.genome_build}):
        return str(rules._controlfreec_get_map_refs_hg38.output.gem)
    else:
        return str(rules._controlfreec_get_map_refs_hg19.output.gem)


rule _controlfreec_symlink_map_refs_unmasked:
    input:
        mappability = _get_unmasked_map_ref
    output:
        mappability = CFG["dirs"]["inputs"] + "references/mappability/unmasked/out100m2_{genome_build}.gem"
    run:
        op.absolute_symlink(input.mappability, output.mappability)

# mappability tracks for hard-masked genomes are generated using GEM
rule _controlfreec_download_GEM:
    output:
        touch(CFG["dirs"]["inputs"] + "references/GEM/.done")
    params:
        dirOut = CFG["dirs"]["inputs"] + "references/GEM"
    conda:
        CFG["conda_envs"]["wget"]
    threads: 1
    shell:
        op.as_one_line("""
        wget https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download -O {params.dirOut}/GEM-lib.tbz2 &&
        bzip2 -dc {params.dirOut}/GEM-lib.tbz2 | tar -xvf - -C {params.dirOut}/
        """)


# grch37_masked and grch38_masked from ensembl have additional information in the chr name lines - need to remove
rule _controlfreec_fix_grch_genomes:
    input:
        reference = reference_files("genomes/{genome_build}_masked/genome_fasta/genome.fa")
    output:
        reference = CFG["dirs"]["inputs"] + "references/mappability/masked/grch_fasta_fixes/{genome_build}_fix.fa"
    threads: 1
    resources: **CFG["resources"]["gem"]
    shell:
        "cat {input.reference} | perl -ne 's/(^\>\S+).+/$1/;print;' > {output.reference} "

# the only masked ref genomes available currently are grch37_masked, grch38_masked, hg19_masked, hg38_masked, and hg19-reddy_masked
# hs37d5 can use grch37 files, since it is grch37+decoys and the decoys aren't used in CN calling
def _get_genome_fasta(wildcards):
    if  "grch37" in str({wildcards.genome_build}) or "hs37d5" in str({wildcards.genome_build}): # covers both cases "grch37" and ones like "grch37-noalt"
        return  str(rules._controlfreec_fix_grch_genomes.output.reference).replace("{genome_build}", "grch37")d
    elif "grch38" in str({wildcards.genome_build}): # covers both cases "grch38" and ones like "grch38-nci"
        return  str(rules._controlfreec_fix_grch_genomes.output.reference).replace("{genome_build}", "grch38")
    elif "hg38" in str({wildcards.genome_build}): # covers both cases "hg38" and ones like "hg38-nci"
        return reference_files("genomes/hg38_masked/genome_fasta/genome.fa")
    else: # must be a flavour of hg19
        return reference_files("genomes/hg19_masked/genome_fasta/genome.fa")

rule _controlfreec_generate_gem_index:
    input:
        software = str(rules._controlfreec_download_GEM.output),
        reference = _get_genome_fasta
    output:
        index = CFG["dirs"]["inputs"] + "references/mappability/masked/{genome_build}.hardmask.all_index.gem"
    params:
        gemDir = CFG["dirs"]["inputs"] + "references/GEM/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin",
        idxpref = CFG["dirs"]["inputs"] + "references/mappability/masked/{genome_build}.hardmask.all_index"
    threads: CFG["threads"]["gem"]
    resources: **CFG["resources"]["gem"]
    log: CFG["logs"]["inputs"] + "gem/{genome_build}/gem_index.stderr.log"
    shell:
        op.as_one_line("""
        PATH=$PATH:{params.gemDir};
        {params.gemDir}/gem-indexer -T {threads} -c dna -i {input.reference} -o {params.idxpref} > {log} 2>&1
        """)

rule _controlfreec_generate_mappability:
    input:
        software = str(rules._controlfreec_download_GEM.output),
        index = str(rules._controlfreec_generate_gem_index.output.index)
    output:
        mappability = CFG["dirs"]["inputs"] + "references/mappability/masked/{genome_build}.hardmask.all.gem.mappability"
    params:
        gemDir = CFG["dirs"]["inputs"] + "references/GEM/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin",
        pref =  CFG["dirs"]["inputs"] + "references/mappability/masked/{genome_build}.hardmask.all.gem",
        kmer = CFG["options"]["kmer"],
        mismatch = CFG["options"]["mismatch"],
        maxEditDistance = CFG["options"]["maxEditDistance"],
        maxBigIndel = CFG["options"]["maxBigIndel"],
        strata = CFG["options"]["strata"]
    threads: CFG["threads"]["gem"]
    resources: **CFG["resources"]["gem"]
    log: CFG["logs"]["inputs"] + "gem/{genome_build}/gem_map.stderr.log"
    shell:
        op.as_one_line("""
        PATH=$PATH:{params.gemDir};
        {params.gemDir}/gem-mappability
            -T {threads}
            -I {input.index}
            -l {params.kmer}
            -m {params.mismatch}
            -t disable
            --mismatch-alphabet ACGNT
            -e {params.maxEditDistance}
            --max-big-indel-length {params.maxBigIndel}
            -s {params.strata}
            -o {params.pref} > {log} 2>&1
        """)

rule _controlfreec_symlink_map:
    input:
        mappability = str(rules._controlfreec_generate_mappability.output.mappability)
    output:
        mappability = CFG["dirs"]["inputs"] + "references/mappability/masked/out100m2_{genome_build}.gem"
    run:
        op.absolute_symlink(input.mappability, output.mappability)


#### Rules for genome reference files needed
# these do not need to be done using the masked reference fasta

# generates file with chromomsome lengths from genome.fa.fai
rule _controlfreec_generate_chrLen:
    input:
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
    output:
        chrLen = CFG["dirs"]["inputs"] + "references/{genome_build}/{genome_build}.len"
    threads: 1
    resources: **CFG["resources"]["gem"]
    shell:
        op.as_one_line("""
        grep -P '^chr[0-9,X,Y]+\t|^[0-9,X,Y]' {input.fai} | awk '{{print $1"\t"$2}}' > {output.chrLen}
        """)

# generates a per-chromosome fasta file
checkpoint _controlfreec_input_chroms:
    input:
        txt = ancient(reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt"))
    output:
        txt = CFG["dirs"]["inputs"] + "references/{genome_build}/main_chromosomes_withY.txt"
    run:
        op.absolute_symlink(input.txt, output.txt)

rule _controlfreec_generate_chrFasta:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        fasta = CFG["dirs"]["inputs"] + "references/{genome_build}/chr/{chromosome}.fa"
    conda:
        CFG["conda_envs"]["controlfreec"]
    threads: CFG["threads"]["gem"]
    resources: **CFG["resources"]["gem"]
    shell:
        "samtools faidx {input.fasta} {wildcards.chromosome} > {output.fasta} "

def _get_chr_fastas(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    with open(checkpoints._controlfreec_input_chroms.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")

    fastas = expand(
        CFG["dirs"]["inputs"] +  "references/{{genome_build}}/chr/{chromosome}.fa",
        chromosome = chrs
        )
    return fastas

# checks that all chroms are accounted for
rule _controlfreec_check_chrFiles:
    input:
        _get_chr_fastas,
        main = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt")
    output:
        touch(CFG["dirs"]["inputs"] + "references/{genome_build}/chr/.all_done")


#### Rule for reference dbsnp file needed for mpileup step
rule _controlfreec_dbsnp_to_bed:
    input:
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        bed = CFG["dirs"]["inputs"] + "references/variation/{genome_build}/dbsnp.common_all-151.bed"
    threads: 1
    resources: **CFG["resources"]["gem"]
    shell:
        op.as_one_line("""
        gunzip -c {input.vcf} | awk {{'printf ("%s\\t%s\\t%s\\n", $1,$2-1,$2)'}} | zgrep -v -h "^#" > {output.bed}
        """)

#### Rules for getting inputs related to the bam files
# Symlinks the input files into the module results directory (under '00-inputs/')
checkpoint _controlfreec_input_bam:
    input:
        bam = ancient(CFG["inputs"]["sample_bam"]),
        bai = ancient(CFG["inputs"]["sample_bai"])
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bai",
        crai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)

rule _controlfreec_mpileup_per_chrom:
    input:
        bam = str(rules._controlfreec_input_bam.output.bam),
        fastaFile = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        bed = str(rules._controlfreec_dbsnp_to_bed.output.bed)
    output:
        pileup = temp(CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.minipileup.pileup.gz")
    conda:
        CFG["conda_envs"]["controlfreec"]
    threads: 1
    resources:
        **CFG["resources"]["mpileup"]
    group: "mpileup_controlfreec"
    log:
        stderr = CFG["logs"]["inputs"] + "mpileup/{seq_type}--{genome_build}/{sample_id}.{chrom}.mpileup.stderr.log",
    shell:
        "samtools mpileup -l {input.bed} -r {wildcards.chrom} -Q 20 -f {input.fastaFile} {input.bam} | gzip -c > {output.pileup} 2> {log.stderr}"

# set-up mpileups per chromosome for BAF calling
def _get_chr_mpileups(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    with open(checkpoints._controlfreec_input_chroms.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")

    mpileups = expand(
        [str(rules._controlfreec_mpileup_per_chrom.output.pileup)],
        chrom = chrs,
        **wildcards
    )
    return mpileups

rule _controlfreec_concatenate_pileups:
    input:
        mpileup = _get_chr_mpileups,
        main = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt")
    output:
        mpileup = temp(CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.bam_minipileup.pileup.gz")
    threads: 1
    resources:
        **CFG["resources"]["cat"]
    group: "mpileup_controlfreec"
    shell:
        "cat {input.mpileup} > {output.mpileup} "


#### Run control-FREEC ####

# set-up controlfreec config files
# these have to be each their own function, bc `params:` can't handle dicts or use unpack()
def _controlfreec_get_optional_step(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    if CFG["options"]["step"] != "":
        cfg_step = "step = " + str(CFG["options"]["step"]) # will replace the entire line, so it's no longer commented out
    else:
        cfg_step = "#step = " + str(CFG["options"]["step"]) # will stay commented out, won't be used, just needed for the code below
    return cfg_step

def _controlfreec_get_optional_window(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    if CFG["options"]["window"] != "":
        cfg_window = "window = " + str(CFG["options"]["window"]) # will replace the entire line, so it's no longer commented out
    else:
        cfg_window = "#window = " + str(CFG["options"]["window"]) # will stay commented out, won't be used, just needed for the code below
    return cfg_window

rule _controlfreec_config_contamAdjTrue:
    input:
        tumour_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.bam_minipileup.pileup.gz",
        normal_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.bam_minipileup.pileup.gz",
        mappability = CFG["dirs"]["inputs"] + "references/mappability/{masked}/out100m2_{genome_build}.gem",
        chrLen = str(rules._controlfreec_generate_chrLen.output.chrLen),
        done = str(rules._controlfreec_check_chrFiles.output),
        dbsnp = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        config = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/contamAdjTrue/config.txt"
    conda:
        CFG["conda_envs"]["controlfreec"]
    threads: 1
    params:
        config = CFG["options"]["configFile"],
        window = _controlfreec_get_optional_window,
        step = _controlfreec_get_optional_step,
        booCon = "TRUE",
        numCon = CFG["options"]["contamination"],
        bedGraphOutput = CFG["options"]["BedGraphOutput"],
        breakPointValue = CFG["options"]["breakPointThreshold"],
        breakPointType = CFG["options"]["breakPointType"],
        coefVar = CFG["options"]["coefficientOfVariation"],
        degree = CFG["options"]["degree"],
        forceGC = CFG["options"]["forceGCcontentNormalization"],
        intercept = CFG["options"]["intercept"],
        minCNAlength = CFG["options"]["minCNAlength"],
        minMapPerWindow = CFG["options"]["minMappabilityPerWindow"],
        minimumSubclonePresence = CFG["options"]["minimalSubclonePresence"],
        noisyData = CFG["options"]["noisyData"],
        readCountThreshold = CFG["options"]["readCountThreshold"],
        ploidy = CFG["options"]["ploidy"],
        naBoo = CFG["options"]["printNA"],
        telocentromeric = CFG["options"]["telocentromeric"],
        uniqBoo = CFG["options"]["uniqueMatch"],
        minimalCoveragePerPosition = CFG["options"]["minimalCoveragePerPosition"],
        minimalQualityPerPosition = CFG["options"]["minQualityPerPosition"],
        shiftInQuality = CFG["options"]["shiftInQuality"],
        threads = CFG["threads"]["controlfreec_run"],
        outdir = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/contamAdjTrue/",
        chrFiles = CFG["dirs"]["inputs"] + "references/{genome_build}/chr/"
    shell:
        "samtoolspath=$(which samtools ) ; "
        "samtoolsPathName=$(echo $samtoolspath) ; "
        "sambambapath=$(which sambamba ) ; "
        "sambambaPathName=$(echo $sambambapath) ; "
        "bedtoolspath=$(which bedtools ) ; "
        "bedtoolsPathName=$(echo $bedtoolspath) ; "
        "sed \"s|BAMFILE|{input.tumour_bam}|g\" {params.config} | "
        "sed \"s|CONTROLFILE|{input.normal_bam}|g\" | "
        "sed \"s|OUTDIR|{params.outdir}|g\" | "
        "sed \"s|DBsnpFile|{input.dbsnp}|g\" | "
        "sed \"s|phredQuality|{params.shiftInQuality}|g\" | "
        "sed \"s|ploidyInput|{params.ploidy}|g\" | "
        "sed \"s|chrLenFileInput|{input.chrLen}|g\" | "
        "sed \"s|chrFilesPath|{params.chrFiles}|g\" | "
        "sed \"s|sambambaPath|$sambambaPathName|g\" | "
        "sed \"s|bedtoolsPath|$bedtoolsPathName|g\" | "
        "sed \"s|samtoolsPath|$samtoolsPathName|g\" | "
        "sed \"s|breakPointValue|{params.breakPointValue}|g\" | "
        "sed \"s|breakPointTypeBe|{params.breakPointType}|g\" | "
        "sed \"s|bedgraphBoo|{params.bedGraphOutput}|g\" | "
        "sed \"s|coefVar|{params.coefVar}|g\" | "
        "sed \"s|#step = stepValue|{params.step}|g\" | "
        "sed \"s|#window = windowSize|{params.window}|g\" | "
        "sed \"s|booCon|{params.booCon}|g\" | "
        "sed \"s|#contamination = numCon|{params.numCon}|g\" | "
        "sed \"s|forceGCvalue|{params.forceGC}|g\" | "
        "sed \"s|interPoly|{params.intercept}|g\" | "
        "sed \"s|numDegree|{params.degree}|g\" | "
        "sed \"s|minCNAvalue|{params.minCNAlength}|g\" | "
        "sed \"s|minMapPerWindow|{params.minMapPerWindow}|g\" | "
        "sed \"s|minimumSubclonePresenceValue|{params.minimumSubclonePresence}|g\" | "
        "sed \"s|minCovPerPos|{params.minimalCoveragePerPosition}|g\" | "
        "sed \"s|minQualPerPos|{params.minimalQualityPerPosition}|g\" | "
        "sed \"s|booNoise|{params.noisyData}|g\" | "
        "sed \"s|rcCountThresold|{params.readCountThreshold}|g\" | "
        "sed \"s|teloValue|{params.telocentromeric}|g\" | "
        "sed \"s|uniqBoo|{params.uniqBoo}|g\" | "
        "sed \"s|naBoo|{params.naBoo}|g\" | "
        "sed \"s|numThreads|{params.threads}|g\" | "
        "sed \"s|referenceFile|{input.mappability}|g\" > {output.config}"


rule _controlfreec_config_contamAdjFalse:
    input:
        tumour_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.bam_minipileup.pileup.gz",
        normal_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.bam_minipileup.pileup.gz",
        mappability = CFG["dirs"]["inputs"] + "references/mappability/{masked}/out100m2_{genome_build}.gem",
        chrLen = str(rules._controlfreec_generate_chrLen.output.chrLen),
        done = str(rules._controlfreec_check_chrFiles.output),
        dbsnp = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        config = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/contamAdjFalse/config.txt"
    conda:
        CFG["conda_envs"]["controlfreec"]
    threads: 1
    params:
        config = CFG["options"]["configFile"],
        window = _controlfreec_get_optional_window,
        step = _controlfreec_get_optional_step,
        booCon = "FALSE",
        numCon = CFG["options"]["contamination"],
        bedGraphOutput = CFG["options"]["BedGraphOutput"],
        breakPointValue = CFG["options"]["breakPointThreshold"],
        breakPointType = CFG["options"]["breakPointType"],
        coefVar = CFG["options"]["coefficientOfVariation"],
        degree = CFG["options"]["degree"],
        forceGC = CFG["options"]["forceGCcontentNormalization"],
        intercept = CFG["options"]["intercept"],
        minCNAlength = CFG["options"]["minCNAlength"],
        minMapPerWindow = CFG["options"]["minMappabilityPerWindow"],
        minimumSubclonePresence = CFG["options"]["minimalSubclonePresence"],
        noisyData = CFG["options"]["noisyData"],
        readCountThreshold = CFG["options"]["readCountThreshold"],
        ploidy = CFG["options"]["ploidy"],
        naBoo = CFG["options"]["printNA"],
        telocentromeric = CFG["options"]["telocentromeric"],
        uniqBoo = CFG["options"]["uniqueMatch"],
        minimalCoveragePerPosition = CFG["options"]["minimalCoveragePerPosition"],
        minimalQualityPerPosition = CFG["options"]["minQualityPerPosition"],
        shiftInQuality = CFG["options"]["shiftInQuality"],
        threads = CFG["threads"]["controlfreec_run"],
        outdir = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/contamAdjFalse/",
        chrFiles = CFG["dirs"]["inputs"] + "references/{genome_build}/chr/"
    shell:
        "samtoolspath=$(which samtools ) ; "
        "samtoolsPathName=$(echo $samtoolspath) ; "
        "sambambapath=$(which sambamba ) ; "
        "sambambaPathName=$(echo $sambambapath) ; "
        "bedtoolspath=$(which bedtools ) ; "
        "bedtoolsPathName=$(echo $bedtoolspath) ; "
        "sed \"s|BAMFILE|{input.tumour_bam}|g\" {params.config} | "
        "sed \"s|CONTROLFILE|{input.normal_bam}|g\" | "
        "sed \"s|OUTDIR|{params.outdir}|g\" | "
        "sed \"s|DBsnpFile|{input.dbsnp}|g\" | "
        "sed \"s|phredQuality|{params.shiftInQuality}|g\" | "
        "sed \"s|ploidyInput|{params.ploidy}|g\" | "
        "sed \"s|chrLenFileInput|{input.chrLen}|g\" | "
        "sed \"s|chrFilesPath|{params.chrFiles}|g\" | "
        "sed \"s|sambambaPath|$sambambaPathName|g\" | "
        "sed \"s|bedtoolsPath|$bedtoolsPathName|g\" | "
        "sed \"s|samtoolsPath|$samtoolsPathName|g\" | "
        "sed \"s|breakPointValue|{params.breakPointValue}|g\" | "
        "sed \"s|breakPointTypeBe|{params.breakPointType}|g\" | "
        "sed \"s|bedgraphBoo|{params.bedGraphOutput}|g\" | "
        "sed \"s|coefVar|{params.coefVar}|g\" | "
        "sed \"s|#step = stepValue|{params.step}|g\" | "
        "sed \"s|#window = windowSize|{params.window}|g\" | "
        "sed \"s|booCon|{params.booCon}|g\" | "
        "sed \"s|#contamination = numCon|{params.numCon}|g\" | "
        "sed \"s|forceGCvalue|{params.forceGC}|g\" | "
        "sed \"s|interPoly|{params.intercept}|g\" | "
        "sed \"s|numDegree|{params.degree}|g\" | "
        "sed \"s|minCNAvalue|{params.minCNAlength}|g\" | "
        "sed \"s|minMapPerWindow|{params.minMapPerWindow}|g\" | "
        "sed \"s|minimumSubclonePresenceValue|{params.minimumSubclonePresence}|g\" | "
        "sed \"s|minCovPerPos|{params.minimalCoveragePerPosition}|g\" | "
        "sed \"s|minQualPerPos|{params.minimalQualityPerPosition}|g\" | "
        "sed \"s|booNoise|{params.noisyData}|g\" | "
        "sed \"s|rcCountThresold|{params.readCountThreshold}|g\" | "
        "sed \"s|teloValue|{params.telocentromeric}|g\" | "
        "sed \"s|uniqBoo|{params.uniqBoo}|g\" | "
        "sed \"s|naBoo|{params.naBoo}|g\" | "
        "sed \"s|numThreads|{params.threads}|g\" | "
        "sed \"s|referenceFile|{input.mappability}|g\" > {output.config}"

# Tries to run with contamAdj = True case first; if it fails then it tries to
# run contamAdj = False; if that also fails, then snakemake receives an error code
# for the rule and the done file isn't made
checkpoint _controlfreec_run:
    input:
        config_contamTrue = str(rules._controlfreec_config_contamAdjTrue.output.config),
        config_contamFalse = str(rules._controlfreec_config_contamAdjFalse.output.config),
        tumour_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.bam_minipileup.pileup.gz",
        normal_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.bam_minipileup.pileup.gz"
    output:
        done = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/done"
    conda:
        CFG["conda_envs"]["controlfreec"]
    threads:
        CFG["threads"]["controlfreec_run"]
    resources:
        **CFG["resources"]["controlfreec_run"]
    log:
        log_contamTrue = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/run.contamTrue.log",
        log_contamFalse = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/run.contamFalse.log"
    shell:
        """
        set +e
        freec -conf {input.config_contamTrue} &> {log.log_contamTrue}
        if [[ $? -ne 0 ]]; then
            freec -conf {input.config_contamFalse} &> {log.log_contamFalse}
            if [[ $? -ne 0 ]]; then
                exit 1
            else
                touch {output.done}
            fi
        else
            touch {output.done}
        fi
        """

def _get_run_result(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    base_dir = os.path.dirname(str(checkpoints._controlfreec_run.get(**wildcards).output[0]))
    CONTAM, = glob_wildcards(base_dir + "/{contamAdj}/{tumour_id}.bam_minipileup.pileup.gz_CNVs").contamAdj
    CONTAM = [CONTAM] if isinstance(CONTAM, str) else CONTAM # makes it a list when only one value is returned by the glob
    if any("True" in c for c in CONTAM):
        files = {
            "info":base_dir + "/contamAdjTrue/{tumour_id}.bam_minipileup.pileup.gz_info.txt",
            "ratios":base_dir + "/contamAdjTrue/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt",
            "CNV":base_dir + "/contamAdjTrue/{tumour_id}.bam_minipileup.pileup.gz_CNVs",
            "BAF":base_dir + "/contamAdjTrue/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt"
        }
    else:
        files = {
            "info":base_dir + "/contamAdjFalse/{tumour_id}.bam_minipileup.pileup.gz_info.txt",
            "ratios":base_dir + "/contamAdjFalse/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt",
            "CNV":base_dir + "/contamAdjFalse/{tumour_id}.bam_minipileup.pileup.gz_CNVs",
            "BAF":base_dir + "/contamAdjFalse/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt"
        }
    return files

rule _controlfreec_symlink_run_result:
    input:
        unpack(_get_run_result)
    output:
        info = CFG["dirs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_info.txt",
        ratios = CFG["dirs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt",
        CNV = CFG["dirs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_CNVs",
        BAF = CFG["dirs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt"
    run:
        op.relative_symlink(input.info, output.info, in_module = True)
        op.relative_symlink(input.ratios, output.ratios, in_module = True)
        op.relative_symlink(input.CNV, output.CNV, in_module = True)
        op.relative_symlink(input.BAF, output.BAF, in_module = True)

rule _controlfreec_calc_sig:
    input:
        info = str(rules._controlfreec_symlink_run_result.output.info),
        ratios = str(rules._controlfreec_symlink_run_result.output.ratios),
        CNV = str(rules._controlfreec_symlink_run_result.output.CNV),
        BAF = str(rules._controlfreec_symlink_run_result.output.BAF)
    output:
        txt = CFG["dirs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_CNVs.p.value.txt"
    params:
        calc_sig = CFG["software"]["FREEC_sig"]
    threads:
        CFG["threads"]["calc_sig"]
    resources:
        **CFG["resources"]["calc_sig"]
    conda:
        CFG["conda_envs"]["controlfreec"]
    log:
        CFG["logs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/calc_sig.log"
    shell:
        "Rscript --vanilla {params.calc_sig} {input.CNV} {input.ratios} {output.txt} &> {log}"


rule _controlfreec_plot:
    input:
        info = str(rules._controlfreec_symlink_run_result.output.info),
        ratios = str(rules._controlfreec_symlink_run_result.output.ratios),
        CNV = str(rules._controlfreec_symlink_run_result.output.CNV),
        BAF = str(rules._controlfreec_symlink_run_result.output.BAF)
    output:
        plot = CFG["dirs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt.png",
        log2plot = CFG["dirs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt.log2.png",
        bafplot = CFG["dirs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt.png"
    params:
        plot = CFG["software"]["FREEC_graph"]
    threads: 1
    resources:
        **CFG["resources"]["plot"]
    conda:
        CFG["conda_envs"]["controlfreec"]
    log:
        CFG["logs"]["calc_sig_and_plot"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/plot.log"
    shell:
        "Rscript --vanilla {params.plot} `grep \"Output_Ploidy\" {input.info} | cut -f 2` {input.ratios} {input.BAF} &> {log}"


rule _controlfreec_freec2bed:
    input:
        info = str(rules._controlfreec_symlink_run_result.output.info),
        ratios = str(rules._controlfreec_symlink_run_result.output.ratios),
        CNV = str(rules._controlfreec_symlink_run_result.output.CNV),
        BAF = str(rules._controlfreec_symlink_run_result.output.BAF)
    output:
        bed = CFG["dirs"]["freec2bed"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bed"
    params:
        freec2bed = CFG["software"]["freec2bed"]
    threads:
        CFG["threads"]["freec2bed"]
    resources:
        **CFG["resources"]["freec2bed"]
    conda:
        CFG["conda_envs"]["controlfreec"]
    log:
        stderr = CFG["logs"]["freec2bed"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/freec2bed.stderr.log"
    shell:
        "ploidy=$(grep Output_Ploidy {input.info} | cut -f 2); "
        "perl {params.freec2bed} -f {input.ratios} -p $ploidy > {output.bed} 2> {log.stderr}"


rule _controlfreec_freec2circos:
    input:
        info = str(rules._controlfreec_symlink_run_result.output.info),
        ratios = str(rules._controlfreec_symlink_run_result.output.ratios),
        CNV = str(rules._controlfreec_symlink_run_result.output.CNV),
        BAF = str(rules._controlfreec_symlink_run_result.output.BAF)
    output:
        circos = CFG["dirs"]["freec2circos"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.circos.bed"
    params:
        freec2circos = CFG["software"]["freec2circos"]
    threads:
        CFG["threads"]["freec2circos"]
    resources:
        **CFG["resources"]["freec2circos"]
    conda:
        CFG["conda_envs"]["controlfreec"]
    log:
        stderr = CFG["logs"]["freec2circos"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/freec2circos.stderr.log"
    shell:
        "ploidy=$(grep Output_Ploidy {input.info} | cut -f 2); "
        "perl {params.freec2circos} -f {input.ratios} -p $ploidy > {output.circos} 2> {log.stderr}"


rule _controlfreec_cnv2igv:
    input:
        cnv = str(rules._controlfreec_calc_sig.output.txt),
        cnv2igv = ancient(CFG["inputs"]["cnv2igv"])
    output:
        seg = CFG["dirs"]["cnv2igv"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.CNVs.seg"
    params:
        opts = CFG["options"]["preserve"]
    conda:
        CFG["conda_envs"]["cnv2igv"]
    threads: 1
    log:
        stderr = CFG["logs"]["cnv2igv"] + "{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}/cnv2igv.stderr.log"
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id} on $(hostname) at $(date)" > {log.stderr};
        python {input.cnv2igv} --mode controlfreec {params.opts} --sample {wildcards.tumour_id}
        {input.cnv} > {output.seg} 2>> {log.stderr}
        """)


def _controlfreec_get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")

# Convert the coordinates of seg file to a different genome build
rule _controlfreec_convert_coordinates:
    input:
        controlfreec_native = str(rules._controlfreec_cnv2igv.output.seg),
        controlfreec_chain = _controlfreec_get_chain
    output:
        controlfreec_lifted = CFG["dirs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.stderr.log"
    threads: 1
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        liftover_minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.liftover_script}
        SEG
        {input.controlfreec_native}
        {output.controlfreec_lifted}
        {input.controlfreec_chain}
        YES
        {params.liftover_minmatch}
        2>> {log.stderr}
        """)

# ensure to request the correct files for each projection and drop wildcards that won't be used downstream
def _controlfreec_prepare_projection(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    tbl = CFG["runs"]
    this_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"].tolist()

    if "38" in this_genome_build[0]:
        hg38_projection = str(rules._controlfreec_cnv2igv.output.seg).replace("{genome_build}", this_genome_build[0])
        grch37_projection = str(rules._controlfreec_convert_coordinates.output.controlfreec_lifted).replace("{genome_build}", this_genome_build[0])
        grch37_projection = grch37_projection.replace("{chain}", "hg38ToHg19")
    else:
        grch37_projection = str(rules._controlfreec_cnv2igv.output.seg).replace("{genome_build}", this_genome_build[0])
        hg38_projection = str(rules._controlfreec_convert_coordinates.output.controlfreec_lifted).replace("{genome_build}", this_genome_build[0])
        hg38_projection = hg38_projection.replace("{chain}", "hg19ToHg38")
    return{
        "grch37_projection": grch37_projection,
        "hg38_projection": hg38_projection
    }


# Fill the missing segments of seg files with neutral regions to complete the genome coverage
rule _controlfreec_fill_segments:
    input:
        unpack(_controlfreec_prepare_projection)
    output:
        grch37_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{masked}/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"),
        hg38_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{masked}/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg")
    log:
        stderr = CFG["logs"]["fill_regions"] + "{seq_type}--projection/{masked}/{tumour_id}--{normal_id}--{pair_status}.{tool}_fill_segments.stderr.log"
    threads: 1
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"]
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        echo "Filling grch37 projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.grch37.bed
        {input.grch37_projection}
        {params.path}src/blacklisted.grch37.bed
        {output.grch37_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        echo "Filling hg38 projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.hg38.bed
        {input.hg38_projection}
        {params.path}src/blacklisted.hg38.bed
        {output.hg38_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        """)

# Normalize chr prefix of the output file
rule _controlfreec_normalize_projection:
    input:
        filled = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{masked}/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg",
        chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
    output:
        projection = CFG["dirs"]["normalize"] + "seg/{seq_type}--projection/{masked}/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    resources:
        **CFG["resources"]["post_controlfreec"]
    threads: 1
    run:
        # read the main chromosomes file of the projection
        chromosomes = pd.read_csv(input.chrom_file, sep = "\t", names=["chromosome"], header=None)
        # handle chr prefix
        if "chr" in chromosomes["chromosome"][0]:
            seg_open = pd.read_csv(input.filled, sep = "\t")
            chrom = list(seg_open['chrom'])
            # avoid cases of chrchr1 if the prefix already there
            for i in range(len(chrom)):
                if 'chr' not in str(chrom[i]):
                    chrom[i]='chr'+str(chrom[i])
            seg_open.loc[:, 'chrom']=chrom
            seg_open.to_csv(output.projection, sep="\t", index=False, na_rep='NA')
        else:
            # remove chr prefix
            seg_open = pd.read_csv(input.filled, sep = "\t")
            seg_open["chrom"] = seg_open["chrom"].astype(str).str.replace('chr', '')
            seg_open.to_csv(output.projection, sep="\t", index=False, na_rep='NA')


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _controlfreec_output_projection:
    input:
        projection = str(rules._controlfreec_normalize_projection.output.projection)
    output:
        projection = CFG["dirs"]["outputs"] + "seg/{seq_type}--projection/{masked}/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    threads: 1
    run:
        op.relative_symlink(input.projection, output.projection, in_module = True)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _controlfreec_output:
    input:
        plot = str(rules._controlfreec_plot.output.plot),
        log2plot = str(rules._controlfreec_plot.output.log2plot),
        CNV = str(rules._controlfreec_calc_sig.output.txt),
        bed = str(rules._controlfreec_freec2bed.output.bed),
        BAFgraph = str(rules._controlfreec_plot.output.bafplot),
        circos = str(rules._controlfreec_freec2circos.output.circos),
        igv = str(rules._controlfreec_cnv2igv.output.seg)
    output:
        plot = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}.ratio.png",
        log2plot = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}.ratio.log2.png",
        CNV = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}.CNVs.txt",
        bed = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}.CNVs.bed",
        BAFgraph = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}.BAF.png",
        circos = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}.circos.bed",
        igv = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{masked}/{tumour_id}--{normal_id}--{pair_status}.CNVs.seg"
    run:
        op.relative_symlink(input.plot, output.plot, in_module = True)
        op.relative_symlink(input.log2plot, output.log2plot, in_module = True)
        op.relative_symlink(input.CNV, output.CNV, in_module = True)
        op.relative_symlink(input.bed, output.bed, in_module = True)
        op.relative_symlink(input.BAFgraph, output.BAFgraph, in_module = True)
        op.relative_symlink(input.circos, output.circos, in_module = True)
        op.relative_symlink(input.igv, output.igv, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _controlfreec_all:
    input:
        expand(
            expand(
            [
                str(rules._controlfreec_output.output.plot),
                str(rules._controlfreec_output.output.log2plot),
                str(rules._controlfreec_output.output.CNV),
                str(rules._controlfreec_output.output.bed),
                str(rules._controlfreec_output.output.BAFgraph),
                str(rules._controlfreec_output.output.circos),
                str(rules._controlfreec_output.output.igv),
                str(rules._controlfreec_run.output.done)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            pair_status=CFG["runs"]["pair_status"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            allow_missing = True
            ),
         masked=["masked", "unmasked"]
        ),
        expand(
            expand(
                [
                    str(rules._controlfreec_output_projection.output.projection)
                ],
                zip,
                seq_type=CFG["runs"]["tumour_seq_type"],
                pair_status=CFG["runs"]["pair_status"],
                tumour_id=CFG["runs"]["tumour_sample_id"],
                normal_id=CFG["runs"]["normal_sample_id"],
                allow_missing = True
            ),
        tool="controlfreec",
        masked=["masked", "unmasked"],
        projection=CFG["requested_projections"]
        )




##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
