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
    subdirectories = ["inputs", "mpileup", "run", "calc_sig", "plot", "freec2bed", "freec2circos", "cnv2igv", "convert_coordinates", "fill_regions", "normalize", "outputs"]
)



# Define rules to be run locally when using a compute cluster
localrules:
    _controlfreec_input_bam,
    _controlfreec_config,
    _controlfreec_plot,
    _controlfreec_output,
    _controlfreec_all


##### RULES #####

#### Rules for mappability reference
# to generate and use hard-masked mappability (i.e. recommended for FFPE genomes) if CFG["options"]["hard_masked"] == True
# to use the default genome's mappability file (downloaded from their website), set it CFG["options"]["hard_masked"] == False
if CFG["options"]["hard_masked"] == True:
    CFG["runs"]["masked"] = "_masked"
else:
    CFG["runs"]["masked"] = ""

wildcard_constraints:
    masked = ".{0}|_masked",
    genome_build = ".+(?<!masked)"

#### generate references ####
# mappability tracks for hg19 and hg38 are available from the source
# even though the hg38 is a .zip, this saves it to the same tar path but unpacks it differently
if CFG["options"]["hard_masked"] == False:
    rule _controlfreec_get_map_refs:
        output:
            tar = temp(CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/out100m2_{genome_build}.tar.gz"),
            gem = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/out100m2_{genome_build}.gem"
        params:
            outdir = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/"
        shell:
            """
            if echo "{wildcards.genome_build}" | grep -q "38"; then
                wget -O {output.tar} http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip && \
                unzip {output.tar} -d {params.outdir} && \
                mv {params.outdir}out100m2_hg38.gem {output.gem}
            else
                wget -O {output.tar} http://xfer.curie.fr/get/7hZIk1C63h0/hg19_len100bp.tar.gz && \
                tar -xvf {output.tar} -C {params.outdir} --wildcards --no-anchored 'out100m2*gem' && \
                mv {params.outdir}out100m2_hg19.gem {output.gem}
                fi
            """

# mappability tracks for hard-masked genomes need to be generated using GEM
rule _download_GEM:
    output:
        touch(CFG["dirs"]["inputs"] + "references/GEM/.done")
    params:
        dirOut = CFG["dirs"]["inputs"] + "references/GEM/"
    conda:
        CFG["conda_envs"]["wget"]
    resources: **CFG["resources"]["gem"]
    shell:
        "wget https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download -O {params.dirOut}/GEM-lib.tbz2 && bzip2 -dc {params.dirOut}/GEM-lib.tbz2 | tar -xvf - -C {params.dirOut}/"

# grch37 and grch38 from ensembl have additional information in header - need to remove
if CFG["options"]["hard_masked"] == True:
    rule _set_up_grch_genomes:
        input:
            reference = reference_files("genomes/{genome_build}{masked}/genome_fasta/genome.fa")
        output:
            reference = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/genome_header.fa"
        resources: **CFG["resources"]["gem"]
        shell:
            "cat {input.reference} | perl -ne 's/(^\>\S+).+/$1/;print;' > {output.reference} "

def get_genome_fasta(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    if  "grch" in str({wildcards.genome_build}):
        return  CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/genome_header.fa"
    else:
        return reference_files("genomes/{genome_build}{masked}/genome_fasta/genome.fa")

if CFG["options"]["hard_masked"] == True:
    rule _generate_gem_index:
        input:
            software = CFG["dirs"]["inputs"] + "references/GEM/.done",
            reference = get_genome_fasta
        output:
            index = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/{genome_build}.hardmask.all_index.gem"
        params:
            gemDir = CFG["dirs"]["inputs"] + "references/GEM/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin",
            idxpref = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/{genome_build}.hardmask.all_index"
        threads: CFG["threads"]["gem"]
        resources: **CFG["resources"]["gem"]
        log: CFG["logs"]["inputs"] + "gem/{genome_build}{masked}/gem_index.stderr.log"
        shell:
            "PATH=$PATH:{params.gemDir}; {params.gemDir}/gem-indexer -T {threads} -c dna -i {input.reference} -o {params.idxpref} > {log} 2>&1 "

if CFG["options"]["hard_masked"] == True:
    rule _generate_mappability:
        input:
            software = CFG["dirs"]["inputs"] + "references/GEM/.done",
            index = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/{genome_build}.hardmask.all_index.gem"
        output:
            mappability = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/{genome_build}.hardmask.all.gem.mappability"
        params:
            gemDir = CFG["dirs"]["inputs"] + "references/GEM/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin",
            pref =  CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/{genome_build}.hardmask.all.gem",
            kmer = CFG["options"]["kmer"],
            mismatch = CFG["options"]["mismatch"],
            maxEditDistance = CFG["options"]["maxEditDistance"],
            maxBigIndel = CFG["options"]["maxBigIndel"],
            strata = CFG["options"]["strata"]
        threads: CFG["threads"]["gem"]
        resources: **CFG["resources"]["gem"]
        log: CFG["logs"]["inputs"] + "gem/{genome_build}{masked}/gem_map.stderr.log"
        shell:
            "PATH=$PATH:{params.gemDir}; {params.gemDir}/gem-mappability -T {threads} -I {input.index} -l {params.kmer} -m {params.mismatch} -t disable --mismatch-alphabet ACGNT -e {params.maxEditDistance} --max-big-indel-length {params.maxBigIndel} -s {params.strata} -o {params.pref} > {log} 2>&1 "

if CFG["options"]["hard_masked"] == True:
    rule _symlink_map:
        input:
            mappability = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/{genome_build}.hardmask.all.gem.mappability"
        output:
            mappability = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/out100m2_{genome_build}.gem"
        resources: **CFG["resources"]["gem"]
        shell:
            "ln -srf {input.mappability} {output.mappability} "

#### Rule for setting chromosome names (chr-prefix or not)
# no chr for grch37 and grch38
# chr for hg19 and hg38
# chromosomes used (i.e. chr1-22,X,Y)


def _controlfreec_get_chr_fastas(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    chrs = []
    for i in range(1, 23):
        chrs.append(str(i))
    chrs.extend(["X", "Y"])
    if str(wildcards.genome_build).startswith("hg"):
        chrs = ["chr" + x for x in chrs]
    fastas = expand(
        CFG["dirs"]["inputs"] +  "references/{{genome_build}}/freec/chr/{chromosome}.fa",
        chromosome = chrs
    )
    return(fastas)



#generates file with chromomsome lengths from genome.fa.fai
rule _controlfreec_generate_chrLen:
    input:
        fai = reference_files("genomes/{genome_build}{masked}/genome_fasta/genome.fa.fai"),
        main = reference_files("genomes/{genome_build}{masked}/genome_fasta/main_chromosomes_withY.txt")
    output:
        chrLen = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/{genome_build}.len"
    resources: **CFG["resources"]["gem"]
    shell:
        op.as_one_line("""
            grep -P '^chr[0-9,X,Y]+\t|^[0-9,X,Y]' {input.fai} | awk '{{print $1"\t"$2}}' > {output.chrLen}
        """)

#generates a per-chromosome fasta file from genome.fa
rule _controlfreec_generate_chrFasta:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        fasta = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/{chromosome}.fa"
    conda:
        CFG["conda_envs"]["controlfreec"]
    resources: **CFG["resources"]["gem"]
    shell:
        "samtools faidx {input.fasta} {wildcards.chromosome} > {output.fasta} "

#checks that all chroms are accounted for
rule _controlfreec_check_chrFiles:
    input:
        _controlfreec_get_chr_fastas
    wildcard_constraints:
        genome_build = "|".join(CFG["runs"]["tumour_genome_build"])
    output:
        touch(CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/.all_done")


rule _controlfreec_dbsnp_to_bed:
    input:
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        bed = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/dbsnp.common_all-151.bed"
    resources: **CFG["resources"]["gem"]
    shell:
        op.as_one_line(""" gunzip -c {input.vcf} | awk {{'printf ("%s\\t%s\\t%s\\n", $1,$2-1,$2)'}} | zgrep -v -h "^#" > {output.bed} """)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _controlfreec_input_bam:
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
    output: # creates a temporary file for mpileup
        pileup = temp(CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.minipileup.pileup.gz")
    conda:
        CFG["conda_envs"]["controlfreec"]
    threads: 1
    resources:
        **CFG["resources"]["mpileup"]
    group: "controlfreec"
    log:
        stderr = CFG["logs"]["inputs"] + "mpileup/{seq_type}--{genome_build}/{sample_id}.{chrom}.mpileup.stderr.log",
    shell:
        "samtools mpileup -l {input.bed} -r {wildcards.chrom} -Q 20 -f {input.fastaFile} {input.bam} | gzip -c > {output.pileup} 2> {log.stderr}"

#### set-up mpileups for BAF calling ####
def _controlfreec_get_chr_mpileups(wildcards):
    chrs = []
    for i in range(1, 23):
        chrs.append(str(i))
    chrs.extend(["X", "Y"])
    if str(wildcards.genome_build).startswith("hg"):
        chrs = ["chr" + x for x in chrs]

    CFG = config["lcr-modules"]["controlfreec"]

    mpileups = expand(
        CFG["dirs"]["mpileup"] + "{{seq_type}}--{{genome_build}}/{{sample_id}}.{chrom}.minipileup.pileup.gz",
        chrom = chrs
    )
    return(mpileups)

rule _controlfreec_concatenate_pileups:
    input:
        mpileup = _controlfreec_get_chr_mpileups
    output:
        mpileup = temp(CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.bam_minipileup.pileup.gz")
    resources:
        **CFG["resources"]["cat"]
    wildcard_constraints:
        genome_build = "|".join(CFG["runs"]["tumour_genome_build"])
    group: "controlfreec"
    shell:
        "cat {input.mpileup} > {output.mpileup} "


#### Run control-FREEC ####

# set-up controlfreec
# these have to be each their own function, bc `params:` can't handle dicts or use unpack()
def _controlfreec_get_optional_step(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    if CFG["options"]["step"] != "":
        cfg_step = "step = " + str(CFG["options"]["step"]) # will replace the entire line, so it's no longer commented out
    else:
        cfg_step = "#step = " + str(CFG["options"]["step"]) # will stay commented out, won't be used, just needed for the code below
    return(cfg_step)

def _controlfreec_get_optional_window(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    if CFG["options"]["window"] != "":
        cfg_window = "window = " + str(CFG["options"]["window"]) # will replace the entire line, so it's no longer commented out
    else:
        cfg_window = "#window = " + str(CFG["options"]["window"]) # will stay commented out, won't be used, just needed for the code below
    return(cfg_window)

rule _controlfreec_config_contamAdjTrue:
    input:
        tumour_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.bam_minipileup.pileup.gz",
        normal_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.bam_minipileup.pileup.gz",
        reference = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/out100m2_{genome_build}.gem",
        chrLen = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/{genome_build}.len",
        done = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/.all_done"
    output:
        config = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjTrue/config.txt"
    conda:
        CFG["conda_envs"]["controlfreec"]
    params:
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
        dbSNP = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz"),
        outdir = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjTrue/",
        chrFiles = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/"
    shell:
        "samtoolspath=$(which samtools ) ; "
        "samtoolsPathName=$(echo $samtoolspath) ; "
        "sambambapath=$(which sambamba ) ; "
        "sambambaPathName=$(echo $sambambapath) ; "
        "bedtoolspath=$(which bedtools ) ; "
        "bedtoolsPathName=$(echo $bedtoolspath) ; "
        "sed \"s|BAMFILE|{input.tumour_bam}|g\" {output.config} | "
        "sed \"s|CONTROLFILE|{input.normal_bam}|g\" | "
        "sed \"s|OUTDIR|{params.outdir}|g\" | "
        "sed \"s|DBsnpFile|{params.dbSNP}|g\" | "
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
        "sed \"s|referenceFile|{input.reference}|g\" > {output.config}"


rule _controlfreec_config_contamAdjFalse:
    input:
        tumour_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.bam_minipileup.pileup.gz",
        normal_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.bam_minipileup.pileup.gz",
        reference = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/out100m2_{genome_build}.gem",
        chrLen = CFG["dirs"]["inputs"] + "references/{genome_build}{masked}/freec/{genome_build}.len",
        done = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/.all_done"
    output:
        config = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjFalse/config.txt"
    conda:
        CFG["conda_envs"]["controlfreec"]
    params:
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
        dbSNP = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz"),
        outdir = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjFalse/",
        chrFiles = CFG["dirs"]["inputs"] + "references/{genome_build}/freec/chr/"
    shell:
        "samtoolspath=$(which samtools ) ; "
        "samtoolsPathName=$(echo $samtoolspath) ; "
        "sambambapath=$(which sambamba ) ; "
        "sambambaPathName=$(echo $sambambapath) ; "
        "bedtoolspath=$(which bedtools ) ; "
        "bedtoolsPathName=$(echo $bedtoolspath) ; "
        "sed \"s|BAMFILE|{input.tumour_bam}|g\" {output.config} | "
        "sed \"s|CONTROLFILE|{input.normal_bam}|g\" | "
        "sed \"s|OUTDIR|{params.outdir}|g\" | "
        "sed \"s|DBsnpFile|{params.dbSNP}|g\" | "
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
        "sed \"s|referenceFile|{input.reference}|g\" > {output.config}"

rule _controlfreec_run:
    input:
        config = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--{contamAdj}/config.txt",
        tumour_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.bam_minipileup.pileup.gz",
        normal_bam = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.bam_minipileup.pileup.gz",
    output:
        done = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--{contamAdj}/{tumour_id}.done"
    conda: CFG["conda_envs"]["controlfreec"]
    threads: CFG["threads"]["controlfreec_run"]
    resources: **CFG["resources"]["controlfreec_run"]
    log:
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--{contamAdj}/run.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--{contamAdj}/run.stderr.log"
    run:
        shell("freec -conf {input.config} > {log.stdout} 2> {log.stderr} || touch {output.done}")

def _get_run_result(wildcards):
    CFG = config["lcr-modules"]["controlfreec"]
    if os.path.exists(CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjTrue/{tumour_id}.bam_minipileup.pileup.gz_CNVs")
        info = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjTrue/{tumour_id}.bam_minipileup.pileup.gz_info.txt",
        ratios = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjTrue/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt",
        CNV = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjTrue/{tumour_id}.bam_minipileup.pileup.gz_CNVs",
        BAF = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjTrue/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt"
    else:
        info = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjFalse/{tumour_id}.bam_minipileup.pileup.gz_info.txt",
        ratios = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjFalse/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt",
        CNV = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjFalse/{tumour_id}.bam_minipileup.pileup.gz_CNVs",
        BAF = CFG["dirs"]["run"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}--contamAdjFalse/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt"
    return{
        "info": info,
        "ratios": ratios,
        "CNV": CNV,
        "BAF": BAF
    }

rule _controlfreec_calc_sig:
    input:
        unpack(_get_run_result),
        str(rules._controlfreec_run.output.done)
    output:
        txt = CFG["dirs"]["calc_sig"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_CNVs.p.value.txt"
    params:
        calc_sig = CFG["software"]["FREEC_sig"]
    threads: CFG["threads"]["calc_sig"]
    resources: **CFG["resources"]["calc_sig"]
    conda: CFG["conda_envs"]["controlfreec"]
    log:
        stdout = CFG["logs"]["calc_sig"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/calc_sig.stdout.log",
        stderr = CFG["logs"]["calc_sig"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/calc_sig.stderr.log"
    shell:
        "Rscript --vanilla {params.calc_sig} {input.CNV} {input.ratios} > {log.stdout} 2> {log.stderr}"


rule _controlfreec_plot:
    input:
        unpack(_get_run_result),
        str(rules._controlfreec_run.output.done)
    output:
        plot = CFG["dirs"]["plot"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt.png",
        log2plot = CFG["dirs"]["plot"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_ratio.txt.log2.png",
        bafplot = CFG["dirs"]["plot"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam_minipileup.pileup.gz_BAF.txt.png"
    params:
        plot = CFG["software"]["FREEC_graph"]
    resources:
        bam = 1
    conda: CFG["conda_envs"]["controlfreec"]
    log:
        stdout = CFG["logs"]["plot"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/plot.stdout.log",
        stderr = CFG["logs"]["plot"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/plot.stderr.log"
    shell:
        "Rscript --vanilla {params.plot} `grep \"Output_Ploidy\" {input.info} | cut -f 2` {input.ratios} {input.BAF} > {log.stdout} 2> {log.stderr} "


rule _controlfreec_freec2bed:
    input:
        unpack(_get_run_result),
        str(rules._controlfreec_run.output.done)
    output:
        bed = CFG["dirs"]["freec2bed"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bed"
    params:
        freec2bed = CFG["software"]["freec2bed"]
    threads: CFG["threads"]["freec2bed"]
    resources: **CFG["resources"]["freec2bed"]
    conda: CFG["conda_envs"]["controlfreec"]
    log:
        stderr = CFG["logs"]["freec2bed"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/freec2bed.stderr.log"
    shell:
        "ploidy=$(grep Output_Ploidy {input.info} | cut -f 2); "
        "perl {params.freec2bed} -f {input.ratios} -p $ploidy > {output.bed} 2> {log.stderr}"


rule _controlfreec_freec2circos:
    input:
        unpack(_get_run_result),
        str(rules._controlfreec_run.output.done)
    output:
        circos = CFG["dirs"]["freec2circos"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.circos.bed"
    params:
        freec2circos = CFG["software"]["freec2circos"]
    threads: CFG["threads"]["freec2circos"]
    resources: **CFG["resources"]["freec2circos"]
    conda: CFG["conda_envs"]["controlfreec"]
    log:
        stderr = CFG["logs"]["freec2circos"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/freec2circos.stderr.log"
    shell:
        "ploidy=$(grep Output_Ploidy {input.info} | cut -f 2); "
        "perl {params.freec2circos} -f {input.ratios} -p $ploidy > {output.circos} 2> {log.stderr}"


rule _controlfreec_cnv2igv:
    input:
        cnv = str(rules._controlfreec_calc_sig.output.txt),
        cnv2igv = ancient(CFG["inputs"]["cnv2igv"])
    output:
        seg = CFG["dirs"]["cnv2igv"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.CNVs.seg"
    params:
        opts = CFG["options"]["preserve"]
    conda:
        CFG["conda_envs"]["cnv2igv"]
    threads: 1
    log:
        stderr = CFG["logs"]["cnv2igv"] + "{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}/cnv2igv.stderr.log"
    group: "controlfreec_post_process"
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
        controlfreec_lifted = CFG["dirs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.stderr.log"
    threads: 1
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        liftover_minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    group: "controlfreec_post_process"
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
        hg38_projection = str(rules._controlfreec_cnv2igv.output.seg).replace("{genome_build}", this_genome_build[0]).replace("{masked}", this_masked[0])
        grch37_projection = str(rules._controlfreec_convert_coordinates.output.controlfreec_lifted).replace("{genome_build}", this_genome_build[0]).replace("{masked}", this_masked[0])
        grch37_projection = grch37_projection.replace("{chain}", "hg38ToHg19")
    else:
        grch37_projection = str(rules._controlfreec_cnv2igv.output.seg).replace("{genome_build}", this_genome_build[0]).replace("{masked}", this_masked[0])
        hg38_projection = str(rules._controlfreec_convert_coordinates.output.controlfreec_lifted).replace("{genome_build}", this_genome_build[0]).replace("{masked}", this_masked[0])
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
        grch37_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"),
        hg38_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg")
    log:
        stderr = CFG["logs"]["fill_regions"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}_fill_segments.stderr.log"
    threads: 1
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"]
    conda:
        CFG["conda_envs"]["bedtools"]
    group: "controlfreec_post_process"
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
        filled = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg",
        chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
    output:
        projection = CFG["dirs"]["normalize"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    resources:
        **CFG["resources"]["post_controlfreec"]
    threads: 1
    group: "controlfreec_post_process"
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
        projection = CFG["dirs"]["outputs"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    threads: 1
    group: "controlfreec_post_process"
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
        ratio = str(rules._controlfreec_run.output.ratio),
        circos = str(rules._controlfreec_freec2circos.output.circos),
        igv = str(rules._controlfreec_cnv2igv.output.seg)
    output:
        plot = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}.ratio.png",
        log2plot = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}.ratio.log2.png",
        CNV = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}.CNVs.txt",
        bed = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}.CNVs.bed",
        BAFgraph = CFG["dirs"]["outputs"] + "png/{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}.BAF.png",
        ratio = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}.ratio.txt",
        circos = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}.circos.bed",
        igv = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}{masked}/{tumour_id}--{normal_id}--{pair_status}.seg"
    group: "controlfreec_post_process"
    run:
        op.relative_symlink(input.plot, output.plot, in_module = True)
        op.relative_symlink(input.log2plot, output.log2plot, in_module = True)
        op.relative_symlink(input.CNV, output.CNV, in_module = True)
        op.relative_symlink(input.bed, output.bed, in_module = True)
        op.relative_symlink(input.BAFgraph, output.BAFgraph, in_module = True)
        op.relative_symlink(input.ratio, output.ratio, in_module = True)
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
                str(rules._controlfreec_output.output.ratio),
                str(rules._controlfreec_output.output.circos),
                str(rules._controlfreec_output.output.igv)
                str(rules._controlfreec_output_projection.output.projection)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            pair_status=CFG["runs"]["pair_status"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            masked=CFG["runs"]["masked"],
            allow_missing=True
            ),
        tool="controlfreec",
        projection=CFG["requested_projections"]
        )



##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
