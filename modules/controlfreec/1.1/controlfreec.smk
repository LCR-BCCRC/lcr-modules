#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Jasper
# Module Author:    Jasper
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`

# `CFG` is a shortcut to `config["lcr-modules"]["controlfreec"]`
CFG = op.setup_module(
    name = "controlfreec",
    version = "1.1",
    subdirectories = ["inputs", "run", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _controlfreec_input_bam,
    _controlfreec_config,
    _controlfreec_output,
    _controlfreec_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _controlfreec_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
        

# generate references
rule _get_map_refs:
    output:
        tar = temp("references/{genome_build}/freec/out100m2_{genome_build}.tar.gz"),
        gem = "references/{genome_build}/freec/out100m2_{genome_build}.gem"
    params:
        provider = "ensembl",
        url = lambda w: {"grch37": "http://xfer.curie.fr/get/7hZIk1C63h0/hg19_len100bp.tar.gz",
                        "hg19": "http://xfer.curie.fr/get/7hZIk1C63h0/hg19_len100bp.tar.gz",
                        "grch38": "http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip",
                        "hg38": "http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip"}[w.genome_build],
        command1 = lambda w: {"grch37": "tar -xvf ",
                            "hg19": "tar -xvf ",
                            "grch38": "unzip ",
                            "hg38": "unzip "}[w.genome_build],
        command2 = lambda w: {"grch37": " --wildcards --no-anchored 'out100m2*gem' && mv out100m2_hg19.gem ",
                            "hg19": " --wildcards --no-anchored 'out100m2*gem' && mv out100m2_hg19.gem ",
                            "grch38": " > ",
                            "hg38": " > "}[w.genome_build]
    shell:
        "wget -O {output.tar} {params.url} "
        "&& {params.command1} {output.tar} {params.command2} {output.gem}"

# no chr for grch37 and grch38
# chr for hg19 and hg38
chromRegions = ("1", "2", "3", "4", "5", "6", "7", "8",
        "9", "10", "11", "12", "13", "14", "15", "16",
        "17", "18", "19", "20", "21", "22", "X", "Y")



rule _generate_chrLen:
    input:
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        chrLen = "references/{genome_build}/freec/{genome_build}.len"
    shell:
        op.as_one_line("""
            grep -P '^chr[0-9,X,Y]+\t|^[0-9,X,Y]' {input.fai} | awk '{{print $1"\t"$2}}' > {output.chrLen}
        """)


rule _generate_chrFasta:
    input:
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        fasta = "references/{genome_build}/freec/chr/{chromosome}.fa"
    conda:
        CFG["conda_envs"]["controlfreec"]
    params:
        chr_prefix = lambda w: {"grch37": " ",
                        "hg19": "chr",
                        "grch38": " ",
                        "hg38": "chr"}[w.genome_build]
    shell:
        "samtools faidx {input.fasta} {params.chr_prefix}{wildcards.chromosome} > {output.fasta} "


rule _check_chrFiles:
    input:
        fasta = expand("references/{{genome_build}}/freec/chr/{chromosome}.fa",
                    chromosome = chromRegions)
    output:
        touch("references/{genome_build}/freec/chr/.all_done")


# set-up controlfreec
rule _controlfreec_config:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        reference = "references/{genome_build}/freec/out100m2_{genome_build}.gem",
        chrLen = "references/{genome_build}/freec/{genome_build}.len",
        done = "references/{genome_build}/freec/chr/.all_done"
    output:
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/config_WGS.txt"
    conda:
        CFG["conda_envs"]["controlfreec"]
    params:
        config = CFG["options"]["configFile"],
        outdir = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/",
        window = CFG["options"]["window"],
        ploidy = CFG["options"]["ploidy"],
        breakPointValue = CFG["options"]["breakPointThreshold"],
        coefVar = CFG["options"]["coefficientOfVariation"],
        numCon = CFG["options"]["contamination"],
        booCon = CFG["options"]["contaminationAdjustment"],
        degree = CFG["options"]["degree"],
        forceGC = CFG["options"]["forceGCcontentNormalization"],
        chrFiles = "references/{genome_build}/freec/chr/",
        minCNAlength = CFG["options"]["minCNAlength"],
        minimumSubclonePresence = CFG["options"]["minimalSubclonePresence"],
        noisyData = CFG["options"]["noisyData"],
        step = CFG["options"]["step"],
        telocentromeric = CFG["options"]["telocentromeric"],
        threads = CFG["options"]["maxThreads"]
    shell:
        "samtoolspath=$(which samtools ) ; "
        "samtoolsPathName=$(echo $samtoolspath) ; "
        "sambambapath=$(which sambamba ) ; "
        "sambambaPathName=$(echo $sambambapath) ; "
        "bedtoolspath=$(which bedtools ) ; "
        "bedtoolsPathName=$(echo $bedtoolspath) ; "
        "sed \"s|BAMFILE|{input.bam}|g\" {params.config} | "
        "sed \"s|OUTDIR|{params.outdir}|g\" | "
        "sed \"s|windowSize|{params.window}|g\" | "
        "sed \"s|ploidyInput|{params.ploidy}|g\" | "
        "sed \"s|chrLenFileInput|{input.chrLen}|g\" | "
        "sed \"s|chrFilesPath|{params.chrFiles}|g\" | "
        "sed \"s|sambambaPath|$sambambaPathName|g\" | "
        "sed \"s|bedtoolsPath|$bedtoolsPathName|g\" | "
        "sed \"s|samtoolsPath|$samtoolsPathName|g\" | "
        "sed \"s|breakPointValue|{params.breakPointValue}|g\" | "
        "sed \"s|coefVar|{params.coefVar}|g\" | "
        "sed \"s|numCon|{params.numCon}|g\" | "
        "sed \"s|booCon|{params.booCon}|g\" | "
        "sed \"s|forceGCvalue|{params.forceGC}|g\" | "
        "sed \"s|numDegree|{params.degree}|g\" | "
        "sed \"s|minCNAvalue|{params.minCNAlength}|g\" | "
        "sed \"s|minimumSubclonePresenceValue|{params.minimumSubclonePresence}|g\" | "
        "sed \"s|booNoise|{params.noisyData}|g\" | "
        "sed \"s|stepValue|{params.step}|g\" | "
        "sed \"s|teloValue|{params.telocentromeric}|g\" | "
        "sed \"s|numThreads|{params.threads}|g\" | "
        "sed \"s|referenceFile|{input.reference}|g\" > {output}"


rule _controlfreec_run:
    input:
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/config_WGS.txt"
    output:
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_info.txt",
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_ratio.txt",
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_CNVs"
    conda: CFG["conda_envs"]["controlfreec"]
    threads: CFG["threads"]["controlfreec_run"]
    resources: mem_mb = CFG["mem_mb"]["controlfreec_run"]
    log:
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/run.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/run.stderr.log"
    shell:
        "freec -conf {input} > {log.stdout} 2> {log.stderr} "


rule _controlfreec_calc_sig:
    input:
        CNVs = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_CNVs",
        ratios = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_ratio.txt",
    output:
        txt = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_CNVs.p.value.txt"
    params:
        calc_sig = CFG["software"]["FREEC_sig"]
    threads: CFG["threads"]["calc_sig"]
    resources: mem_mb = CFG["mem_mb"]["calc_sig"]
    conda: CFG["conda_envs"]["controlfreec"]
    log:         
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/calc_sig.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/calc_sig.stderr.log"
    shell:
        "cat {params.calc_sig} | R --slave --args {input.CNVs} {input.ratios} > {log.stdout} 2> {log.stderr}"


rule _controlfreec_plot:
    input:
        ratios = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_ratio.txt",
        # BAF = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_BAF.txt",
        info = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_info.txt"
    output:
        plot = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_ratio.txt.png",
        log2plot = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_ratio.txt.log2.png"
    params:
        plot = CFG["software"]["FREEC_graph"]
    threads: CFG["threads"]["plot"]
    resources: mem_mb = CFG["mem_mb"]["plot"]
    conda: CFG["conda_envs"]["controlfreec"]
    log: 
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/plot.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/plot.stderr.log"
    shell:
        "cat {params.plot} | R --slave --args `grep \"Output_Ploidy\" {input.info} | cut -f 2` {input.ratios} > {log.stdout} 2> {log.stderr} "


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _controlfreec_output:
    input:
        plot = str(rules._controlfreec_plot.output.plot),
        log2plot = str(rules._controlfreec_plot.output.log2plot),
        CNV = str(rules._controlfreec_calc_sig.output.txt)
    output:
        plot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/plots/{sample_id}.bam_ratio.txt.png",
        log2plot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/log2plots/{sample_id}.bam_ratio.txt.log2.png",
        CNV = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/CNV/{sample_id}.bam_CNVs.p.value.txt"
    run:
        op.relative_symlink(input.plot, output.plot)
        op.relative_symlink(input.log2plot, output.log2plot)
        op.relative_symlink(input.CNV, output.CNV)


# Generates the target sentinels for each run, which generate the symlinks
rule _controlfreec_all:
    input:
        expand(
            [
                str(rules._controlfreec_output.output.plot),
                str(rules._controlfreec_output.output.log2plot),
                str(rules._controlfreec_output.output.CNV)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            sample_id=CFG["runs"]["tumour_sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
