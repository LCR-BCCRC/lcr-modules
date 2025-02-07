"""A pipeline to use a panel of normals to collect variant calls that are
common SNPs and artifacts, and then use this panel to filter out these
variants from the tumor samples. This is done to reduce the number of false
positives in the tumor samples.
"""
import os
import glob

# setup some helpful globals
SAMPLESHEET = config["lcr-modules"]["_shared"]["samples"]
BAM_DIR = config["lcr-modules"]["make_pon"]["bam_dir"]
WORKDIR = config["lcr-modules"]["make_pon"]["workdir"]
# find all bams in the bam directory
BAMS = glob.glob(os.path.join(BAM_DIR, "*.bam"))
LOGSDIR = os.path.join(WORKDIR, "logs")
###### setup functions 

def divide_targe_bed(in_bed: str, chromosomes: list) -> dict:
    """Divide the target bed into chromosomes.

    Write out dividided targets to individual beds to be used by
    mutect2. Return a list of the outpaths to passed to mutect2.
    """
    indir = os.path.join(WORKDIR,"00-inputs/targets")
    # if it doesnt exist, make it
    if not os.path.exists(indir):
        os.makedirs(indir)
    # divide beds
    bed = pd.read_csv(in_bed, sep="\t", header=None)
    # get number of columns
    if bed.shape[1] == 4:
        bed.columns = ["chrom", "start", "end", "name"]
    elif bed.shape[1] == 3:
        bed.columns = ["chrom", "start", "end"]
    else:
        raise ValueError("Bed file must have 3 or 4 columns")
    bed["chrom"] = bed["chrom"].astype(str)
    chroms = chromosomes
    beds = []
    for chrom in chroms:
        chrmbed = bed[bed["chrom"] == chrom]
        # if no targets for chrom, skip
        if chrmbed.empty:
            continue
        beds.append(chrmbed)

    outpaths = {}
    for i, b in enumerate(beds):
        outpath = os.path.join(indir, f"{chroms[i]}.bed")
        b.to_csv(outpath, sep="\t", header=False, index=False)
        outpaths[chroms[i]] = outpath

    return outpaths

def get_genome_chromosomes(genome_build)-> list:
    """Return a list of chromosomes for the genome build.
    """
    if genome_build.lower() in ["hg38","hg19-reddy", "grch38"]:
        return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    else:
        return [str(i) for i in range(1, 23)] + ["X", "Y"]

# if a bed file of targets provided than divide those by chromosome, if not just pass list of chromosomes
if config["lcr-modules"]["make_pon"]["target_bed"]:
    TARGET_INTERVALS = divide_targe_bed(config["lcr-modules"]["make_pon"]["target_bed"], 
                                get_genome_chromosomes(config["lcr-modules"]["make_pon"]["genome_build"]))
else:
    TARGET_INTERVALS = get_genome_chromosomes(config["lcr-modules"]["make_pon"]["genome_build"])

####### RULES
localrules:
    symlink_inputs,
    generate_samples_map


def find_norm_bam(wildcards):
    # find normal bam in list of BAMS
    bam = [b for b in BAMS if f"{wildcards.normal_id}" in os.path.basename(b)]
    if len(bam) == 0:
        raise FileNotFoundError(f"No bam found for {wildcards.normal_id}")
    if len(bam) > 1:
        raise FileNotFoundError(f"Multiple bams found for {wildcards.normal_id}\n {bam}")
    return bam

rule symlink_inputs:
    input:
        inbam = find_norm_bam
    output:
        bam = os.path.join(WORKDIR, "00-inputs/bams", "{normal_id}.bam")
    shell:
        "ln -s {input.inbam} {output.bam}" 

rule make_pon_mutect2_germline:
    input:
        bam = rules.symlink_inputs.output.bam,
    output:
        vcf = temp(os.path.join(WORKDIR,"01-Normals/" ,"{seq_type}/{capture_space}--{genome_build}/{pon}/{normal_id}/{normal_id}.{interval}.vcf.gz")),
        tbi = temp(os.path.join(WORKDIR,"01-Normals/" ,"{seq_type}/{capture_space}--{genome_build}/{pon}/{normal_id}/{normal_id}.{interval}.vcf.gz.tbi")),
        stats = temp(os.path.join(WORKDIR,"01-Normals/" ,"{seq_type}/{capture_space}--{genome_build}/{pon}/{normal_id}/{normal_id}.{interval}.vcf.gz.stats")),
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8), 
        fasta = config["lcr-modules"]["make_pon"]["ref_fasta"],
        gnomad = config["lcr-modules"]["make_pon"]["gnomad"],
        # interval = TARGET_INTERVALS["{interval}"]
        interval = lambda wildcards: TARGET_INTERVALS[wildcards.interval]
    log: os.path.join(LOGSDIR, "norm_mut2","{seq_type}/{capture_space}--{genome_build}/{pon}/{normal_id}/{interval}.log") 
    conda: config["lcr-modules"]["make_pon"]["gatk_env"]
    resources: **config["lcr-modules"]["make_pon"]["resources"]["mutect"]
    threads: 1
    shell:
        """
            gatk Mutect2 --java-options "-Xmx{params.mem_mb}m" --genotype-germline-sites true --genotype-pon-sites true --interval-padding 150 --max-mnp-distance 0 \
            --germline-resource {params.gnomad} -R {params.fasta} -L {params.interval} -I {input.bam} -O {output.vcf} > {log} 2>&1
        """

def mutect2_normal_get_chr_vcf(wildcards):
    return expand(expand(str(rules.make_pon_mutect2_germline.output.vcf), zip,
    seq_type = wildcards.seq_type,
    genome_build = wildcards.genome_build,
    capture_space = wildcards.capture_space,
    normal_id = wildcards.normal_id,
    pon = wildcards.pon, allow_missing = True),
   interval = list(TARGET_INTERVALS.keys()), )

rule make_pon_mutect2_concat_norm_vcf:
    input: 
        vcf = mutect2_normal_get_chr_vcf,
    output: 
        vcf = os.path.join(WORKDIR,"02-Normals/" , "{seq_type}/{capture_space}--{genome_build}/{normal_id}-{pon}/{normal_id}.vcf.gz"),
        tbi = os.path.join(WORKDIR,"02-Normals/" , "{seq_type}/{capture_space}--{genome_build}/{normal_id}-{pon}/{normal_id}.vcf.gz.tbi")
    resources: 
        **config["lcr-modules"]["make_pon"]["resources"]["concatenate_vcf"]
    conda:
        config["lcr-modules"]["make_pon"]["bcftools_env"]
    threads: 1
    log:
        os.path.join(LOGSDIR, "concat_norm_vcf", "{seq_type}/{capture_space}--{genome_build}/{pon}/{normal_id}.log")
    shell: 
        """
            bcftools concat {input.vcf} -Oz -o {output.vcf} && 
            tabix -p vcf {output.vcf} &> {log}
        """

def mutect2_get_chr_stats_norm(wildcards):
    return expand(expand(str(rules.make_pon_mutect2_germline.output.stats),zip,
    seq_type = wildcards.seq_type, 
    genome_build = wildcards.genome_build, 
    capture_space = wildcards.capture_space, 
    normal_id = wildcards.normal_id,
    pon = wildcards.pon, allow_missing = True
    ), interval = list(TARGET_INTERVALS.keys()))

rule make_pon_mutect2_merge_stats_norms:
    input:
        stats = mutect2_get_chr_stats_norm
    output:
        stats = os.path.join(WORKDIR,"02-Normals/" ,"{seq_type}/{capture_space}--{genome_build}/{normal_id}-{pon}/{normal_id}.stats")
    log:
        os.path.join(LOGSDIR, "tum_merge_stats","{seq_type}--{genome_build}/mutect2/{capture_space}/{normal_id}-{pon}/stats.log")
    conda: config["lcr-modules"]["make_pon"]["gatk_env"]
    resources: **config["lcr-modules"]["make_pon"]["resources"]["mutect"]
    shell:
        """
        gatk MergeMutectStats $(for i in {input.stats}; do echo -n "-stats $i "; done) -O {output.stats} > {log} 2>&1
        """

def get_normals_vcfs(wildcards):
    return expand(
        str(rules.make_pon_mutect2_concat_norm_vcf.output.vcf), zip,
        normal_id = SAMPLESHEET["sample_id"], 
        seq_type = SAMPLESHEET["seq_type"], 
        genome_build = SAMPLESHEET["genome_build"], 
        capture_space = SAMPLESHEET["capture_space"],
        pon = SAMPLESHEET["pon_name"]
    )


rule generate_samples_map:
    input:
        normal = get_normals_vcfs
    output:
        map_sample = os.path.join(WORKDIR, "{seq_type}/{capture_space}--{genome_build}/{pon}_samples_map.txt"),
    shell:
        """
            for samples in {input.normal}
            do
                name=$(basename $samples)
                name=${{name/.vcf.gz}}
                echo -e "$name\t$samples" >> {output.map_sample}
            done
        """

rule create_genomedb:
    input:
        sample_map = rules.generate_samples_map.output.map_sample
    output:
        somaticdb = directory(os.path.join(WORKDIR, "03-GDB", "{seq_type}/{capture_space}--{genome_build}/{pon}/somaticdb/")),
        done = touch(os.path.join(WORKDIR, "03-GDB", "{seq_type}/{capture_space}--{genome_build}/{pon}/makedb.done"))
    conda:
        config["lcr-modules"]["make_pon"]["gatk_env"]
    params:
        pon_path = os.path.join(WORKDIR, "03-GDB", "{seq_type}/{capture_space}--{genome_build}/{pon}/somaticdb/"),
        ref_fasta = config["lcr-modules"]["make_pon"]["ref_fasta"],
        target_bed = config["lcr-modules"]["make_pon"]["target_bed"]
    log:
        os.path.join(LOGSDIR, "create_genomedb", "{seq_type}/{capture_space}--{genome_build}/{pon}/create_genomedb.log")
    resources:
        mem_mb = 6000
    threads: 1
    shell:
        """
        gatk GenomicsDBImport --java-options "-Xmx{resources.mem_mb}m" --genomicsdb-workspace-path {params.pon_path} --batch-size 50 \
        --reference {params.ref_fasta} --merge-input-intervals TRUE \
        --sample-name-map {input.sample_map} --overwrite-existing-genomicsdb-workspace TRUE -L {params.target_bed} --interval-padding 150 > {log} 2>&1
        """

rule somaticPON:
    input:
        db_done = rules.create_genomedb.output.done
    output:
        pon = os.path.join(WORKDIR, "99-PON", "{seq_type}/{capture_space}--{genome_build}/{pon}/somaticvcf/{capture_space}_pon.vcf.gz"),
    conda:
        config["lcr-modules"]["make_pon"]["gatk_env"]
    resources: **config["lcr-modules"]["make_pon"]["resources"]["mutect"]
    params:
        ref_fasta = config["lcr-modules"]["make_pon"]["ref_fasta"],
        pon_path = "gendb://" + os.path.join(WORKDIR, "03-GDB", "{seq_type}/{capture_space}--{genome_build}/{pon}/somaticdb/")
    log:
        os.path.join(LOGSDIR, "somaticPON", "{seq_type}/{capture_space}--{genome_build}/{pon}/somaticPON.log")
    threads: 1
    shell:
        """
        gatk CreateSomaticPanelOfNormals -V {params.pon_path} -R {params.ref_fasta} -O {output.pon}  -OVI TRUE > {log} 2>&1
        """

rule decompress_vcf:
    input:
        vcf = rules.somaticPON.output.pon
    output:
        vcf = os.path.join(WORKDIR, "99-PON", "{seq_type}/{capture_space}--{genome_build}/{pon}/somaticvcf/{capture_space}_pon.vcf"),
    shell:
        "gzip -d {input.vcf}"


rule vcf2maf_annotate:
    input:
        vcf = rules.decompress_vcf.output.vcf
    output:
        vep_vcf = temp(os.path.join(WORKDIR, "99-PON", "{seq_type}/{capture_space}--{genome_build}/{pon}/somaticvcf/{capture_space}_pon.vep.vcf")),
        maf = os.path.join(WORKDIR, "99-PON", "{seq_type}/{capture_space}--{genome_build}/{pon}/somaticvcf/{capture_space}_pon.maf"),
    params:
        custom_enst = os.path.join(config["lcr-modules"]["make_pon"]["repo_path"], "resources/custom_enst.hg38.txt"),
        vep_data = config["lcr-modules"]["make_pon"]["vep_data"],
        centre = config["lcr-modules"]["make_pon"]["centre"],
        ref_fasta = config["lcr-modules"]["make_pon"]["ref_fasta"],
        ref_ver = config["lcr-modules"]["make_pon"]["genome_build"]
    resources:
        mem_mb = 8000
    threads: 2
    conda:
        "envs/vcf2maf.yaml"
    log:
        os.path.join(LOGSDIR, "vcf2maf","{seq_type}/{capture_space}--{genome_build}/{pon}.vcf2maf.log")
    shell:
        """
        vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} \
        --vep-path $CONDA_PREFIX/bin/ \
        --vep-data {params.vep_data} --vep-forks {threads} \
        --custom-enst {params.custom_enst} --ref-fasta {params.ref_fasta} \
        --species homo_sapiens --ncbi-build {params.ref_ver} \
        --maf-center {params.centre} 2> {log} >> {log}
        """

rule makeSomaticPON:
    input:
        expand(str(rules.vcf2maf_annotate.output.maf), zip, seq_type = SAMPLESHEET["seq_type"], 
                                                    capture_space = SAMPLESHEET["capture_space"], 
                                                    genome_build = SAMPLESHEET["genome_build"], 
                                                    pon = SAMPLESHEET["pon_name"])