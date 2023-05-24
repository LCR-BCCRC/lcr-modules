##### HEADER #####


include: "reference_files_header.smk"


##### SEQUENCE AND INDICES #####


rule get_genome_fasta_download:
    input: 
        fasta = rules.download_genome_fasta.output.fasta
    output: 
        fasta = "genomes/{genome_build}/genome_fasta/genome.fa"
    wildcard_constraints:
        genome_build = ".+(?<!masked)"
    conda: CONDA_ENVS["coreutils"]
    shell:
        "ln -srf {input.fasta} {output.fasta}"


rule index_genome_fasta:
    input: 
        fasta = rules.get_genome_fasta_download.output.fasta
    output: 
        fai = "genomes/{genome_build}/genome_fasta/genome.fa.fai"
    log: 
        "genomes/{genome_build}/genome_fasta/genome.fa.fai.log"
    wildcard_constraints:
        genome_build = ".+(?<!masked)"
    conda: CONDA_ENVS["samtools"]
    shell:
        "samtools faidx {input.fasta} > {log} 2>&1"


rule create_bwa_index:
    input: 
        fasta = rules.get_genome_fasta_download.output.fasta
    output: 
        prefix = touch("genomes/{genome_build}/bwa_index/bwa-{bwa_version}/genome.fa")
    log: 
        "genomes/{genome_build}/bwa_index/bwa-{bwa_version}/genome.fa.log"
    conda: CONDA_ENVS["bwa"]
    resources:
        mem_mb = 20000
    shell:
        "bwa index -p {output.prefix} {input.fasta} > {log} 2>&1"


rule create_gatk_dict:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta,
        fai = rules.index_genome_fasta.output.fai
    output:
        dict = "genomes/{genome_build}/genome_fasta/genome.dict"
    log:
        "genomes/{genome_build}/gatk_fasta/genome.dict.log"
    conda: CONDA_ENVS["gatk"]
    resources:
        mem_mb = 20000
    shell:
        op.as_one_line(""" 
        gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict} > {log} 2>&1
        """)


rule create_star_index:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta,
        gtf = get_download_file(rules.download_gencode_annotation.output.gtf)
    output: 
        index = directory("genomes/{genome_build}/star_index/star-{star_version}/gencode-{gencode_release}/overhang-{star_overhang}")
    log: 
        "genomes/{genome_build}/star_index/star-{star_version}/gencode-{gencode_release}/overhang-{star_overhang}.log"
    conda: CONDA_ENVS["star"]
    threads: 12
    resources:
        mem_mb = 42000
    shell:
        op.as_one_line("""
        mkdir -p {output.index}
            &&
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.index}
        --genomeFastaFiles {input.fasta} --sjdbOverhang {wildcards.star_overhang}
        --sjdbGTFfile {input.gtf} --outTmpDir {output.index}/_STARtmp
        --outFileNamePrefix {output.index}/ > {log} 2>&1
        """)


rule get_liftover_chains:
    input:
        chains = rules.download_liftover_chains.output.chains
    output:
        chains = "genomes/{genome_build}/chains/{version}/{chain_version}.over.chain"
    shell:
        "ln -srf {input.chains} {output.chains}"

rule get_sdf_refs: 
    input: 
        sdf = ancient(rules.download_sdf.output.sdf)
    output: 
        sdf = directory("genomes/{genome_build}/sdf")
    wildcard_constraints: 
        genome_build = "|".join(SDF_GENOME_BUILDS)
    shell: 
        "ln -srfT {input.sdf} {output.sdf}"


rule get_masked_genome_fasta_download:
    input: 
        fasta = rules.download_masked_genome_fasta.output.fasta
    output: 
        fasta = "genomes/{genome_build}/genome_fasta/genome.fa"
    wildcard_constraints:
        genome_build = ".+_masked"
    conda: CONDA_ENVS["coreutils"]
    shell:
        "ln -srf {input.fasta} {output.fasta}"


rule index_masked_genome_fasta:
    input: 
        fasta = rules.get_masked_genome_fasta_download.output.fasta
    output: 
        fai = "genomes/{genome_build}/genome_fasta/genome.fa.fai"
    log: 
        "genomes/{genome_build}/genome_fasta/genome.fa.fai.log"
    wildcard_constraints:
        genome_build = ".+_masked"
    conda: CONDA_ENVS["samtools"]
    shell:
        "samtools faidx {input.fasta} > {log} 2>&1"


##### METADATA #####


rule store_genome_build_info:
    output: 
        version = "genomes/{genome_build}/version.txt",
        provider = "genomes/{genome_build}/provider.txt"
    params:
        version = lambda w: config["genome_builds"][w.genome_build]["version"],
        provider = lambda w: config["genome_builds"][w.genome_build]["provider"]
    shell: 
        op.as_one_line("""
        echo "{params.version}" > {output.version}
            &&
        echo "{params.provider}" > {output.provider}
        """)


rule get_main_chromosomes_download:
    input: 
        txt = get_download_file(rules.download_main_chromosomes.output.txt),
        chrx = get_download_file(rules.download_chromosome_x.output.txt),
        fai = rules.index_genome_fasta.output.fai
    output: 
        txt = "genomes/{genome_build}/genome_fasta/main_chromosomes.txt",
        bed = "genomes/{genome_build}/genome_fasta/main_chromosomes.bed",
        chrx = "genomes/{genome_build}/genome_fasta/chromosome_x.txt",
        patterns = temp("genomes/{genome_build}/genome_fasta/main_chromosomes.patterns.txt")
    conda: CONDA_ENVS["coreutils"]
    shell: 
        op.as_one_line("""
        sed 's/^/^/' {input.txt} > {output.patterns}
            &&
        egrep -w -f {output.patterns} {input.fai}
            |
        cut -f1 > {output.txt}
            &&
        egrep -w -f {output.patterns} {input.fai}
            |
        awk 'BEGIN {{FS=OFS="\t"}} {{print $1,  0, $2}}' > {output.bed}
            &&
        ln -srf {input.chrx} {output.chrx}
        """)


rule get_main_chromosomes_withY_download:
    input: 
        txt = get_download_file(rules.download_main_chromosomes_withY.output.txt),
        fai = rules.index_genome_fasta.output.fai
    output: 
        txt = "genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt",
        bed = "genomes/{genome_build}/genome_fasta/main_chromosomes_withY.bed",
        patterns = temp("genomes/{genome_build}/genome_fasta/main_chromosomes.patterns.txt")
    conda: CONDA_ENVS["coreutils"]
    shell: 
        op.as_one_line("""
        sed 's/^/^/' {input.txt} > {output.patterns}
            &&
        egrep -w -f {output.patterns} {input.fai}
            |
        cut -f1 > {output.txt}
            &&
        egrep -w -f {output.patterns} {input.fai}
            |
        awk 'BEGIN {{FS=OFS="\t"}} {{print $1,  0, $2}}' > {output.bed}
        """)

##### ANNOTATIONS #####


rule get_gencode_download: 
    input:
        gtf = get_download_file(rules.download_gencode_annotation.output.gtf)
    output:
        gtf = "genomes/{genome_build}/annotations/gencode_annotation-{gencode_release}.gtf"
    conda: CONDA_ENVS["coreutils"]
    shell:
        "ln -srf {input.gtf} {output.gtf}"


rule calc_gc_content:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta
    output:
        wig = "genomes/{genome_build}/annotations/gc_wiggle.window_{gc_window_size}.wig.gz"
    log:
        "genomes/{genome_build}/annotations/gc_wiggle.window_{gc_window_size}.wig.gz.log"
    conda: CONDA_ENVS["sequenza-utils"]
    shell:
        op.as_one_line("""
        sequenza-utils gc_wiggle --fasta {input.fasta} -w {wildcards.gc_window_size} -o -
            |
        gzip -c > {output.wig}
        """)
rule get_jabba_gc_rds:
    input:
        bed = get_download_file(rules.download_ucsc_gc.output.bed),
        txt = get_download_file(rules.download_ucsc_chrom_sizes.output.txt)
    output:
        rds = "genomes/{genome_build}/annotations/jabba/gc1000.rds"
    conda: CONDA_ENVS["rtracklayer"]
    shell:
        op.as_one_line("""
        Rscript
            -e 'library(rtracklayer); gr <- import("{input.bed}")'
            -e 'names(mcols(gr))[1] <- "score"; gr[gr$score == "nan",]$score <- 0'
            -e 'gr$score <- as.numeric(gr$score)'
            -e 'sizes <- read.delim("{input.txt}", stringsAsFactors = FALSE, col.names = c("chrom", "length"), header = FALSE)'
            -e 'sizes <- unlist(split(as.numeric(sizes$length), sizes$chrom))[levels(seqnames(gr))]'
            -e 'seqlengths(gr) <- sizes'
            -e 'saveRDS(gr, "{output.rds}")'
        """)

rule get_jabba_map_rds:
    input:
        bed = get_download_file(rules.download_ucsc_map.output.bed),
        txt = get_download_file(rules.download_ucsc_chrom_sizes.output.txt)
    output:
        rds = "genomes/{genome_build}/annotations/jabba/map1000.rds"
    conda: CONDA_ENVS["rtracklayer"]
    shell:
        op.as_one_line("""
        Rscript
            -e 'library(rtracklayer); gr <- import("{input.bed}")'
            -e 'names(mcols(gr))[1] <- "score"; gr[gr$score == "nan",]$score <- 0'
            -e 'gr$score <- as.numeric(gr$score)'
            -e 'sizes <- read.delim("{input.txt}", stringsAsFactors = FALSE, col.names = c("chrom", "length"), header = FALSE)'
            -e 'sizes <- unlist(split(as.numeric(sizes$length), sizes$chrom))[levels(seqnames(gr))]'
            -e 'seqlengths(gr) <- sizes'
            -e 'saveRDS(gr, "{output.rds}")'
        """)

rule get_par_rds:
    input:
        bed = get_download_file(rules.download_par_bed.output.bed),
        txt = get_download_file(rules.download_ucsc_chrom_sizes.output.txt)
    output:
        rds = "genomes/{genome_build}/annotations/jabba/PAR_{genome_build}.rds"
    conda: CONDA_ENVS["rtracklayer"]
    shell:
        op.as_one_line(""" 
        Rscript 
            -e 'library(rtracklayer); gr <- import("{input.bed}")'
            -e 'sizes <- read.delim("{input.txt}", stringsAsFactors = FALSE, col.names = c("chrom", "length"), header = FALSE)'
            -e 'sizes <- unlist(split(as.numeric(sizes$length), sizes$chrom))[levels(seqnames(gr))]'
            -e 'seqlengths(gr) <- sizes'
            -e 'saveRDS(gr, "{output.rds}")'
        """)




##### VARIATION #####


rule get_dbsnp_download: 
    input:
        vcf = get_download_file(rules.download_dbsnp_vcf.output.vcf),
        fai = str(rules.index_genome_fasta.output.fai),
        bed = str(rules.get_main_chromosomes_withY_download.output.bed)
    output:
        vcf = "genomes/{genome_build}/variation/dbsnp.common_all-{dbsnp_build}.vcf.gz",
        tmpfile = temp("genomes/{genome_build}/variation/dbsnp.common_all-{dbsnp_build}.vcf.tmp")
    conda: CONDA_ENVS["bcftools"]
    shell:
        op.as_one_line("""
        zgrep -v '##contig' {input.vcf} > {output.tmpfile} &&
        bcftools reheader --fai {input.fai} {output.tmpfile} | bcftools view -T {input.bed} -O z -o {output.vcf} &&
        bcftools index -t {output.vcf}
        """)

##### PICARD METRICS

rule create_rRNA_interval:
    input:
        gtf = rules.get_gencode_download.output.gtf
    output: 
        rrna_int = "genomes/{genome_build}/rrna_intervals/rRNA_int_gencode-{gencode_release}.txt"
    log: 
        "genomes/{genome_build}/rrna_intervals/rRNA_int_gencode-{gencode_release}.log"
    conda: CONDA_ENVS["picard"]
    shell:
        op.as_one_line("""
        grep 'gene_type "rRNA"' {input.gtf} |
        awk '$3 == "transcript"' |
        cut -f1,4,5,7,9 |
        perl -lane '
            /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
            print join "\t", (@F[0,1,2,3], $1)
        ' | 
        sort -k1V -k2n -k3n >> {output.rrna_int}
        &&
        chmod a-w {output.rrna_int}
        """)


rule create_refFlat:
    input:
        gtf = rules.get_gencode_download.output.gtf
    output:
        txt = "genomes/{genome_build}/annotations/refFlat_gencode-{gencode_release}.txt"
    log: "genomes/{genome_build}/annotations/gtfToGenePred-{gencode_release}.log"
    conda: CONDA_ENVS["ucsc-gtftogenepred"]
    threads: 4
    resources:
        mem_mb = 6000
    shell:
        op.as_one_line("""
        gtfToGenePred -genePredExt -geneNameAsName2 
        {input.gtf} {output.txt}.tmp 
        2> {log} &&
        paste <(cut -f 12 {output.txt}.tmp) <(cut -f 1-10 {output.txt}.tmp) > {output.txt}
        """)

##### ENCODE #####

rule get_blacklist_download: 
    input: 
        bed = get_download_file(rules.download_blacklist.output.bed)
    output: 
        bed = "genomes/{genome_build}/encode/encode-blacklist.{genome_build}.bed"
    conda: CONDA_ENVS["coreutils"]
    shell:
        "ln -srf {input.bed} {output.bed}"

##### REPEATMASKER #####

rule get_repeatmasker_download:
    input: 
        bed = get_download_file(rules.download_repeatmasker.output.bed)
    output: 
        bed = "genomes/{genome_build}/repeatmasker/repeatmasker.{genome_build}.bed"
    conda: CONDA_ENVS["coreutils"]
    shell: 
        "ln -srf {input.bed} {output.bed}"

# salmon index
rule _download_salmon_script:
    output: 
        script = "downloads/scripts/salmon/generateDecoyTranscriptome.sh"
    log: 
        "downloads/scripts/salmon/log"
    params: 
        url = "https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh"
    shell:
        op.as_one_line("""
        curl -L {params.url} > {output.script} 2> {log}
            &&
        chmod a-w {output.script}
        """)


rule _create_transcriptome_fasta:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta,
        gtf = get_download_file("downloads/gencode-33/gencode.annotation.{version}.gtf")
    output: 
        fasta = "genomes/{genome_build}/salmon_index/salmon-{salmon_version}/transcriptome.fa"
    log: 
        "genomes/{genome_build}/salmon_index/salmon-{salmon_version}/transcriptome.log"
    conda: CONDA_ENVS["gffread"]
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        op.as_one_line("""
        gffread -F
        -w {output.fasta}
        -g {input.fasta}
        {input.gtf}
        > {log} 2>&1
            &&
        chmod a-w {output.fasta}
        """)

rule create_salmon_index:
    input:
        fasta = rules._create_transcriptome_fasta.output.fasta,
        gtf = get_download_file("downloads/gencode-33/gencode.annotation.{version}.gtf")
    output: 
        index = directory("genomes/{genome_build}/salmon_index/salmon-{salmon_version}/index")
    log: 
        "genomes/{genome_build}/salmon_index/salmon-{salmon_version}/log"
    conda: CONDA_ENVS["salmon"]
    threads: 8
    resources:
        mem_mb = 12000
    shell:
        op.as_one_line("""
        mkdir {output.index}
            &&
        chmod u+w {output.index}
            &&
        salmon index
        --threads {threads}
        -t {input.fasta}
        -i {output.index}
        > {log} 2>&1
        """)

rule get_af_only_gnomad_vcf:
    input:
        vcf = get_download_file(rules.download_af_only_gnomad_vcf.output.vcf),
        fai = str(rules.index_genome_fasta.output.fai),
        bed = str(rules.get_main_chromosomes_withY_download.output.bed)
    output:
        vcf = "genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz",
        tmpfile = temp("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.tmp")
    conda: CONDA_ENVS["bcftools"]
    shell:
        op.as_one_line("""
        zgrep -v '##contig' {input.vcf} > {output.tmpfile} &&
        bcftools reheader --fai {input.fai} {output.tmpfile} | bcftools view -T {input.bed} -O z -o {output.vcf} &&
        bcftools index -t {output.vcf}
        """)

rule normalize_af_only_gnomad_vcf:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta,
        vcf = str(rules.get_af_only_gnomad_vcf.output.vcf)
    output:
        vcf = "genomes/{genome_build}/variation/af-only-gnomad.normalized.{genome_build}.vcf.gz",
        vcf_index = "genomes/{genome_build}/variation/af-only-gnomad.normalized.{genome_build}.vcf.gz.tbi"
    conda: CONDA_ENVS["bcftools"]
    shell:
        op.as_one_line("""
        bcftools view {input.vcf} | grep -v "_alt" | bcftools norm -m -any -f {input.fasta} | bgzip -c > {output.vcf}
            &&
        bcftools index -t {output.vcf}
            &&
        touch {output.vcf_index}
        """)

rule get_mutect2_pon:
    input:
        vcf = get_download_file(rules.download_mutect2_pon.output.vcf),
        fai = str(rules.index_genome_fasta.output.fai),
        bed = str(rules.get_main_chromosomes_withY_download.output.bed)
    output:
        vcf = "genomes/{genome_build}/gatk/mutect2_pon.{genome_build}.vcf.gz",
        tmpfile = temp("genomes/{genome_build}/gatk/mutect2_pon.{genome_build}.vcf.tmp")
    conda: CONDA_ENVS["bcftools"]
    log:
        "genomes/{genome_build}/gatk/mutect2_pon.{genome_build}.vcf.log"
    shell:
        op.as_one_line(""" 
        zgrep -v '##contig' {input.vcf} > {output.tmpfile} &&
        bcftools reheader --fai {input.fai} {output.tmpfile} | bcftools view -T {input.bed} -O z -o {output.vcf} 2> {log} &&
        bcftools index -t {output.vcf}
        """)

rule get_mutect2_small_exac:
    input:
        vcf = get_download_file(rules.download_mutect2_small_exac.output.vcf),
        fai = str(rules.index_genome_fasta.output.fai),
        bed = str(rules.get_main_chromosomes_withY_download.output.bed)
    output:
        vcf = "genomes/{genome_build}/gatk/mutect2_small_exac.{genome_build}.vcf.gz",
        tmpfile = temp("genomes/{genome_build}/gatk/mutect2_small_exac.{genome_build}.vcf.tmp")
    log:
        "genomes/{genome_build}/gatk/mutect2_pon.{genome_build}.vcf.log"
    conda: CONDA_ENVS["bcftools"]
    shell:
        op.as_one_line(""" 
        zgrep -v '##contig' {input.vcf} > {output.tmpfile} &&
        bcftools reheader --fai {input.fai} {output.tmpfile} | bcftools view -T {input.bed} -O z -o {output.vcf} 2> {log} &&
        bcftools index -t {output.vcf}
        """)


### FOR HANDLING CAPTURE SPACE/TARGETED SEQUENCING DATA ###
# Added by Chris
# What this *should* do (if I have written these rules correctly) is obtain a capture space BED file
# check the contig names against the reference (and replace them as necessary), pad, sort,merge, and generate a bgzip/tabix
# indexed version of the BED as well as an interval list.

def _check_capspace_provider(w):
    # Checks to determine if both the capture space BED and associated reference genome
    # are chr prefixed or not

    # If this is specified as the "default", make sure we obtain the relevent capture space
    default_key = "default-" + w.genome_build
    if default_key not in config["capture_space"] and w.capture_space == default_key:
        # aka if the user hasn't explicitly specified a default using a default name
        capture_space = _get_default_capspace(w)
    else:
        capture_space = w.capture_space

    genome_provider = config["genome_builds"][w.genome_build]["provider"]
    bed_provider = config["capture_space"][capture_space]["provider"]
    
    # If the providers match (i.e. they share the same prefix), just use the downloaded version
    if genome_provider == bed_provider:
        return {'bed': expand(rules.download_capspace_bed.output.capture_bed, capture_space=capture_space, genome_build=w.genome_build)}
    else:
        # Prompt the prefix to be converted
        chr_status = "chr" if genome_provider == "ucsc" else "no_chr"
        return {'bed': expand(rules.add_remove_chr_prefix_bed.output.converted_bed, capture_space = capture_space, genome_build = w.genome_build, chr_status= chr_status)}


rule get_capspace_bed_download:
    input:
        unpack(_check_capspace_provider)
    output:
        capture_bed = "genomes/{genome_build}/capture_space/{capture_space}.bed"
    conda: CONDA_ENVS["coreutils"]
    shell:
        "ln -srf {input.bed} {output.capture_bed}"


rule sort_and_pad_capspace:
    input:
       capture_bed = rules.get_capspace_bed_download.output.capture_bed,
       fai = rules.index_genome_fasta.output.fai
    output:
        intermediate_bed = temp("genomes/{genome_build}/capture_space/{capture_space}.intermediate.bed"),
        padded_bed = "genomes/{genome_build}/capture_space/{capture_space}.padded.bed"
    params:
        padding_size = config['capture_params']['padding_size']  # Default to 200. I would be warry of changing
    log:
        "genomes/{genome_build}/capture_space/{capture_space}.padded.bed.log"
    conda: CONDA_ENVS["bedtools"]
    shell:
        op.as_one_line("""
        cat {input.fai} | cut -f 1-2 | perl -ane 'print "$F[0]\\t0\\t$F[1]\\n"' | bedtools intersect -wa -a {input.capture_bed} -b stdin > {output.intermediate_bed}
            &&
        bedtools slop -b {params.padding_size} -i {output.intermediate_bed} -g {input.fai} | bedtools sort | bedtools merge > {output.padded_bed} 2> {log}
        """)

rule check_capspace_contigs:
    input:
        bed = rules.get_capspace_bed_download.output.capture_bed,
        fai = rules.index_genome_fasta.output.fai
    output:
        contig_log = "genomes/{genome_build}/capture_space/{capture_space}.check_contigs.log"
    run:
        # Parse BED contigs
        bed_contigs = set()
        with open(input.bed) as f:
            for line in f:
                line = line.rstrip()
                contig = line.split("\t")[0]
                if contig not in bed_contigs:
                    bed_contigs.add(contig)

        # Parse fai contigs
        fai_contigs = set()
        with open(input.fai) as f:
            for line in f:
                line = line.rstrip()
                contig = line.split('\t')[0]
                if contig not in fai_contigs:
                    fai_contigs.add(contig)

        # Check the BED file for contigs that are not in the reference genome
        missing_contigs = list(x for x in bed_contigs if x not in fai_contigs)
        with open(output.contig_log, "w") as o:
            if len(missing_contigs) == 0:
                o.write("No contigs missing from reference")
            else:
                o.write("The following contigs were missing from the reference\n")
                o.write("\n".join(missing_contigs))
            o.write("\n")


rule compress_index_capspace_bed:
    input:
        capture_bed = rules.sort_and_pad_capspace.output.padded_bed
    output:
        bgzip_bed = "genomes/{genome_build}/capture_space/{capture_space}.padded.bed.gz",
        tabix = "genomes/{genome_build}/capture_space/{capture_space}.padded.bed.gz.tbi"
    log:
        "genomes/{genome_build}/capture_space/{capture_space}.padded.bed.gz.log"
    conda: CONDA_ENVS["tabix"]
    shell:
        op.as_one_line("""
        bgzip -c {input.capture_bed} > {output.bgzip_bed}
            &&
        tabix -p bed {output.bgzip_bed}
        """)


rule create_interval_list:
    input:
        bed = rules.sort_and_pad_capspace.output.padded_bed,
        sd = rules.create_gatk_dict.output.dict
    output:
        interval_list = "genomes/{genome_build}/capture_space/{capture_space}.padded.interval_list"
    log:
        "genomes/{genome_build}/capture_space/{capture_space}.padded.interval_list.log"
    conda: CONDA_ENVS["gatk"]
    shell:
        "gatk BedToIntervalList --INPUT {input.bed} -SD {input.sd} -O {output.interval_list} > {log} 2>&1"

##### SigProfiler #####

rule download_sigprofiler_genome:
    output:
        complete = "downloads/sigprofiler_prereqs/{sigprofiler_build}.installed"
    conda: CONDA_ENVS["sigprofiler"]
    shell:
        op.as_one_line("""
        python -c 'from SigProfilerMatrixGenerator import install as genInstall;
        genInstall.install("{wildcards.sigprofiler_build}", rsync = False, bash = True)'
            &&
        touch {output.complete}
        """)

def get_sigprofiler_genome(wildcards):
    sigprofiler_build = ''
    if wildcards.genome_build in ['grch37','hg19','hs37d5']:
        sigprofiler_build = "GRCh37"
    elif wildcards.genome_build in ['grch38','grch38-legacy','hg38','hg38-panea']:
        sigprofiler_build = "GRCh38"
    return("downloads/sigprofiler_prereqs/" + sigprofiler_build + ".installed")

rule install_sigprofiler_genome:
    input:
        get_sigprofiler_genome
    output:
        complete = "genomes/{genome_build}/sigprofiler_genomes/{genome_build}.installed"
    run:
        op.relative_symlink(input, output.complete)


##### Non-B DNA structure #####

rule combine_gff_nonb:
    input:
        complete = str(rules.download_nonb_dna_gff.output.complete)
    params:
        gff = lambda w: expand("downloads/nonb_dna/gff/{genome_build}/chr{chr}_{Btype}.gff", chr=[str(i) for i in range(1,22)] + ["MT","X","Y"], Btype=["Z","MR","IR","GQ","DR","APR"], allow_missing=True)
    output:
        gff = "genomes/{version}/nonb_dna/gff/{genome_build}.gff",
        complete = "downloads/nonb_dna/gff/{version}.{genome_build}.combine.complete"
    shell:
        op.as_one_line("""
            tail -n +2 -q {params.gff} > {output.gff} &&
            touch {output.complete}
        """)

rule gff_nonb_cvbio:
    input:
        gff = str(rules.combine_gff_nonb.output.gff)
    output:
        bed = "genomes/{version}/nonb_dna/{genome_build}.bed"
    conda:
        CONDA_ENVS["bedops"]
    shell:
        op.as_one_line("""
            gff2bed < {input.gff} > {output.bed}
        """)