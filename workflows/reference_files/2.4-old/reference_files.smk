#!/usr/bin/env snakemake


##### MODULES #####


import os
import re
import tarfile
import urllib.request
from collections import defaultdict

import yaml
import snakemake as smk
from snakemake.logging import logger

import oncopipe as op


##### CONFIG #####
localrules: download_genome_fasta,
            download_main_chromosomes, download_gencode_annotation,
            hardlink_download, update_contig_names,
            get_genome_fasta_download, index_genome_fasta,
            get_main_chromosomes_download, create_bwa_index,
            get_gencode_download, create_star_index


# Check for genome builds
assert "genome_builds" in config and len(config["genome_builds"]) > 0, (
    "No `genome_builds` in snakemake configuration."
)

# Switch between case for version names
VERSION_UPPER = {
    "grch37": "GRCh37",
    "GRCh37": "GRCh37",
    "grch38": "GRCh38",
    "GRCh38": "GRCh38",
}

# Check genome build versions, providers, and genome_fasta
possible_versions = list(VERSION_UPPER.keys())
possible_providers = ["ensembl", "ucsc", "gencode", "ncbi"]
for build_name, build_info in config["genome_builds"].items():
    assert "version" in build_info and build_info["version"] in possible_versions, (
        f"`version` not set for `{build_name}` or `version` not among {possible_versions}."
    )
    assert "provider" in build_info and build_info["provider"] in possible_providers, (
        f"`provider` not set for `{build_name}` or `provider` not among {possible_providers}."
    )
    assert "genome_fasta_url" in build_info, f"`genome_fasta_url` not set for `{build_name}`."
    url_code = urllib.request.urlopen(build_info["genome_fasta_url"]).getcode()
    assert url_code == 200, (
        f"Pinging `genome_fasta_url` for {build_name} returned HTTP code {url_code} "
        f"(rather than 200): \n{build_info['genome_fasta_url']}"
    )    


##### TOOLS #####


CONDA_ENVS = { pkg: pkg_info["conda_env"] for pkg, pkg_info in config["tools"].items() }

TOOL_VERSIONS = { pkg: pkg_info["version"] for pkg, pkg_info in config["tools"].items() }


##### WILDCARD CONSTRAINTS


wildcard_constraints:
    genome_build = "|".join(config["genome_builds"].keys()),
    version = "|".join(VERSION_UPPER.keys()),
    bwa_version = TOOL_VERSIONS["bwa"],
    star_version = TOOL_VERSIONS["star"],
    gencode_release = "|".join(config["wildcard_values"]["gencode_release"]),


##### CHROMOSOME MAPPINGS #####


# Define method for download ChromosomeMappings repository
def download_chrom_mappings(chrom_mappings_dir):
    chrom_mappings_url = "https://github.com/BrunoGrandePhD/ChromosomeMappings/archive/master.tar.gz"
    tar_file_path = "ChromosomeMappings.tar.gz"
    os.makedirs("./", exist_ok=True)
    urllib.request.urlretrieve(chrom_mappings_url, tar_file_path)
    zip_file = tarfile.open(tar_file_path, "r:gz")
    zip_file.extractall("./")
    zip_file.close()
    os.remove(tar_file_path)
    os.rename("ChromosomeMappings-master", chrom_mappings_dir)

# Download ChromosomeMappings if not already done
CHROM_MAPPINGS_DIR = "chrom_mappings"
if not os.path.exists(CHROM_MAPPINGS_DIR):
    download_chrom_mappings(CHROM_MAPPINGS_DIR)

# Find chromosome mappings for human genome
CHROM_MAPPINGS_FILES = os.listdir(CHROM_MAPPINGS_DIR)
CHROM_MAPPINGS_FILES = [ f for f in CHROM_MAPPINGS_FILES if f.startswith("GRCh") ]

TO_PROVIDERS = set()
CHROM_MAPPINGS = defaultdict(lambda: defaultdict(set))
for chrom_map_file in CHROM_MAPPINGS_FILES:
    chrom_map_file_split = re.split("[_2.]", chrom_map_file)
    version, from_provider, to_provider, _ = [ x.lower() for x in chrom_map_file_split ]
    CHROM_MAPPINGS[version][to_provider].add(from_provider)
    # Include the provider itself since you can always just use a file as is
    CHROM_MAPPINGS[version][to_provider].add(to_provider)
    TO_PROVIDERS.add(to_provider)


##### DOWNLOAD #####


rule download_genome_fasta:
    output: 
        fasta = "downloads/genome_fasta/{genome_build}.fa"
    log: 
        "downloads/genome_fasta/{genome_build}.fa.log"
    params: 
        url = lambda w: config["genome_builds"][w.genome_build]["genome_fasta_url"]
    shell:
        op.as_one_line("""
        curl -L {params.url} > {output.fasta} 2> {log}
            &&
        chmod a-w {output.fasta}
        """)


rule download_main_chromosomes:
    input:
        mapping = lambda w: f"{CHROM_MAPPINGS_DIR}/{VERSION_UPPER[w.version]}_ensembl2ucsc.txt"
    output:
        txt = "downloads/main_chromosomes/main_chromosomes.{version}.txt"
    params:
        provider = "ensembl"
    shell:
        op.as_one_line("""
        egrep -w "^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|MT)" {input.mapping}
            |
        cut -f1 > {output.txt}
            &&
        chmod a-w {output.txt}
        """)


rule download_gencode_annotation:
    output:
        gtf = "downloads/gencode-{gencode_release}/gencode.annotation.{version}.gtf"
    params:
        provider = "ucsc"
    run:
        url_parts = [
            "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human",
            f"release_{wildcards.gencode_release}"
        ]
        release_fmt = f"v{wildcards.gencode_release}"
        version = VERSION_UPPER[wildcards.version]
        if version == "GRCh37":
            url_parts.append("GRCh37_mapping")
            release_fmt += "lift37"
        url_parts.append(f"gencode.{release_fmt}.annotation.gtf.gz")
        url = "/".join(url_parts)
        urllib.request.urlretrieve(url, output.gtf + ".gz")
        shell("gunzip {output.gtf}.gz")
        shell("chmod a-w {output.gtf}")


##### FUNCTIONS #####


def get_matching_download_rules(file):
    ignored_rules = ["download_genome_fasta"]
    rule_names = [ r for r in dir(rules) if r.startswith("download_")]
    rule_names = [ r for r in rule_names if r not in ignored_rules ]
    rule_list = [ getattr(rules, name) for name in rule_names ]
    matching_rules = []
    for rule_name, r in zip(rule_names, rule_list):
        # Ensure provider is specified
        assert "provider" in r.params._names, (
                f"The `{rule_name}` download rule doesn't have a `provider` param."
            )
        # At least one output file should produce
        num_matches = []
        for output_file in r.output:
            assert "{version}" in output_file, (
                f"The `{rule_name}` download rule doesn't have a `{{version}}` "
                f"wildcard in the output file ('{output_file}')."
            )
            matches = smk.io.glob_wildcards(output_file, [file])
            num_matches.append(len(matches[0]))
        if any(num > 0 for num in num_matches):
            matching_rules.append(r)
    return matching_rules


def hardlink_same_provider(wildcards):

    genome_build = wildcards.genome_build
    suffix = wildcards.suffix
    output_file = f"genomes/{genome_build}/{suffix}"
    raw_download_file = f"downloads/{suffix}"

    version = config["genome_builds"][genome_build]["version"]
    to_provider = config["genome_builds"][genome_build]["provider"]
    from_provider_options = CHROM_MAPPINGS[version][to_provider]
    
    dependencies = []
    matching_rules = get_matching_download_rules(raw_download_file)

    for r in matching_rules:
        # The provider must be among the ones we can convert from
        if r.params.provider not in from_provider_options:
            continue
        dependencies.append(r)

    if len(dependencies) == 0:
        msg = f"Could not find rule to generate {output_file}."
        raise Exception(msg)

    if len(dependencies) > 1:
        msg = f"Found conflicting rules to generate {output_file}."
        raise Exception(msg)

    raw_download_root, raw_download_ext = os.path.splitext(raw_download_file)
    target_download_file = f"{raw_download_root}.{to_provider}{raw_download_ext}"
    return target_download_file


def get_cvbio_params(field):

    def get_cvbio_params_custom(wildcards, input, output):

        dependencies = []
        matching_rules = get_matching_download_rules(input.before)
        assert len(matching_rules) == 1, "Not just one matching rule."

        dependency = matching_rules[0]
        version = VERSION_UPPER[wildcards.version]
        from_provider = dependency.params.provider
        to_provider = wildcards.to_provider
        file_ext = wildcards.ext

        cvbio_params = {
            "mapping": f"{CHROM_MAPPINGS_DIR}/{version}_{from_provider}2{to_provider}.txt",
            "from_provider": from_provider,
            "comment": config["cvbio_config"][file_ext]["comment"],
            "columns": config["cvbio_config"][file_ext]["columns"],
            "skip": config["cvbio_config"][file_ext]["skip"],
            "delimiter": config["cvbio_config"][file_ext]["delimiter"]
        }
        return cvbio_params[field]

    return get_cvbio_params_custom


def get_download_file(file):
    file_root, file_ext = os.path.splitext(file)
    assert file_ext.lstrip(".") in config["cvbio_config"], (
        f"`{file_ext}` files are not yet configured for cvbio. See `config['cvbio_config']` variable."
    )
    def get_download_file_custom(wildcards):
        genome_build = wildcards.genome_build
        version = config["genome_builds"][genome_build]["version"]
        provider = config["genome_builds"][genome_build]["provider"]
        download_file = file.replace("{version}", version).replace("{provider}", provider)
        download_file = "genomes/{genome_build}/" + download_file
        return download_file
    return get_download_file_custom


##### SHARED #####


rule hardlink_download:
    input: hardlink_same_provider
    output: "genomes/{genome_build}/downloads/{suffix}"
    shell: "ln -f {input} {output}"


rule update_contig_names:
    input:
        before = "downloads/{parent_dir}/{prefix}.{version}.{ext}"
    output:
        after = "downloads/{parent_dir}/{prefix}.{version}.{to_provider}.{ext}"
    log:
        "downloads/{parent_dir}/{prefix}.{version}.{to_provider}.{ext}.log"
    params:
        mapping = get_cvbio_params("mapping"),
        from_provider = get_cvbio_params("from_provider"),
        comment = get_cvbio_params("comment"),
        columns = get_cvbio_params("columns"),
        skip = get_cvbio_params("skip"),
        delimiter = get_cvbio_params("delimiter")
    wildcard_constraints:
        ext = "|".join(config["cvbio_config"].keys()),
        to_provider = "|".join(TO_PROVIDERS)
    conda: CONDA_ENVS["cvbio"]
    shell:
        op.as_one_line("""
        if [[ '{wildcards.to_provider}' == '{params[from_provider]}' ]]; then
            ln -f {input.before} {output.after};
        else
            cvbio UpdateContigNames --in {input.before} --out {output.after}
            --mapping {params.mapping} --comment-chars '{params.comment}'
            --columns {params.columns} --skip-missing {params.skip}
            --delimiter '{params.delimiter}' > {log} 2>&1
                &&
            chmod a-w {output.after}; 
        fi
        """)


##### GENOME BUILDS #####


rule get_genome_fasta_download:
    input: 
        fasta = rules.download_genome_fasta.output.fasta
    output: 
        fasta = "genomes/{genome_build}/genome_fasta/genome.fa"
    shell:
        "ln -f {input.fasta} {output.fasta}"


rule index_genome_fasta:
    input: 
        fasta = rules.get_genome_fasta_download.output.fasta
    output: 
        fai = "genomes/{genome_build}/genome_fasta/genome.fa.fai"
    log: 
        "genomes/{genome_build}/genome_fasta/genome.fa.fai.log"
    conda: CONDA_ENVS["samtools"]
    shell:
        op.as_one_line("""
        samtools faidx {input.fasta} > {log} 2>&1
            &&
        chmod a-w {output.fai}
        """)



rule get_main_chromosomes_download:
    input: 
        txt = get_download_file(rules.download_main_chromosomes.output.txt),
        fai = rules.index_genome_fasta.output.fai
    output: 
        txt = "genomes/{genome_build}/genome_fasta/main_chromosomes.txt",
        bed = "genomes/{genome_build}/genome_fasta/main_chromosomes.bed",
        patterns = temp("genomes/{genome_build}/genome_fasta/main_chromosomes.patterns.txt")
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
        chmod a-w {output.txt} {output.bed}
        """)


rule create_bwa_index:
    input: 
        fasta = rules.get_genome_fasta_download.output.fasta
    output: 
        prefix = touch("genomes/{genome_build}/bwa_index/bwa-{bwa_version}/genome.fa")
    log: 
        "genomes/{genome_build}/bwa_index/bwa-{bwa_version}/genome.fa.log"
    conda: CONDA_ENVS["bwa"]
    shell:
        op.as_one_line("""
        bwa index -p {output.prefix} {input.fasta} > {log} 2>&1
            &&
        find `dirname {output.prefix}` -type f -exec chmod a-w {{}} \;
        """)


rule get_gencode_download: 
    input:
        gtf = get_download_file(rules.download_gencode_annotation.output.gtf)
    output:
        gtf = "genomes/{genome_build}/annotations/gencode_annotation-{gencode_release}.gtf"
    shell:
        "ln -f {input.gtf} {output.gtf}"


rule create_star_index:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta,
        gtf = get_download_file(rules.download_gencode_annotation.output.gtf)
    output: 
        index = directory("genomes/{genome_build}/star_index/star-{star_version}/gencode-{gencode_release}/overhang-{star_overhang}")
    log: 
        "genomes/{genome_build}/star_index/star-{star_version}/gencode-{gencode_release}/overhang-{star_overhang}/log"
    conda: CONDA_ENVS["star"]
    threads: 12
    shell:
        op.as_one_line("""
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.index}
        --genomeFastaFiles {input.fasta} --sjdbOverhang {wildcards.star_overhang}
        --sjdbGTFfile {input.gtf} --outTmpDir {output.index}/_STARtmp > {log} 2>&1
        --outFileNamePrefix {output.index}/
            &&
        find {output.index} -type f -exec chmod a-w {{}} \;
        """)


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


rule create_salmon_index:
    input:
        fasta = rules.get_genome_fasta_download.output.fasta,
        gtf = get_download_file("downloads/gencode-33/gencode.annotation.{version}.gtf")
    output: 
        index = directory("genomes/{genome_build}/salmon_index/salmon-{salmon_version}")
    log: 
        "genomes/{genome_build}/salmon_index/salmon-{salmon_version}/log"
    conda: CONDA_ENVS["salmon"]
    threads: 8
    shell:
        op.as_one_line("""
        salmon index
        -t {input.fasta}
        -i {output.index}
        > {log} 2>&1
            &&
        chmod a-w {output.index}
        """)

