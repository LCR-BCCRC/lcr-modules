#!/usr/bin/env snakemake


localrules: _vcf2maf_input, _vcf2maf_output, _vcf2maf_all

##### RULES #####


rule _vcf2maf_run:
    input:
        vcf = CFG["dirs"]["vcf_maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf"
    output:
        maf = CFG["dirs"]["vcf_maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.maf"
    log:
        stdout = CFG["logs"]["vcf_maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf2maf.stdout.log",
        stderr = CFG["logs"]["vcf_maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf2maf.stderr.log"
    params:
        opts = "--filter-vcf 0", #CFG["options"]["vcf2maf"],
        fasta  = md.get_reference(CFG, "genome_fasta"),
        vep = md.get_reference(CFG, "vep_dir")
    conda:
        CFG["conda_envs"].get("vcf2maf") or "envs/vcf2maf-1.6.17.yaml"
    threads:
        CFG["threads"].get("vcf2maf") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("vcf2maf") or 4000
    shell:
        md.as_one_line("""
        vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} 
        --ref-fasta {params.fasta} 
        --vep-data {params.vep} 
        --vep-path $(dirname $(which vep))/../share/variant-effect-predictor* 
        {params.opts} 
        > {log.stdout} 2> {log.stderr}
        """)
