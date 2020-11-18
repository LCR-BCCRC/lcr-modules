#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import glob 
import string
import copy

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["starfish"]`
CFG = op.setup_module(
    name = "starfish",
    version = "2.0",
    subdirectories = ["inputs", "starfish", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _starfish_input_vcf,
    _starfish_rename_output, 
    _starfish_output_vcf,
    _starfish_output_union,
    _starfish_output_venn,
    _starfish_all,

##### GLOBAL VARIABLES #####

# All variant callers from config
callers = list(CFG["inputs"]["vcf"].keys())
callers = [caller.lower() for caller in callers]

# The name of the union VCF file
union_vcf = copy.deepcopy(callers)
union_vcf.append("union")
union_vcf = "_".join(union_vcf)

assert len(callers) < 6, (
    "Starfish can only handle a maximum of 6 input VCFs."
)

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _starfish_input_vcf:
    input:
        vcf = lambda w: config["lcr-modules"]["starfish"]["inputs"]["vcf"][w.caller] 
    output:
        vcf = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{caller}.vcf.gz", 
        tbi = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{caller}.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf), 
        op.relative_symlink(input.vcf + ".tbi", output.tbi)


# Run Starfish
rule _starfish_run:
    input:
        vcfs = expand(
            CFG["dirs"]["inputs"] + "vcf/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}.{caller}.vcf.gz", 
            caller = callers
            ),
        tbis = expand(
            CFG["dirs"]["inputs"] + "vcf/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}.{caller}.vcf.gz.tbi", 
            caller = callers
            ),
        reference = reference_files("genomes/{genome_build}/sdf"),
        starfish_script = CFG["inputs"]["starfish_script"]
    output:
        complete = touch(CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish.complete"), 
        venn = CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/venn.pdf"
    log:
        stdout = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish_run.stdout.log",
        stderr = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish_run.stderr.log"
    params:
        opts = CFG["options"]["starfish_run"]
    conda:
        CFG["conda_envs"]["starfish"]
    threads:
        CFG["threads"]["starfish_run"]
    resources:
        **CFG["resources"]["starfish_run"]
    shell:
        op.as_one_line("""
        {input.starfish_script} --sdf {input.reference} -O $(dirname {output.complete})
        --names {callers} --threads {threads}
        {params.opts} --venn --verbose
        -V {input.vcfs} > {log.stdout} 2> {log.stderr} 
        """)

checkpoint _starfish_rename_output: 
    input: 
        complete = CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish.complete"
    output: 
        renamed = touch(CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/starfish.renamed")
    run: 
        # Get the starfish letter names for the variant callers
        names = list(string.ascii_uppercase[0:(len(callers))])
        # Generate a dictionary of letter names and variant callers
        name_dict = {names[i]: callers[i] for i in range(len(callers))}
        outdir = os.path.dirname(input.complete)
        for src in glob.glob(outdir + "/[A-Z]*.vcf.gz*"): 
            bname = os.path.basename(src)
            lhs, rhs = bname.split(".", 1)
            lhs = "_".join(lhs)
            for key in list(name_dict.keys()): 
                lhs = lhs.replace(key, name_dict[key])
            dest = os.path.join(outdir, ".".join([lhs, rhs]))
            os.rename(src, dest) 


def _starfish_get_output(wildcards): 
    CFG = config["lcr-modules"]["starfish"]
    # Get the path to the output directory from the checkpoint _starfish_rename_output
    # **wildcards unpacks the wildcards object from the rule where this input function is called
    checkpoint_output = os.path.dirname(checkpoints._starfish_rename_output.get(**wildcards).output.renamed)
    # Obtain the vcf_name wildcards using glob_wildcards and expand that into the vcfs object
    vcfs = expand(
        CFG['dirs']['starfish'] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf.gz",
        **wildcards,  
        vcf_name = glob_wildcards(os.path.join(checkpoint_output, "{vcf_name, [a-z]\S+}.vcf.gz")).vcf_name # Constrain to omit 2, 2+, 2-, etc. VCF files. 
        )
    return vcfs

# Create a union VCF
rule _starfish_union: 
    input: 
        _starfish_get_output
    output: 
        union = CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{union_vcf}.vcf.gz", 
        union_tbi = CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{union_vcf}.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{union_vcf}.stdout.log",
        stderr = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{union_vcf}.stderr.log"
    params: 
        options = CFG["options"]["starfish_union"]
    conda: 
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["starfish_union"]
    resources:
        **CFG["resources"]["starfish_union"]
    wildcard_constraints: 
        union_vcf = union_vcf
    shell: 
        op.as_one_line("""
        bcftools merge {params.options} --force-samples 
        -Oz -o {output.union} 
        {input} && 
        tabix -p vcf {output.union}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _starfish_output_vcf:
    input:
        vcf = CFG['dirs']['starfish'] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf.gz"
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz", 
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz.tbi"
    wildcard_constraints: 
        vcf_name = "\S+(?<!union)"
    run:
        op.relative_symlink(input.vcf, output.vcf), 
        op.relative_symlink(input.vcf + ".tbi", output.tbi)

rule _starfish_output_union:
    input:
        vcf = CFG['dirs']['starfish'] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{union_vcf}.vcf.gz", 
        tbi = CFG['dirs']['starfish'] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{union_vcf}.vcf.gz.tbi",
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{union_vcf}.vcf.gz", 
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{union_vcf}.vcf.gz.tbi"
    wildcard_constraints: 
        union_vcf = union_vcf
    run:
        op.relative_symlink(input.vcf, output.vcf), 
        op.relative_symlink(input.tbi, output.tbi)

rule _starfish_output_venn: 
    input: 
        venn = str(rules._starfish_run.output.venn)
    output: 
        venn = CFG["dirs"]["outputs"] + "venn/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.venn.pdf"
    run: 
        op.relative_symlink(input.venn, output.venn)

def _starfish_get_output_target(wildcards): 
    CFG = config["lcr-modules"]["starfish"]
    # Get the path to the output directory from the checkpoint _starfish_rename_output
    # **wildcards unpacks the wildcards object from the rule where this input function is called
    checkpoint_output = os.path.dirname(checkpoints._starfish_rename_output.get(**wildcards).output.renamed)
    # Obtain the vcf_name wildcards using the glob_wildcards function
    vcfs = expand(
        CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf.gz",
        **wildcards,
        vcf_name = glob_wildcards(os.path.join(checkpoint_output, "{vcf_name}.vcf.gz")).vcf_name
        )
    return vcfs


rule _starfish_dispatch: 
    input: 
        expand(CFG["dirs"]["outputs"] + "vcf/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}.{union_vcf}.vcf.gz", union_vcf = union_vcf),
        _starfish_get_output_target        
    output: 
        touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _starfish_all:
    input:
        expand(
            [
                str(rules._starfish_output_venn.output.venn), 
                str(rules._starfish_dispatch.output)
            ],
            zip,  
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
