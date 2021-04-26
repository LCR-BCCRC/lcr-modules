#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["pti"]`
CFG = op.setup_module(
    name = "pti",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "pti", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _pti_input_maf,
    _pti_step_2,
    _pti_output_txt,
    _pti_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _pti_input_maf:
    input:
        maf = CFG["inputs"]["sample_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
    run:
        op.relative_symlink(input.maf, output.maf)


# Create the required input file
rule _pti_generate_input:
    input:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
    output:
        tsv = CFG["dirs"]["pti"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--{pair_status}.pti_input.tsv"
    log:
        stderr = CFG["logs"]["pti"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/_pti_generate_input.stderr.log"
    params:
        **CFG["options"]["generate_input"]
    conda:
        CFG["conda_envs"]["pti"]
    threads:
        CFG["threads"]["generate_input"]
    resources:
        mem_mb = CFG["mem_mb"]["generate_input"]
    run:
        # Load the MAF file as a pandas data frame
        MAF = pd.read_csv(input.maf, sep = "\t", header = {params.header_row}, low_memory=False)
        # Filter the MAF file to include only coding/UTR/splicing variants in known genes
        MAF = op.discard_samples(MAF, Hugo_Symbol = "Unknown")
        MAF = op.filter_samples(MAF, Variant_Classification = ["3'UTR", "5'UTR", "Frame_Shift_Del", "Frame_Shift_Ins", "Intron", "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Region", "Splice_Site"])
        # Generate lists of data for each column in final table and store in a dictionary
        table_data = {'uniq_mutation_id': [], 'var_count': [], 'ref_count': [], 'gene': []}
        table_data['uniq_mutation_id'] = (MAF['Chromosome'].astype(str) + "_" + MAF['Start_Position'].astype(str) + "_" + MAF['End_Position'].astype(str)).tolist()
        table_data['var_count'] = MAF['t_alt_count'].tolist()
        table_data['ref_count'] = MAF['t_ref_count'].tolist()
        table_data['gene'] = MAF['Hugo_Symbol'].tolist()
        # Assemble table
        TABLE = pd.DataFrame.from_dict(table_data)
        # Write to file
        TABLE.to_csv(output.tsv, sep = "\t", header = True, index = False)


# Retrieve all time point maf files per patient
def get_maf_per_patient(wildcards): 
    PATIENT = op.filter_samples(CFG['runs'], tumour_patient_id = wildcards.patient_id)
    mafs = expand(
        str(rules._pti_input_maf.output.maf), 
        zip, 
        tumour_id = PATIENT['tumour_sample_id'], 
        normal_id = PATIENT['normal_sample_id'], 
        pair_status = PATIENT['pair_status'], 
        allow_missing = TRUE
    )
    return(mafs)

      
# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _pti_run:
    input:
        txt = str(rules._pti_step_1.output.txt)
    output:
        txt = CFG["dirs"]["pti"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.txt"
    log:
        stderr = CFG["logs"]["pti"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.txt} > {output.txt} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _pti_output_txt:
    input:
        txt = str(rules._pti_step_2.output.txt)
    output:
        txt = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.txt"
    run:
        op.relative_symlink(input.txt, output.txt)


# Generates the target sentinels for each run, which generate the symlinks
rule _pti_all:
    input:
        expand(
            [
                str(rules._pti_output_txt.output.txt),
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
