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
    subdirectories = ["inputs", "generate_tsv", "pti", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _pti_input_maf,
    _pti_generate_input,
    _pti_output_txt,
    _pti_all,

# Ensure multiple tumour samples per patient

RUNS = CFG['runs']
PATIENTS = pd.DataFrame(RUNS['tumour_patient_id'].value_counts())

DROPPED_PATIENTS = PATIENTS.loc[PATIENTS.tumour_patient_id == 1].index.tolist()
if len(DROPPED_PATIENTS) >= 1: 
    print("The following patients were dropped because they have only a single time point: \n")
    print(DROPPED_PATIENTS)

PATIENTS = PATIENTS.loc[PATIENTS.tumour_patient_id > 1].index.tolist()
assert len(PATIENTS) > 0, (
    "There are 0 patients with more than one time point. \n"
    "PTI can only be run with multiple time points per patient. \n"
)

# Subset runs table to only multi-timepoint samples
RUNS = op.filter_samples(RUNS, tumour_patient_id = PATIENTS)
CFG['runs'] = RUNS

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _pti_input_maf:
    input:
        maf = CFG["inputs"]["sample_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
    run:
        op.relative_symlink(input.maf, output.maf)

def get_maf_per_timepoint(wildcards): 
    PATIENT = op.filter_samples(config["lcr-modules"]["pti"]["runs"], tumour_time_point = wildcards.time_point, tumour_patient_id = wildcards.patient_id)
    maf = expand(
        str(rules._pti_input_maf.output.maf), 
        zip, 
        tumour_id = PATIENT['tumour_sample_id'], 
        normal_id = PATIENT['normal_sample_id'], 
        pair_status = PATIENT['pair_status'], 
        allow_missing = True
    )
    return(maf)


# Create the required input file
rule _pti_generate_input:
    input:
        get_maf_per_timepoint
    output:
        tsv = CFG["dirs"]["generate_tsv"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}_{time_point}.pti_input.tsv"
    log:
        stderr = CFG["logs"]["generate_tsv"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.{time_point}._pti_generate_input.stderr.log"
    params:
        **CFG["options"]["generate_input"]
    threads:
        CFG["threads"]["generate_input"]
    resources:
        mem_mb = CFG["mem_mb"]["generate_input"]
    run:
        # Load the MAF file as a pandas data frame
        MAF = pd.read_csv(input[0], sep = "\t", header = int(params.header_row), low_memory=False)
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
        TABLE.to_csv(str(output.tsv), sep = "\t", header = True, index = False)


# Retrieve all time point tsv files per patient
def get_tsv_per_patient(wildcards): 
    PATIENT = op.filter_samples(config["lcr-modules"]["pti"]["runs"], tumour_patient_id = wildcards.patient_id)
    tsvs = expand(
        str(rules._pti_generate_input.output.tsv), 
        zip, 
        time_point = PATIENT['tumour_time_point'], 
        allow_missing = True
    )
    return(tsvs)

      
# Example variant filtering rule (single-threaded; can be run on cluster head node)
rule _pti_run:
    input:
        get_tsv_per_patient
    output:
        complete = touch(CFG["dirs"]["pti"] + "{seq_type}--{genome_build}/{patient_id}/pti.complete")
    log:
        stderr = CFG["logs"]["pti"] + "{seq_type}--{genome_build}/{patient_id}/pti_run.stderr.log"
    conda: 
        CFG["conda_envs"]["pti"]
    params:
        **CFG["options"]["run_pti"]
    shell:
        op.as_one_line("""
        python {params.script} --AF {params.AF} --drivers {params.drivers} -i $(dirname {input[0]}) -o $(dirname {output.complete})
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
# rule _pti_output_txt:
#     input:
#         txt = str(rules._pti_step_2.output.txt)
#     output:
#         txt = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.txt"
#     run:
#         op.relative_symlink(input.txt, output.txt)

PATIENT_TABLE = CFG["runs"].drop_duplicates(subset = "tumour_patient_id", ignore_index = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _pti_all:
    input:
        expand(
            [
                str(rules._pti_run.output.complete),
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=PATIENT_TABLE["tumour_seq_type"],
            genome_build=PATIENT_TABLE["tumour_genome_build"],
            patient_id=PATIENT_TABLE["tumour_patient_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
