#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton and Shaghayegh Soudi
# Module Author:    Laura Hilton and Shaghayegh Soudi
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import hashlib


# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["phylowgs"]`
CFG = op.setup_module(
	name = "phylowgs",
	version = "1.0",
	subdirectories = ["inputs", "preprocess_battenberg", "preprocess_inputs", "multievolve", "results",  "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
	_phylowgs_input_vcf,
	_phylowgs_input_battenberg,
	_phylowgs_output_html,
	_phylowgs_all


# Generate a de-duplicated table of patient_ids etc. 
PATIENTS = CFG["runs"][["tumour_patient_id", "normal_patient_id", "tumour_genome_build", "tumour_seq_type", "tumour_sex"]].drop_duplicates(subset = None, ignore_index = True)

# Obtain the path to the phylowgs conda environment
md5hash = hashlib.md5()
if workflow.conda_prefix: 
	conda_prefix = workflow.conda_prefix
else: 
	conda_prefix = os.path.abspath(".snakemake/conda")
md5hash.update(conda_prefix.encode())
f = open(CFG['conda_envs']['phylowgs'], 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()
PHYLO = conda_prefix + "/" + h[:8] + "/share/phylowgs/"


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _phylowgs_input_vcf:
	input:
		vcf = CFG["inputs"]["vcf"],
		tbi = CFG["inputs"]["tbi"]
	output:
		vcf = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched.vcf.gz",
		tbi = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched.vcf.gz.tbi"
	run:
		op.absolute_symlink(input.vcf, output.vcf)
		op.absolute_symlink(input.tbi, output.tbi)


rule _phylowgs_input_battenberg:
	input:
		cellularity = CFG["inputs"]["cellularity"],
		subclones = CFG["inputs"]["subclones"]
	output:
		cellularity = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched.cellularity_ploidy.txt",
		subclones = CFG["dirs"]["inputs"] + "battenberg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched.subclones.txt"
	run:
		op.absolute_symlink(input.cellularity, output.cellularity)
		op.absolute_symlink(input.subclones, output.subclones)


# Preprocess the battenberg file to match requirements
rule _phylowgs_preprocess_battenberg: 
	input:
		cellularity = str(rules._phylowgs_input_battenberg.output.cellularity),
		subclones = str(rules._phylowgs_input_battenberg.output.subclones)
	output:
		txt = CFG["dirs"]["preprocess_battenberg"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--matched.cnvs.txt"
	log:
		stderr = CFG["logs"]["preprocess_battenberg"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--matched.preprocess_battenberg.stderr.log", 
		stdout = CFG["logs"]["preprocess_battenberg"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--matched.preprocess_battenberg.stdout.log"
	params:
		script = PHYLO + "parser/parse_cnvs.py"
	conda: 
		 CFG["conda_envs"]["phylowgs"]
	threads: 
		CFG["threads"]["create_inputs"]
	resources: 
		**CFG["resources"]["create_inputs"]
	shell:
		op.as_one_line("""
		cellularity=$(tail -n +2 {input.cellularity} | cut -f 1); 
		python2 {params.script} -f battenberg-smchet -c $cellularity --cnv-output {output.txt} {input.subclones} 
		2> {log.stderr} > {log.stdout}
		""")

# Expand the input files to create a command-line argument for create_phylowgs_inputs.py
def create_phylowgs_inputs_cli(wildcards): 
	CFG = config["lcr-modules"]["phylowgs"]
	PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id)
	cnvs = expand(
		"--cnvs {time_point}=" + CFG['dirs']['preprocess_battenberg'] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--matched.cnvs.txt", 
		zip, 
		time_point = PATIENT["tumour_time_point"], 
		tumour_id = PATIENT["tumour_sample_id"], 
		normal_id = PATIENT["normal_sample_id"], 
		seq_type = PATIENT["tumour_seq_type"], 
		genome_build = PATIENT["tumour_genome_build"], 
		patient_id = PATIENT["tumour_patient_id"]
	)
	vcf_types = expand(
		"--vcf-type {time_point}=" + CFG['options']['create_inputs']['vcf_type'], 
		zip, 
		time_point = PATIENT["tumour_time_point"]
	)
	vcfs = expand(
		"{time_point}=" + CFG['dirs']['inputs'] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched.vcf.gz", 
		zip, 
		time_point = PATIENT["tumour_time_point"], 
		tumour_id = PATIENT["tumour_sample_id"], 
		normal_id = PATIENT["normal_sample_id"], 
		seq_type = PATIENT["tumour_seq_type"], 
		genome_build = PATIENT["tumour_genome_build"], 
		patient_id = PATIENT["tumour_patient_id"]
	)
	cli = cnvs + vcf_types + vcfs
	cli = " ".join([str(elem) for elem in cli])
	return(cli)

# Input function to pull in input VCF and preprocessed CNV data
def create_phylowgs_inputs(wildcards): 
	CFG = config["lcr-modules"]["phylowgs"]
	PATIENT = op.filter_samples(CFG["runs"], tumour_patient_id = wildcards.patient_id)
	inputs = expand(
		[
			CFG["dirs"]["preprocess_battenberg"] + "{seq_type}--{genome_build}/{patient_id}/{tumour_id}--{normal_id}--matched.cnvs.txt",
			CFG['dirs']['inputs'] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched.vcf.gz"
		], 
		zip, 
		tumour_id = PATIENT["tumour_sample_id"], 
		normal_id = PATIENT["normal_sample_id"], 
		allow_missing = True
	)
	return(inputs)


# Preprocess vcf and battenberg inputs together
rule _phylowgs_create_inputs:
	input:
		create_phylowgs_inputs
	output:
		ssms = CFG["dirs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/ssm_data.txt", 
		cnvs = CFG["dirs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/cnv_data.txt", 
		params = CFG["dirs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/params.json"
	log:
		stderr = CFG["logs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/create_phylowgs_inputs.stderr.log",
		stdout = CFG["logs"]["preprocess_inputs"] + "{seq_type}--{genome_build}/{patient_id}/create_phylowgs_inputs.stdout.log"
	params:
		cli = create_phylowgs_inputs_cli,
		opts = CFG["options"]["create_inputs"]["opts"], 
		sex = lambda w: config["lcr-modules"]["phylowgs"]["switches"]["sex"][op.filter_samples(PATIENTS, tumour_patient_id = w.patient_id)["tumour_sex"].values[0]] if op.filter_samples(PATIENTS, tumour_patient_id = w.patient_id)["tumour_sex"].values[0] in config["lcr-modules"]["phylowgs"]["switches"]["sex"].keys() else "auto",
		script = PHYLO + "parser/create_phylowgs_inputs.py"
	conda:
		CFG["conda_envs"]["phylowgs"]
	threads:
		CFG["threads"]["create_inputs"]
	resources:
		**CFG["resources"]["create_inputs"]
	shell:
		op.as_one_line("""
		python2 {params.script}
		--output-cnvs {output.cnvs} 
		--output-variants {output.ssms} 
		--output-params {output.params} 
		--sex {params.sex}
		{params.opts}
		{params.cli}
		2> {log.stderr} > {log.stdout}
		""")

# Run multievolve to sample trees and reconstruct phylogeny
rule _phylowgs_multievolve: 
	input: 
		**rules._phylowgs_create_inputs.output
	output: 
		trees = CFG["dirs"]["multievolve"] + "{seq_type}--{genome_build}/{patient_id}/trees.zip"
	log:
		stderr = CFG["logs"]["multievolve"] + "{seq_type}--{genome_build}/{patient_id}/multievolve.stderr.log",
		stdout = CFG["logs"]["multievolve"] + "{seq_type}--{genome_build}/{patient_id}/multievolve.stdout.log"
	params: 
		script = PHYLO + "multievolve.py",
		opts = CFG["options"]["multievolve"]
	conda:
		CFG["conda_envs"]["phylowgs"]
	threads: 
		CFG["threads"]["multievolve"]
	resources: 
		**CFG["resources"]["multievolve"]
	shell: 
		op.as_one_line("""
		python2 {params.script} 
		{params.opts}
		-n {threads} 
		-O $(dirname {output.trees}) 
		--ssms {input.ssms} 
		--cnvs {input.cnvs} 
		2> {log.stderr} > {log.stdout}
		""")

# Write the results	
rule _phylowgs_write_results: 
	input: 
		str(rules._phylowgs_multievolve.output.trees)
	output: 
		muts = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.muts.json",
		summ = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.summ.json",
		mutass = CFG["dirs"]["results"] + "{seq_type}--{genome_build}/{patient_id}/{patient_id}.mutass.zip"
	params: 
		script = PHYLO + "write_results.py", 
		opts = CFG["options"]["write_results"]
	conda: 
		CFG["conda_envs"]["phylowgs"]
	threads: 
		CFG["threads"]["write_results"]
	resources: 
		**CFG["resources"]["write_results"]
	shell: 
		op.as_one_line("""
		python2 {params.script} 
		{params.opts}
		{wildcards.patient_id}
		{input}
		{output.summ}.gz
		{output.muts}.gz
		{output.mutass} && 
		gunzip -f $(dirname {output.mutass})/*.gz 
		""")
	
# Symlinks the final output files to the witness directory in preparation for HTTP browsing
rule _phylowgs_output_html:
	input:
		mutass = str(rules._phylowgs_write_results.output.mutass)
	output:
		complete = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{patient_id}.phylo_complete"
	params: 
		witness = PHYLO + "witness/data/{patient_id}"
	run:
		op.absolute_symlink(os.path.split(input.mutass)[0], params.witness)
		f = open(output.complete, "a")
		f.write("To view PhyloWGS results, navigate to " + PHYLO + "witness\n")
		f.write("Run the following commands: \n")
		f.write("python2 index_data.py\n")
		f.write("python2 -m SimpleHTTPServer\n")
		f.write("On a local machine you will be able to view your results in a browser at http://localhost:8000\n")
		f.write("For a remote machine, launch the following command in a terminal: \n")
		f.write("ssh -N -L localhost:8000:localhost:8000 <ssh_config>\n")
		f.write("Now you can view your results in a browser at http://localhost:8000\n")
		f.close()


# Generates the target sentinels for each run, which generate the symlinks
rule _phylowgs_all:
	input:
		expand(
			[
				str(rules._phylowgs_output_html.output.complete),
			],
			zip,  # Run expand() with zip(), not product()
			seq_type=PATIENTS["tumour_seq_type"],
			genome_build=PATIENTS["tumour_genome_build"],
			patient_id=PATIENTS["tumour_patient_id"]
		)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
