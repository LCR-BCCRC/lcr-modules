##### PYTHON MODULES #####
import os
import oncopipe as op

##### SETUP VARIABLES #####

# the contents of config/general.yaml will be loaded by header.smk
config["tool_names"] = "sv_repair"
config["pipeline_name"]= "sv_repair"

include: "header_2.1.smk"


configfile: f"config/sv_repair.yaml"

OUTPUT_DIR = "results/gambl/sv_repair/"

##### SPECIFY SAMPLES #####

tool_name = "sv_repair"
if tool_name in PIPELINE_SAMPLES["GENOMES"]:
    print("using custom set")
    config['lcr-modules'][tool_name]['samples'] = PIPELINE_SAMPLES["GENOMES"][tool_name]
else:
    config['lcr-modules'][tool_name]['samples'] = PIPELINE_SAMPLES["GENOMES"]['shared']

print(config['lcr-modules'][tool_name]['samples'])


##### RULES #####

rule input_bedpe:
	input:
		sample_bedpe =  "results/gambl/manta-2.3/04-bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.bedpe",
		tumour_bam = "data/genome_bams/{tumour_id}.{genome_build}.bam",
		tumour_bai = "data/genome_bams/{tumour_id}.{genome_build}.bam.bai"
	output:
		sample_bedpe = OUTPUT_DIR + "01-input_bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.bedpe",
		tumour_bam = OUTPUT_DIR + "01-input_bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.{genome_build}.bam",
		tumour_bai = OUTPUT_DIR + "01-input_bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.{genome_build}.bam.bai"
	run:
		op.absolute_symlink(input.sample_bedpe,output.sample_bedpe),
		op.absolute_symlink(input.tumour_bam,output.tumour_bam),
		op.absolute_symlink(input.tumour_bai,output.tumour_bai)

rule add_windows_bedpe:
	input:
		bedpe = str(rules.input_bedpe.output.sample_bedpe)
	output:
		bedpe = OUTPUT_DIR + "02-window_bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.bedpe"
	log:
		stdout = OUTPUT_DIR + "00-logs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/add_windows.stdout.log"
	params:
		window = 1000
	conda:
		"../../config/envs/sv_repair.yaml"
	script:
		"../../src/R/bedpe.windows.R"
		

rule bedpe_to_bed:
	input:
		bedpe = str(rules.add_windows_bedpe.output.bedpe)
	output:
		bed = OUTPUT_DIR + "03-bed/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.bed"
	log:
		stdout = OUTPUT_DIR + "00-logs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/bedpe.to.bed.stdout.log"
	conda:
		"../../config/envs/sv_repair.yaml"
	script:
		"../../src/R/bedpe.to.bed.R"

rule get_fasta:
	conda:
		"../../config/envs/bedtools.yaml"		
	input:
		bed = str(rules.bedpe_to_bed.output.bed),
		genome_fasta = "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/ref/lcr-modules-references-STABLE/genomes/{genome_build}/genome_fasta/genome.fa"
	log:
		stdout = OUTPUT_DIR + "00-logs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/get_fasta.stdout.log"
	output:
		reference_fasta = OUTPUT_DIR + "04-fasta/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/reference.fa"
	shell:
		op.as_one_line("""
			bedtools getfasta -s -fo {output.reference_fasta} -fi {input.genome_fasta} -bed {input.bed} 
		""")

# rule get_reads:
# 	conda:
		
# 	input:
		#genome_bam = str(rules.input_bedpe.output.tumour_bam),
# 	log:
		
# 	params:
		
# 	output:
 		
# 	shell:


# rule align_reads:
# 	conda:
		
# 	input:
		
# 	log:
		
# 	params:
		
# 	output:
 		
# 	shell:

rule all:
    input:
        expand(
            [
                str(rules.get_fasta.output.reference_fasta),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=["tumour_seq_type"],
            genome_build=["tumour_genome_build"],
            tumour_id=["tumour_sample_id"],
            normal_id=["normal_sample_id"],
            pair_status=["pair_status"])