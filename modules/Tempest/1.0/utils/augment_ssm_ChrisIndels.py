#!/usr/bin/env python
# Inputs (via Snakemake):
# 1. one tumour bam (or cram) file or both a tumour and normal bam file
# 2. an index MAF file that will be updated with missing rows.
# The index MAF should have variants called from analysis of the bam file(s) supplied
# 3. One ore more additional MAFs with variants called from the same patient using other bam files
# Output: a new MAF file that contains all variants from the index MAF plus rows for every variant unique to the additional MAFs
# Important: this script will determine the tumour_ref_count and tumour_alt_count for new MAF rows using the supplied bam
# It uses the pysam pileup engine to calculate the VAF of SNVs, indels based on the bam (i.e. trusting the alignment as-is)

#If you are using tabular data that isn't a MAF, the minimum set of required columns are:
#Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Type, t_ref_count, t_alt_count, t_depth
# e.g. chr2 1234 1235 A T SNP 0 0 0
import argparse
from pathlib import Path
import pysam
import os
import pandas as pd
import multiprocessing
from multiprocessing import set_start_method
verbose = False #turn off print statements
min_map_qual = 1 #minimum mapping quality of a read to use it
padding = 100

#restrict processing to only one chromosome, if desired. For parallelization.
#this list should be automatically pulled from the MAFs but I worry that if one or more chromosomes is missing from the index MAF this might lead to issues

#normal_bam = ""


#TODO: Normal bam is not implemented yet! Currently this just extracts read counts for the mutation site in the tumour bam
#augmented_maf = "./00-14595_tumorC.grch37.augmented.maf"
#tumour_bam_file = "data/genome_bams/00-14595_tumorC.grch37.bam"
#normal_bam_file = "data/genome_bams/00-14595_normal.grch37.bam"
#index_maf_file = "results/gambl/slms-3_vcf2maf_current/99-outputs/genome--grch37/00-14595_tumorC--00-14595_normal--matched_slms-3.final.maf"

#add_maf_files = ["results/gambl/slms-3_vcf2maf_current/99-outputs/genome--grch37/00-14595_tumorA--00-14595_normal--matched_slms-3.final.maf",
#                "results/gambl/slms-3_vcf2maf_current/99-outputs/genome--grch37/00-14595_tumorB--00-14595_normal--matched_slms-3.final.maf"]
#remove the above lines and uncomment the next few lines to get this snakemake-ready


#special case: if index_maf_file is among add_maf files, zero all t_ref_count and t_alt_count values

def get_genome_chromosomes(genome_build: str)-> list:
    """Return a list of chromosomes for the genome build.
    """
    if genome_build.lower() in ["hg38","hg19-reddy", "grch38"]:
        return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    else:
        return [str(i) for i in range(1, 23)] + ["X", "Y"]

#normal_bam_file = snakemake.input["normal_bam"]
# I hope this will work but I've never tried it with >1 input file so this may need to be tweaked

def get_args():

    parser = argparse.ArgumentParser(description="Annotates allele counts")
    parser.add_argument("-m", "--maf", required=True, help="Input MAF file")
    parser.add_argument("-b", "--bam", metavar="BAM", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", metavar="MAF", required=True, help="Output MAF file")
    parser.add_argument("-g", "--genome", metavar="GENOME", required=True, help="Genome build")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use")

    return parser.parse_args()



# Input: two data frames from MAF files, an index maf and a supplemental MAF
# Output: just the rows from the supplemental MAF that are not present in the index MAF
def get_missing_maf_rows(this_index_maf,other_maf,this_chromosome):
    #subset on chromosome
    this_index_maf = this_index_maf[this_index_maf.Chromosome == this_chromosome]
    other_maf = other_maf[other_maf.Chromosome == this_chromosome]
    missing_from_index_maf = other_maf[other_maf.Start_Position.isin(this_index_maf.Start_Position) == False]
    return(this_index_maf)


#generate a pileup and get the ref and alt read count for each variant (row) in the missing rows
# add some padding to each region
# Deletions are handled in a somewhat loosey-goosey way. Basically, the count is the maximum number of reads supporting a deletion at any site
# within the region start-end. This is because there were some examples of a few less reads supporting the deletion near the end (likely the gap in the alignment was shifted)
# DNPs are handled in a super lazy way that should be fine most of the time. The allele is only checked for matching at the first base in the DNP
def get_read_support(chrom,start,end,ref_base,alt_base,bamfile,mut_class,min_base_qual=10, min_mapping_qual=1):
    original_start = start
    actual_start = start - 1 # because of zero vs 1-based indexing
    actual_end = end
    del_positions=[]
    if mut_class == "DEL":

        del_length = len(ref_base)
        del_positions = list(range(actual_start,end))
    if mut_class == "DNP":
        alt_base = alt_base[0]
        ref_base = ref_base[0]
    start = int(start) - padding
    end = int(end) + padding
    chrom = str(chrom)
    ref_count = 0
    alt_count = 0
    total_depth = 0
    #handle the case in which the bam file is a cram file. This information must be provided by the calling function
    #realname = os.readlink(bamfile)
    realname = str(Path(bamfile).resolve())
    #print(realname)
    if realname.endswith("cram"):
        samfile = pysam.AlignmentFile(realname,'rc')
    else:
        samfile = pysam.AlignmentFile(bamfile,'rb')
    read_bases = {}
    del_support = {}
    max_depth = 0

    for pos in del_positions:
        del_support[pos]=0
    #print(f"{chrom} {start} {end}")
    for pileupcolumn in samfile.pileup(chrom, start, end, ignore_overlaps=True, min_base_quality=min_base_qual,min_mapping_quality=min_mapping_qual,truncate=True):
        if mut_class == "INS":
            #I think this is double-counting overlapping reads but at least it does it consistently
            if pileupcolumn.pos == actual_start:
                seqs = pileupcolumn.get_query_sequences(add_indels=True)
                for base in seqs:
                    if "+" in base:
                        alt_count+=1
                    else:
                        ref_count+=1
                    total_depth+=1
        elif mut_class == "DEL":
            if pileupcolumn.pos in del_positions:
                seqs = pileupcolumn.get_query_sequences(add_indels=True)
                ref_count = 0
                if len(seqs)>max_depth:
                    max_depth = len(seqs)
                for base in seqs:
                    if base == "*":
                        del_support[pileupcolumn.pos]+=1
                    else:
                        ref_count+=1

        elif pileupcolumn.pos == actual_start:
            #print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
            for pileupread in pileupcolumn.pileups:
                if mut_class == "SNP" or mut_class == "DNP":
                    if not pileupread.is_del and not pileupread.is_refskip:
                        # query position is None if is_del or is_refskip is set.

                        thisname = pileupread.alignment.query_name
                        this_base = pileupread.alignment.query_sequence[pileupread.query_position]
                        this_qual = pileupread.alignment.query_qualities[pileupread.query_position]
                        if thisname in read_bases:
                            read_bases[thisname].append(this_base)
                        else:
                            read_bases[thisname] = [this_base]
    #now count up bases and deal with disagreement between read 1 and 2 from the same fragment
    if mut_class == "DEL":
        for position in del_support:
            if del_support[position] > alt_count:
                alt_count = del_support[position]
        ref_count = max_depth - alt_count
        total_depth = max_depth
    if mut_class == "SNP" or mut_class == "DNP":
        for readname in read_bases:
            bases = read_bases[readname]
            #at most there should be two from the same read (fragment/pair etc)
            if len(bases)==1:
                total_depth+=1
                if bases[0] == ref_base:
                    ref_count+=1
                elif bases[0] == alt_base:
                    alt_count+=1
            elif bases[0] == bases[1]:
                total_depth+=1
                if bases[0] == ref_base:
                    ref_count+=1
                elif bases[0] == alt_base:
                    alt_count+=1
            else:
                total_depth+=1
                if bases[0] == ref_base or bases[1] == ref_base:
                    ref_count+=1
                elif bases[0] == alt_base and bases[1] == alt_base:
                    alt_count+=1


    #return as strings to prevent pandas from turning them into floats and breaking the MAF in the process
    return({'tumor_depth':str(total_depth),'tumor_ref_count':str(ref_count),'tumor_alt_count':str(alt_count)})

#TODO/WARNING: check if the bam file is a cram file and supply this information to the pileup function. This might break if a cram is provided currently.
# there is code that does this in augment_sv.py but the implementation is more complex here so I haven't gotten to it
def subset_and_run(full_index_maf, bamfile, this_chromosome, other_maf_files):
    missing_maf_rows = []

    other_maf = pd.read_csv(other_maf_files,sep="\t", comment='#',dtype={'Chromosome':'string','t_ref_count':'string','t_alt_count':'string','t_depth':'string'},low_memory=False)
    to_process_maf = get_missing_maf_rows(full_index_maf,other_maf,this_chromosome)

    for existing_rows in missing_maf_rows:
        current_dim = to_process_maf.shape[0]
        this_dim = existing_rows.shape[0]
        to_process_maf = get_missing_maf_rows(existing_rows,to_process_maf,this_chromosome)
        #this should remove any overlaping rows
        new_dim = to_process_maf.shape[0]
        if verbose:
            print(f"subtracting from a maf with {current_dim} using {this_dim} now at {new_dim}")
    missing_maf_rows.append(to_process_maf)

    missing_maf = pd.concat(missing_maf_rows)
    total_to_process = missing_maf.shape[0]
    if verbose:
        print(f"will work on {total_to_process} new MAF rows using {bamfile}")
    chroms = missing_maf.Chromosome.tolist()
    start_positions = missing_maf.Start_Position.tolist()
    end_positions = missing_maf.End_Position.tolist()
    ref_alleles = missing_maf.Reference_Allele.tolist()
    alt_alleles = missing_maf.Tumor_Seq_Allele2.tolist()
    mutation_class = missing_maf.Variant_Type.tolist()

    ref_counts = []
    alt_counts = []
    depths = []
    l = len(alt_alleles)
    for i,start_position in enumerate(start_positions):
        read_counts = get_read_support(chrom=chroms[i],
                        start=start_position,
                        end=end_positions[i],
                        ref_base=ref_alleles[i],
                        alt_base=alt_alleles[i],
                        bamfile=bamfile,
                        mut_class=mutation_class[i])

        ref_counts.append(read_counts['tumor_ref_count'])
        alt_counts.append(read_counts['tumor_alt_count'])
        depths.append(read_counts['tumor_depth'])
        if not i % 500 and verbose:
            if i > 0:
                percent_complete = int(100 * i / l)
                print(f"{i} of {l} rows ({percent_complete}%) complete (chromosome in this thread is {this_chromosome})")
    #add new values to MAF, overwriting what is there
    missing_maf.t_depth = depths
    missing_maf.t_ref_count = ref_counts
    missing_maf.t_alt_count = alt_counts
    return(missing_maf)

def main(args=None):

    if args is None:
        args = get_args()
    # load index MAF into a data frame using pandas
    index_maf = pd.read_csv(args.maf ,sep="\t",   comment='#',dtype={'Chromosome':'string','t_ref_count':'string','t_alt_count':'string','t_depth':'string'},low_memory=False)

    #run in parallel for all chromosomes
    chromosomes = get_genome_chromosomes(args.genome)
    pool_args = []
    for chrom in chromosomes:
        pool_args.append([index_maf, args.bam, chrom, args.maf])

    cool_pool = multiprocessing.Pool(processes = args.threads)
    results = cool_pool.starmap(subset_and_run, pool_args)
    cool_pool.close()
    cool_pool.join() #Keeping Chris' variable naming

    full_augmented_maf_merged = pd.concat(results)
    #final step: update the Tumor_Sample_Barcode to refer to the index sample

    full_augmented_maf_merged.to_csv(args.output,sep="\t",index=False)

main()