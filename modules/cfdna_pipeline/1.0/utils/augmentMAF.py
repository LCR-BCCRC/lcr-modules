"""
This script is used to add variants complete with read depths found in additional mafs to an index maf.

The code looks for variants in using start position, chromosome, and alt allele that are found in a list of
additional maf files but not in the index maf file. If found, the script will update the index maf
file with the new variants with allele counts > 0 from the index bam file. If no new variants are found
it will return the index maf file as is.

The original script is augment_ssm by many authors. So here is a simplified version that is more
readable and maintainable. Also easier to incorporate into pipelines with the CLI, or import into your
own python scripts and use the class AugmentMAF. The core logic on what the script does is the same, just
in some prettier syntax and more documentation.

todo:
    I did not touch the "_get_read_support" function, didn't want to break it. But it could use a review
    and maybe some refactoring.

Author: Kurt Yakimovich et al.
"""

import pandas as pd
from pathlib import Path 
import pysam
import os
import multiprocessing
import argparse
import shutil

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_id',required=True,type=str,help='tumour sample id')
    parser.add_argument('--threads',required=False,type=int,default=24,help='number of threads to use, default is 24')
    parser.add_argument('--index_maf',required=True,type=str,help='index maf file')
    parser.add_argument('--index_bam',required=True,type=str,help='index bam file')
    parser.add_argument('--add_maf_files',required=False,nargs='*',type=str,help='additional maf files')
    parser.add_argument('--genome_build',required=True,type=str,help='genome build')
    parser.add_argument('--alt_count_min',required=True,type=int,help='minimum alt count reads needed to keep variant')
    parser.add_argument('--output',required=True,type=str,help='output file')
    return parser.parse_args()


class AugmentMAF(object):
    """Add variants complete with read depths found in additional mafs to an index maf.

    Useful if you have multiple samples from the same patient and want to see if there is supporting
    reads for a variant in one sample that is not called in another sample.

    Attributes:
        sample_id (str): sample id
        index_maf (pd.DataFrame): index maf file
        index_bam (str): index bam file
        add_maf_files (list): list of additional maf files
        genome_build (str): genome build
        output (str): output file
    
    Methods:
        read_maf: read a maf file
        get_genome_chromosomes: get a list of chromosomes for the genome build
        write_output: write the augmented maf to a file
        augment_maf: loop through each chromosome in parallel and check for variants in the index bam from
                    the add maf files. If found, update the index maf file with the new variant information.
        filter_augmented_vars: filter variants in the augmented maf based on the minimum alt depth and alt frequency
        _get_missing_maf_rows: compare other maf to index maf and return rows that are missing from the index maf
                            based on chromosome and start position of variants.
        _subset_and_run: subset variants per chr and find support in index sample bam file. 
        _get_read_support: get read support for a variant in a bam file.

    """

    def __init__(self, sample_id: str, index_maf:str,
                        index_bam: str, add_maf_files: list, genome_build: str,
                        output: str , threads: int = 6, min_alt_count: int= 3):
        super(AugmentMAF, self).__init__()
        # user provided inputs
        self.sample_id = sample_id
        self.threads = threads
        self.index_maf = self.read_maf(index_maf)
        self.index_bam = index_bam
        self.add_maf_files = add_maf_files
        self.genome_build = genome_build
        self.output = output
        self.min_alt_count = min_alt_count
        # computed variables
        self.chromosomes = self.get_genome_chromosomes()
        self.master_maf = self.run_or_not()
        self.augmented_vaf = self.augment_maf() # the magix

    def read_maf(self, maf_file: str) -> pd.DataFrame:
        return pd.read_csv(maf_file, sep="\t", comment='#',
        dtype={'Chromosome':'string','t_ref_count':'string','t_alt_count':'string','t_depth':'string', "Tumor_Sample_Barcode":"string"},
        low_memory=False)

    def get_genome_chromosomes(self)-> list:
        """Return a list of chromosomes for the genome build.
        """
        if self.genome_build.lower() in ["hg38","hg19-reddy", "grch38"]:
            return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        else:
            return [str(i) for i in range(1, 23)] + ["X", "Y"]

    def write_output(self, output_file: str = None) -> None:
        """Write the augmented maf to a file.

        Args:
            output_file (str): output file name
        """
        if output_file:
            out_name = output_file
        elif self.output:
            out_name = self.output
        else:
            raise print("No output file name provided ...")
        self.augmented_vaf.to_csv(out_name, sep="\t", index=False)

    def augment_maf(self):
        """Loop through each chromosome in parallel
        and check for variants in the index bam from
        the add maf files. If found, update the index
        maf file with the new variant information.


        """
        # test if master_maf is a df or not
        if self.master_maf is None:
            print("No additional maf files provided ... returning index maf")
            # add standard column anyways
            self.index_maf["origin_samples"] = ""
            return self.index_maf

        # if there are no new variants to add, return the index maf
        new_vars = self.master_maf.loc[self.master_maf["variant_source"] == "additional_maf"].shape[0]
        if new_vars == 0:
            print("No new variants to add to index maf ... returning index maf")
            # add standard column anyways
            self.index_maf["origin_samples"] = ""
            return self.index_maf
        elif new_vars > 0:
            print(f"Found {new_vars} new variants to add to index maf ...")

        # run in parallel for all chromosomes
        # run in parallel
        cool_pool = multiprocessing.Pool(processes = self.threads)
        results = cool_pool.starmap(self._subset_and_run, [[chrm] for chrm in self.chromosomes])
        cool_pool.close()
        cool_pool.join() # Keeping Chris' variable naming

        # concatenate mafs and write them out to a file
        augmented_maf_merged = pd.concat(results)

        # add autmented vaf to index maf, reset sample_id in all rows, and remove duplicates
        full_augmented_maf_merged = pd.concat([self.index_maf, augmented_maf_merged])
        full_augmented_maf_merged["Tumor_Sample_Barcode"] = self.sample_id
        full_augmented_maf_merged = full_augmented_maf_merged.drop_duplicates().reset_index(drop=True)

        # drop key column
        if "key" in full_augmented_maf_merged.columns:
            full_augmented_maf_merged.drop(columns=["key"], inplace=True)
        # recalculate AF
        # droop any rows with NA in Hugo_Symbol or t_alt_count, t_depth
        full_augmented_maf_merged = full_augmented_maf_merged.dropna(subset=["Hugo_Symbol","t_alt_count","t_depth"]).copy()
        full_augmented_maf_merged["AF"] = full_augmented_maf_merged["t_alt_count"].astype(int) / full_augmented_maf_merged["t_depth"].astype(int)
        # insert AF column at 43 column
        full_augmented_maf_merged.insert(43, "AF", full_augmented_maf_merged.pop("AF")) # make it easier to find goodness sake

        return full_augmented_maf_merged

    def _get_missing_maf_rows(self):
        """Compare other maf to index maf and 
        return rows that are missing from the index maf
        based on chromosome and start position of variants.

        """

        additional_mafs = []
        for maf in self.add_maf_files:
            other_maf = self.read_maf(maf)
            other_maf["variant_source"] = "additional_maf"
            additional_mafs.append(other_maf)
        # add column saying what samples vars came from
        additional_vars = self.mark_origin(pd.concat(additional_mafs))

        # concatenate all mafs with index maf and drop duplicate keys, but keep index maf version
        all_mafs = pd.concat([self.index_maf, additional_vars]).reset_index(drop=True).copy()
        # add variant key
        all_mafs["variant_key"] = all_mafs["Chromosome"] + "_" + all_mafs["Start_Position"].astype(str) + "_" + all_mafs["Tumor_Seq_Allele2"].astype(str)
        # keep only the first instance of a variant, so if its in the index maf, keep that one
        all_mafs = all_mafs.drop_duplicates(subset="variant_key", keep="first").copy()
        all_mafs.drop(columns=["variant_key"], inplace=True)

        return all_mafs

    def run_or_not(self):
        """Check if there are additional maf files to work with."""
        if len(self.add_maf_files) > 0:
            return self._get_missing_maf_rows()
        else:
            return None

    def mark_origin(self, add_mafs: pd.DataFrame) -> pd.DataFrame:
        """Mark the origin sample(s) of each variant
        if it is form an additional maf file. 

        Args:
            add_mafs (pd.DataFrame): additional maf files dataframe

        Returns:
            add_mafs (pd.DataFrame): additional maf files dataframe with origin samples column
        """
        # add variant key
        add_mafs["variant_key"] = add_mafs["Chromosome"] + "_" + add_mafs["Start_Position"].astype(str) + "_" + add_mafs["Tumor_Seq_Allele2"].astype(str)
        # create dict with variant_key as key and a comma separated list of samples it was found in as value
        variant_key_dict = add_mafs.groupby("variant_key")["Tumor_Sample_Barcode"].agg(lambda x: ",".join(list(x))).to_dict()
        # create a new column with the origin samples
        add_mafs["origin_samples"] = add_mafs["variant_key"].map(variant_key_dict)

        return add_mafs


    def _subset_and_run(self, this_chromosome: str) -> pd.DataFrame:
        """Subset variants per chr and find support in index sample bam file.

        Args:
            this_chromosome (str): chromosome to work on

        Returns:
            missing_maf (pd.DataFrame): maf with updated read counts for
                                        variants not in index but in additional
                                        maf files.
        """
        # subset the maf to only the chromosome we are working on
        print(f"Working on chromosome {this_chromosome} ...")
        missing_maf = self.master_maf[self.master_maf["Chromosome"] == this_chromosome].copy()
        # find only rows with the source "additional_maf" aka were not called in the index maf
        missing_maf = missing_maf[missing_maf["variant_source"] == "additional_maf"].copy()

        print(f"Found these variants to augment: {missing_maf}")
        # loop through variants in missing vaf and get support

        for var in missing_maf.itertuples():
            # get read support
            read_counts = self._get_read_support(chrom=var.Chromosome,
                            start= int(var.Start_Position),
                            end=int(var.End_Position),
                            ref_base=var.Reference_Allele,
                            alt_base=var.Tumor_Seq_Allele2,
                            bamfile=self.index_bam,
                            mut_class=var.Variant_Type)
            # take values from returned read_counts dict and update missing maf
            missing_maf.at[var.Index, "t_depth"] = read_counts["tumor_depth"]
            missing_maf.at[var.Index, "t_ref_count"] = read_counts["tumor_ref_count"]
            missing_maf.at[var.Index, "t_alt_count"] = read_counts["tumor_alt_count"]

        print(f"Finished working on chromosome {this_chromosome} ...")
        # only return variants with a t_alt_count > min_alt_count
        missing_maf = missing_maf[(missing_maf["t_alt_count"].astype(int) > self.min_alt_count) & (missing_maf["t_depth"].astype(int) > self.min_alt_count)].copy()
        return missing_maf


    def _get_read_support(self, chrom: str, start: int, end: int, ref_base: str, alt_base: str, 
                        bamfile: str, mut_class: str, min_base_qual=10, min_mapping_qual=1,
                        padding: int = 100) -> dict:
        """Get read support for a variant in a bam file.

        generate a pileup and get the ref and alt read count for each variant (row) in the missing
        rows add some padding to each region.

        Deletions are handled in a somewhat loosey-goosey way. Basically, the count is the maximum number 
        of reads supporting a deletion at any site within the region start-end. This is because there were
        some examples of a few less reads supporting the deletion near the end (likely the gap in the alignment
        was shifted). DNPs are handled in a super lazy way that should be fine most of the time.
        The allele is only checked for matching at the first base in the DNP.


        Args:
            chrom (str): chromosome of variant
            start (int): start position of variant
            end (int): end position of variant
            ref_base (str): reference allele
            alt_base (str): alt allele
            bamfile (str): path to bam file
            mut_class (str): mutation class (SNP, DNP, INS, DEL)
            min_base_qual (int): minimum base quality
            min_mapping_qual (int): minimum mapping quality
            padding (int): padding around the start location of the variant.

        Returns:
            dict: dictionary with read support values, tumor_depth, tumor_ref_count, tumor_alt_count.
        """

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
        for pileupcolumn in samfile.pileup(chrom, start, end, ignore_overlaps=False, 
                                                            min_base_quality=min_base_qual,
                                                            min_mapping_quality=min_mapping_qual,
                                                            truncate=True):
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
        return {'tumor_depth':str(total_depth),'tumor_ref_count':str(ref_count),'tumor_alt_count':str(alt_count)}

def main():
    args = get_args()
    # if no additional maf files are provided, print warning and do nothing
    # initialize the class
    augment_maf = AugmentMAF(args.sample_id, args.index_maf, args.index_bam, 
                            args.add_maf_files, args.genome_build, args.output, 
                            threads=args.threads, min_alt_count=args.alt_count_min)
    # write the output
    augment_maf.write_output()

if __name__ == '__main__':
    main()
