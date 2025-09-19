"""
This script is used to add variants complete with read depths found in additional mafs to an index (original) maf.

The code looks for variants in using start position, chromosome, and alt allele that are found in a list of
additional maf files but not in the index maf file. If found, the script will update the index maf
file with the new variants with allele counts > 0 from the index bam file. If no new variants are found
it will return the index maf file as is.

This script is a derivitive of the "augment_ssm" by the Morin Lab. 
Key features added here are CLI support, and packaing of the code into a class that
has been refactored for readability and maintainability.

The core logic on what the script does is the same, just
in some prettier syntax and more documentation.

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
        """Main coordinator function for getting read support for a variant."""
        
        # Setup common parameters
        setup_params = self._setup_variant_params(start, end, ref_base, mut_class, padding)
        
        # Open BAM file
        samfile = self._open_bam_file(bamfile)
        
        # Process variant based on type
        if mut_class == "INS":
            ref_count, alt_count, total_depth = self._process_insertion(
                samfile, chrom, setup_params, min_base_qual, min_mapping_qual)
        elif mut_class == "DEL":
            ref_count, alt_count, total_depth = self._process_deletion(
                samfile, chrom, setup_params, ref_base, min_base_qual, min_mapping_qual)
        elif mut_class in ["SNP", "DNP"]:
            ref_count, alt_count, total_depth = self._process_snp_dnp(
                samfile, chrom, setup_params, ref_base, alt_base, mut_class, 
                min_base_qual, min_mapping_qual)
        else:
            raise ValueError(f"Unsupported mutation class: {mut_class}")
        
        samfile.close()
        
        return self._format_results(ref_count, alt_count, total_depth)

    def _setup_variant_params(self, start: int, end: int, ref_base: str, mut_class: str, padding: int) -> dict:
        """Setup common parameters for variant processing."""
        actual_start = start - 1  # Convert to 0-based indexing
        actual_end = end
        
        # Calculate padded region
        padded_start = start - padding
        padded_end = end + padding
        
        params = {
            'original_start': start,
            'actual_start': actual_start,
            'actual_end': actual_end,
            'padded_start': padded_start,
            'padded_end': padded_end,
            'variant_length': len(ref_base)
        }
        
        if mut_class == "DEL":
            params['del_positions'] = list(range(actual_start, end))
        
        return params

    def _open_bam_file(self, bamfile: str):
        """Open BAM or CRAM file."""
        realname = str(Path(bamfile).resolve())
        if realname.endswith("cram"):
            return pysam.AlignmentFile(realname, 'rc')
        else:
            return pysam.AlignmentFile(bamfile, 'rb')

    def _process_insertion(self, samfile, chrom: str, params: dict, min_base_qual: int, min_mapping_qual: int) -> tuple:
        """Process insertion variants."""
        ref_count = alt_count = total_depth = 0
        
        for pileupcolumn in samfile.pileup(chrom, params['padded_start'], params['padded_end'], 
                                        ignore_overlaps=False, min_base_quality=min_base_qual,
                                        min_mapping_quality=min_mapping_qual, truncate=True):
            if pileupcolumn.pos == params['actual_start']:
                seqs = pileupcolumn.get_query_sequences(add_indels=True)
                for base in seqs:
                    if "+" in base:
                        alt_count += 1
                    else:
                        ref_count += 1
                    total_depth += 1
        
        return ref_count, alt_count, total_depth

    def _process_deletion(self, samfile, chrom: str, params: dict, ref_base: str, 
                        min_base_qual: int, min_mapping_qual: int) -> tuple:
        """Process deletion variants."""
        del_support = {pos: 0 for pos in params['del_positions']}
        max_depth = alt_count = 0
        
        for pileupcolumn in samfile.pileup(chrom, params['padded_start'], params['padded_end'],
                                        ignore_overlaps=False, min_base_quality=min_base_qual,
                                        min_mapping_quality=min_mapping_qual, truncate=True):
            if pileupcolumn.pos in params['del_positions']:
                seqs = pileupcolumn.get_query_sequences(add_indels=True)
                if len(seqs) > max_depth:
                    max_depth = len(seqs)
                for base in seqs:
                    if base == "*":
                        del_support[pileupcolumn.pos] += 1
        
        # Get maximum deletion support across all positions
        for position in del_support:
            if del_support[position] > alt_count:
                alt_count = del_support[position]
        
        ref_count = max_depth - alt_count
        total_depth = max_depth
        
        return ref_count, alt_count, total_depth

    def _process_snp_dnp(self, samfile, chrom: str, params: dict, ref_base: str, alt_base: str, 
                        mut_class: str, min_base_qual: int, min_mapping_qual: int) -> tuple:
        """Process SNP and DNP variants."""
        read_bases = {}
        
        for pileupcolumn in samfile.pileup(chrom, params['padded_start'], params['padded_end'],
                                        ignore_overlaps=False, min_base_quality=min_base_qual,
                                        min_mapping_quality=min_mapping_qual, truncate=True):
            if pileupcolumn.pos == params['actual_start']:
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        read_name = pileupread.alignment.query_name
                        sequence = self._extract_variant_sequence(pileupread, mut_class, params['variant_length'])
                        
                        if sequence:  # Only add if we successfully extracted sequence
                            if read_name in read_bases:
                                read_bases[read_name].append(sequence)
                            else:
                                read_bases[read_name] = [sequence]
        
        return self._count_read_support(read_bases, ref_base, alt_base)

    def _extract_variant_sequence(self, pileupread, mut_class: str, variant_length: int) -> str:
        """Extract the appropriate sequence from a read based on variant type."""
        query_pos = pileupread.query_position
        query_sequence = pileupread.alignment.query_sequence
        
        if mut_class == "SNP":
            return query_sequence[query_pos]
        elif mut_class == "DNP":
            # Extract full DNP sequence
            if query_pos + variant_length <= len(query_sequence):
                return query_sequence[query_pos:query_pos + variant_length]
            else:
                return None  # Can't extract full sequence
        else:
            return None

    def _count_read_support(self, read_bases: dict, ref_base: str, alt_base: str) -> tuple:
        """Count read support from collected read bases."""
        ref_count = alt_count = total_depth = 0
        
        for read_name in read_bases:
            bases = read_bases[read_name]
            
            # Handle cases where we have 1 or 2 observations from the same fragment
            if len(bases) == 1:
                total_depth += 1
                if bases[0] == ref_base:
                    ref_count += 1
                elif bases[0] == alt_base:
                    alt_count += 1
            elif bases[0] == bases[1]:
                # Both reads from pair agree
                total_depth += 1
                if bases[0] == ref_base:
                    ref_count += 1
                elif bases[0] == alt_base:
                    alt_count += 1
            else:
                # Reads from pair disagree
                total_depth += 1
                if bases[0] == ref_base or bases[1] == ref_base:
                    ref_count += 1
                elif bases[0] == alt_base and bases[1] == alt_base:
                    alt_count += 1
        
        return ref_count, alt_count, total_depth

    def _format_results(self, ref_count: int, alt_count: int, total_depth: int) -> dict:
        """Format results as strings to prevent pandas float conversion."""
        return {
            'tumor_depth': str(total_depth),
            'tumor_ref_count': str(ref_count),
            'tumor_alt_count': str(alt_count)
        }


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
