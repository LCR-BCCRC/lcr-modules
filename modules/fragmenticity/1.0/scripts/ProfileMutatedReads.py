"""Find mutated and wildtype reads and collect metadata on them.

optionally exclude chip genes

Mutation info is pulled from provided maf file.

Data collected:
5' and 3' end motiff count (4 bps)
tlength
mutated or wildtype
position of mutation
mutation p. notation

Only a distribution of read lengths is returned.

This script is designed to process a single sample at a time. For
parallel processing see accompanying smk. Or write one yourself.
"""
import pandas as pd
import pysam as ps
import glob
import argparse
import tqdm

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', type=str,help='BAM file to process')
    parser.add_argument('--maf', help='MAF file to process')
    parser.add_argument('--output', help='Output file name')
    parser.add_argument("--exclude_chip_genes", action='store_true',help="Exclude chip genes from the analysis")
    parser.add_argument("--save_wildtype", action='store_true',
                        help="Save wildtype read length distribution as well")
    return parser.parse_args()

class ProfileReads(object):
    """"""

    def __init__(self, maf_path, bam_path, exclude_chip_genes=False):
        super(ProfileReads, self).__init__()
        self.mutations = pd.read_csv(maf_path, sep='\t')
        self.exclude_chip_genes = exclude_chip_genes
        self.chip_genes = ['TET2', 'DNMT3A', 'ASXL1', 'SF3B1','IDH2'] # keep TP53 cause it commonly mutated

        self.bam = bam_path # only open when you need it
        self.read_metadata = self.profile_reads()

    def profile_reads(self):
        """Loop through mutations and find reads that contain them.

        Returns a dataframe with all reads that contain mutations.
        """

        # if no mutations return empty df
        if self.mutations.empty:
            return pd.DataFrame()

        # if exclude chip genes is set, filter out chip genes
        if self.exclude_chip_genes:
            print("Excluding chip genes from analysis...")
            self.mutations = self.mutations[~self.mutations['Hugo_Symbol'].isin(self.chip_genes)].copy()
            if self.mutations.empty:
                print("No mutations left after excluding chip genes. Exiting.")
                return pd.DataFrame()

        # loop through mutation in maf file df
        all_mutation_reads = []
        for mutation in self.mutations.itertuples():
            # skip indels
            if len(mutation.Reference_Allele) > 1 or len(mutation.Tumor_Seq_Allele2) > 1 or mutation.Reference_Allele == "-" or mutation.Tumor_Seq_Allele2 == "-":
                all_mutation_reads.append(pd.DataFrame()) # if indel just return empty df
                continue
            print(f"Processing mutation: {mutation.Hugo_Symbol} {mutation.Start_Position}:{mutation.End_Position} {mutation.Reference_Allele}->{mutation.Tumor_Seq_Allele2}")
            # find reads that contain mutation
            reads = self.find_mutated_read(mutation)
            all_mutation_reads.append(reads)

        return pd.concat(all_mutation_reads)
        
    def find_mutated_read(self, mutation):
        """Search all reads in bam file for mutation and collect metadata.

        This function needs to be chunked up into smaller functions.
        """

        bamfile = ps.AlignmentFile(self.bam, 'rb')
        # get read info
        reads = bamfile.fetch(mutation.Chromosome, mutation.Start_Position-1, mutation.End_Position)

        all_reads = []
        # add status bar 
        pbar = tqdm.tqdm(reads)
        pbar.set_description("Processing reads")
        for read in reads:
            # add mutation info to read dict, one read per dict
            read_info = {}
            read_info["chromosome"] = mutation.Chromosome
            read_info["Hugo_Symbol"] = mutation.Hugo_Symbol
            read_info["mut_start"] = mutation.Start_Position
            read_info["mut_end"] = mutation.End_Position
            read_info["ref"] = mutation.Reference_Allele
            read_info["alt"] = mutation.Tumor_Seq_Allele2
            # check if read is mutated or wildtype or other
            # PYSAM RETURNS 0 BASED POSITIONS, wherease MAF is 1 based
            seqdf = pd.DataFrame(read.get_aligned_pairs(), columns=["read_pos (0 based)", "ref_pos (0 based)"])
            # drop rows with NaN
            seqdf = seqdf.dropna()
            # make sure each col type is int
            seqdf = seqdf.astype(int)
            #print(seqdf.loc[(seqdf["ref_pos (0 based)"] > mutation.Start_Position-3) & (seqdf["ref_pos (0 based)"] < mutation.End_Position+3), :])
            # find mutation position
            mut_pos_df = seqdf[seqdf["ref_pos (0 based)"] == mutation.Start_Position-1].reset_index(drop=True).copy()
            # if df is empty skip read
            if mut_pos_df.empty:
                continue
            #print(mut_pos_df)
            mut_read_pos = int(mut_pos_df.iloc[0]["read_pos (0 based)"])
            # get read allele aka query
            # if the read is shorter than the mutation position skip it, the position was probably soft clipped
            if len(read.query_alignment_sequence) < mut_read_pos+1:
                continue
            query_allele = read.query_alignment_sequence[mut_read_pos]
            # print(f"Mutation position: {mut_pos}")
            # print(read.query_alignment_sequence)
            read_info["observed"] = query_allele

            if read_info["observed"] == mutation.Reference_Allele:
                mut = "ref"
            elif read_info["observed"] == mutation.Tumor_Seq_Allele2:
                mut = "alt"
            else:
                mut = "other"

            read_info["mutated"] = mut
            # get read length
            read_info["tlength" ]= abs(read.template_length)
            read_info["rlength"] = len(read.query_sequence)
            # read_info["qlen"] = read.infer_read_length()
            # get 5' end motiff
            read_info["start_motif"] = read.query_sequence[:4]
            # get 3' end motiff
            read_info["end_motif"] = read.query_sequence[-4:]
            all_reads.append(read_info)
            pbar.update(1)

        bamfile.close()
        pbar.close()

        return pd.DataFrame(all_reads)

    def summarize_read_length_distribution(self):
        """
        Summarize the read length distribution by mutation status.
        Returns separate distributions for mutated and wildtype reads.
        """
        if self.read_metadata.empty:
            return pd.DataFrame(), pd.DataFrame()

        # Separate mutated (alt) and wildtype (ref) reads
        mutated_reads = self.read_metadata[self.read_metadata['mutated'] == 'alt']
        wildtype_reads = self.read_metadata[self.read_metadata['mutated'] == 'ref']
        
        # Create length distributions for each
        if not mutated_reads.empty:
            mutated_distribution = mutated_reads.groupby("tlength").size().reset_index(name='count')
            mutated_distribution.rename(columns={'tlength': 'length'}, inplace=True)
            mutated_distribution = mutated_distribution.sort_values(by="length").reset_index(drop=True)
        else:
            mutated_distribution = pd.DataFrame(columns=['length', 'count'])
        
        if not wildtype_reads.empty:
            wildtype_distribution = wildtype_reads.groupby("tlength").size().reset_index(name='count')
            wildtype_distribution.rename(columns={'tlength': 'length'}, inplace=True)
            wildtype_distribution = wildtype_distribution.sort_values(by="length").reset_index(drop=True)
        else:
            wildtype_distribution = pd.DataFrame(columns=['length', 'count'])
        
        return mutated_distribution, wildtype_distribution

    def summarize_read_length_distribution_by_gene(self):
        """
        Summarize read length distribution with gene information.
        Returns a single DataFrame with length, count, and gene columns.
        """
        if self.read_metadata.empty:
            return pd.DataFrame(columns=['length', 'count', 'gene'])

        # Get only mutated reads and group by gene and length
        mutated_reads = self.read_metadata[self.read_metadata['mutated'] == 'alt']
        
        if mutated_reads.empty:
            return pd.DataFrame(columns=['length', 'count', 'gene'])
        
        # Group by gene and length, count occurrences
        gene_length_distribution = mutated_reads.groupby(['Hugo_Symbol', 'tlength']).size().reset_index(name='count')
        gene_length_distribution.rename(columns={'tlength': 'length', 'Hugo_Symbol': 'gene'}, inplace=True)
        gene_length_distribution = gene_length_distribution.sort_values(['gene', 'length']).reset_index(drop=True)
        
        # calc relative frequency
        gene_length_distribution['relative_frequency'] = gene_length_distribution['count'] / gene_length_distribution.groupby('gene')['count'].transform('sum')

        return gene_length_distribution

def main():
    args = get_args()
    print('\033[95m' + "Running ProfileMutatedReads.py...")
    print("Processing BAM file: " + args.bam)
    print("Processing MAF file: " + args.maf + '\033[0m')
    
    read_metadata = ProfileReads(args.maf, args.bam, args.exclude_chip_genes)

    # Generate gene-specific length distributions
    gene_distribution = read_metadata.summarize_read_length_distribution_by_gene()
    
    # Save the gene-specific distribution
    gene_output = args.output.replace(".tsv", "_by_gene.tsv")
    gene_distribution.to_csv(gene_output, sep='\t', index=False)
    print(f"Wrote gene-specific distribution: {gene_output}")
    

    # Generate separate length distributions for mutated and wildtype reads
    mutated_distribution, wildtype_distribution = read_metadata.summarize_read_length_distribution()
    
    # Save mutated distribution to main output
    print(f'\033[95mWriting mutated read length distribution to: {args.output}\033[0m')
    mutated_distribution.to_csv(args.output, sep='\t', index=False)
    
    # Optionally save wildtype distribution
    if args.save_wildtype:
        wildtype_output = args.output.replace(".tsv", "_wildtype_read_length_distribution.tsv")
        print(f"Saving wildtype read length distribution to: {wildtype_output}")
        wildtype_distribution.to_csv(wildtype_output, sep='\t', index=False)
    
    print("Done!" + '\033[92m' + " ✔" + '\033[0m')

if __name__ == '__main__':
    main()