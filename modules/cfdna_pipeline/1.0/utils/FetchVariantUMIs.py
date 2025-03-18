"""Script for fetching the UMI tags of reads that were used to call variants.
(FYI that is more difficult than it sounds)

The UMI family sizes for each read will be added to the maf file as two new columns.

UMI_mean: the mean UMI family size of the variant reads
UMI_max: the maximum UMI family size of the variant reads
STR: the proportion of reads that are have a short tandem repeat region near the variant site
    in the reference genome.


The is script first:
1. Recalculates the AF column based on the new set of variants.
2. Loops through each variant in MAF to determine the UMI family size of the variant reads.
3. Adds the UMI family sizes mean and max to the MAF file.
4. Adds a STR column to the MAF file, which is the proportion of reads that are in a short tandem repeat region.

To determine the UMI family size of the variant reads, the script:
1. Loops through each read that overlaps with the variant region.
2. Checks if the variant is in a short tandem repeat region.
3. Checks if the read has the variant allele positions (start and end)
4. Determines the type of variant (SNV, INS, DEL, DNV)
5. Expands the read sequence by the cigar string, adding gaps as "-", and includes softclipped bases.
6. Determines if the mutant allele is in the read sequence
7. If the allele is in the read, the UMI family size is added to the list of UMI family sizes.


Exclaimers:
Variant callers make a lot of decisions under the hood, so it can be hard to determine what reads were used to make the call.
For example, some insertions use reads that have additional bases in the insert, but the call will only involve a subset of those bases.

When pysam/samtools opens a bam file it converts it to a sam file, which is 1-based. So need to adjust for 0-based variant positions made by caller.

Some callers, like SAGE, sometimes use softclipped bases to call insertions. This script will include those bases in the read sequence.
"""
import pysam
import pandas as pd
import argparse
import re
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_maf',required=True,type=str,help='')
    parser.add_argument('--input_bam',required=True,type=str,help='')
    parser.add_argument('--output_maf',required=True,type=str,help='')
    return parser.parse_args()

def recalculate_af(df:pd.DataFrame) -> pd.DataFrame:
    """Recalculate the AF column based on the new set of variants.

    Always good to be sure! Especially if this is run downstream
    of a script like augment_ssm!

    Args:
        df (pd.DataFrame): DataFrame of variants.
    Returns:
        pd.DataFrame: DataFrame of variants with recalculated AF column.
    """
    # calculate the AF
    df["AF"] = df["t_alt_count"] / df["t_depth"]
    return df.copy()

def UMI_Fam_size(bamfile:str, chrom:str, start:int, end:int, allele: str, ref_allele: str) -> tuple:
    """Get the UMI tags of reads that overlap with the variant region.

    Use get_aligned_pairs to get the position of the variant in the read.
    Also tracks if the variant is in a short tandem repeat region.

    Pysam shifts bam files into sam files, which is 1-based. So need to adjust
    for 0-based variant positions.

    Args:
        bamfile (str): Path to the bam file.
        chrom (str): Chromosome of the variant.
        start (int): Start position of the variant.
        end (int): End position of the variant.
        allele (str): Allele of the variant.
        ref_allele (str): Reference allele of the variant.
    Returns:
        list: UMI family sizes of the variant reads.
        float: Proportion of reads that are in a short tandem repeat region.
        
    """
    umi_depths = []
    processed_UMIs = set()
    STR = []
    var_type = __determine_var_type(allele, ref_allele)

    if var_type == "OTHER":
        print("Variant type not recognized, skipping")
        return [],0

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam.fetch(chrom, start-1, end):

            try:
                # check UMI tag to see if mate has been processed
                um = read.get_tag("MI")
                if um in processed_UMIs:
                    continue
                else:
                    processed_UMIs.add(um)

            except KeyError: # skip read if no UMI tag
                continue

            cig = read.cigarstring
            if cig is None:
                continue

            # get expanded seq with indels and softclips
            full_seq = expand_seq_wcigar(read.query_sequence, cig)

            apairsdf = pd.DataFrame(read.get_aligned_pairs(matches_only=False ,with_seq=True), columns=["readpos", "refpos", "refbase"])
            # add full seq to 
            apairsdf["readbase"] = list(full_seq)
            # make refbase "-" if None
            apairsdf.fillna({"refbase":"-"}, inplace=True)
            # change base colums to uppercase
            apairsdf["refbase"] = apairsdf["refbase"].str.upper()
            apairsdf["readbase"] = apairsdf["readbase"].str.upper()

            # Ignore reads pulled that done have proper alignment for var type, varies for vartype
            if var_type == "INS" and start+len(allele) not in apairsdf["refpos"].values:
                continue
            if (var_type == "DEL" and (end-1 not in apairsdf["refpos"].values or
                start-1 not in apairsdf["refpos"].values)):
                continue
            # need previous position to find alleles for deletions and snvs, dnvs
            if var_type in ["SNV", "DNV", "INS", "TNV"] and start-1 not in apairsdf["refpos"].values:
                continue
            if var_type in ["DNV", "TNV"] and start not in apairsdf["refpos"].values:
                continue

            # get position of variant in reference
            ref_pos =  int(start - apairsdf["refpos"].min())

            # check if short tandem repeat region
            STR_STATUS = determine_STR_status(ref_pos, allele, read.get_reference_sequence())
            if STR_STATUS:
                STR.append(1)
            else:
                STR.append(0)
            # insertions
            if var_type == "INS":
                # get mutated read UMI family size
                row_index_start = apairsdf[apairsdf["refpos"] == start-1].index
                row_index_end = apairsdf[apairsdf["refpos"] == start].index
                # get the dataframe of the rows between the start and end indexes
                allele_df = apairsdf.loc[row_index_start[0]+1:row_index_end[0]-1].copy()
            # snvs
            elif var_type == "SNV":
                # get mutated read UMI family size
                row_index_start = apairsdf[apairsdf["refpos"] == start-1].index
                # get the dataframe of the rows between the start and end indexes
                allele_df = apairsdf.loc[row_index_start[0]:row_index_start[0]].copy()
            # dnvs and tnvs
            elif var_type == "DNV" or var_type == "TNV":
                # get mutated read UMI family size
                row_index_start = apairsdf[apairsdf["refpos"] == start-1].index
                row_index_end = apairsdf[apairsdf["refpos"] == start].index
                # get the dataframe of the rows between the start and end indexes
                allele_df = apairsdf.loc[row_index_start[0]:row_index_end[0]].copy()
            # deletions
            elif var_type == "DEL":
                # get mutated read UMI family size
                row_index_start = apairsdf[apairsdf["refpos"] == start-1].index
                row_index_end = apairsdf[apairsdf["refpos"] == end-1].index
                # get the dataframe of the rows between the start and end indexes
                allele_df = apairsdf.loc[row_index_start[0]:row_index_end[0]].copy()

            # join the readbase column to get the full allele
            read_allele = "".join(allele_df["readbase"].tolist())
            # sometimes INS called are only part of the insert in a individual read
            if read_allele == allele or (var_type == "INS" and allele in read_allele):
                if read.has_tag("cD"):
                    umi_depth = read.get_tag("cD")
                    umi_depths.append(umi_depth)

    if len(STR) == 0:
        str_score = 0
    else:
        str_score = sum(STR)/len(STR)
    return umi_depths, str_score

def __determine_var_type(alt_allele:str, ref_alele: str)-> str:
    """Determine the type of variant.

    Based on alt and ref alleles. This support function is to
    standardize annotations used to support the UMI family size
    determination.
    """
    # SNV
    if len(alt_allele) == 1 and len(ref_alele) == 1 and "-" not in alt_allele and ref_alele != "-":
        return "SNV"
    # insertion
    if "-" in ref_alele:
        return "INS"
    # deletion
    if "-" in alt_allele:
        return "DEL"
    # dinucleotide variant
    if len(alt_allele) == 2 and len(ref_alele) == 2:
        return "DNV"
    # trinucleotide variant
    if len(alt_allele) == 3 and len(ref_alele) == 3:
        return "TNV"
    else:
        print(f"Variant type not recognized: {ref_alele} -> {alt_allele}, sorry about it")
        return "OTHER"

def determine_STR_status(start: int, allele: str, seq: str) -> bool:
    """Determine if a variant occurs in a short tandem repeat.

    If the variant is in a short tandem repeat, the variant is likely to be a false positive.
    Only considers consecutive repeating patterns as STRs.

    Args:
        start (int): start position of the variant.
        allele (str): variant sequence
        seq (str): full reference sequence where the read aligns to.

    Returns:
        bool: is the variant in a short tandem repeat?
    """
    # Handle empty allele case
    if not allele:
        return False
    
    # Check for special regex characters and skip with warning
    symbols = ['+', '*', '?', '.']
    if any(sym in seq for sym in symbols):
        print(f"Warning: Skipping STR check for seq with special characters: '{seq}'")
        return np.nan
    
    # Use a larger window for better detection
    window_size = max(10, len(allele) * 4)  # Larger window for better context
    start_buffed = max(0, start - window_size)
    end_buffed = min(len(seq), start + window_size)  # Prevent index errors
    
    ref_allele = seq[start:start+len(allele)]
    seq_by_var = seq[start_buffed:end_buffed].upper()
    
    # If single base, check for runs of 3+ of the same base
    if len(ref_allele) == 1:
        pattern = re.compile(f"{ref_allele}{{{3},}}")
        if pattern.search(seq_by_var):
            return True
        return False

    # Special homopolymer check - detect cases like "AAAAAA" in "GCGCAAAAAAAAGCGC"
    # if allele is homopolymer longer than 3 bases, return True
    # that means it must appear in the ref sequence, so need to check it
    if len(set(allele)) == 1 and len(allele) >= 3:
        return True

    # Check if the allele appears consecutively in the window
    if len(ref_allele) >= 2:
        pattern = re.compile(f"({re.escape(ref_allele)}){{{2},}}")
        if pattern.search(seq_by_var):
            return True
    
    # Check for consecutive repeats of half patterns
    if len(ref_allele) > 2:
        half = len(ref_allele) // 2
        first_half = ref_allele[:half]
        second_half = ref_allele[half:]
        
        # Must have 3+ consecutive repeats of the half
        if len(first_half) > 0:
            pattern = re.compile(f"({re.escape(first_half)}){{{3},}}")
            if pattern.search(seq_by_var):
                return True
        
        if len(second_half) > 0:
            pattern = re.compile(f"({re.escape(second_half)}){{{3},}}")
            if pattern.search(seq_by_var):
                return True
    
    return False

def add_umi_support(inmaf:pd.DataFrame, bamfile:str) -> pd.DataFrame:
    """Get the UMI family sizes for each variant in the maf file.

    Args:
        inmaf (pd.DataFrame): DataFrame of variants.
        bamfile (str): Path to the bam file.
    Returns:
        pd.DataFrame: DataFrame of variants with UMI family sizes.
    """
    var_count = 1
    var_total = inmaf.shape[0]

    for index, row in inmaf.iterrows():
        chrom = str(row["Chromosome"])
        start = row["Start_Position"]
        end = row["End_Position"]
        ref_allele = row["Reference_Allele"]
        allele = row["Tumor_Seq_Allele2"]
        t_depth = row["t_depth"]
        # allele needs to be the exact number of positions
        if allele == "-":
            allele = "-" * len(ref_allele)

        print(f"Processing variant {var_count} of {var_total}")
        print(f"Chrom: {chrom}, Start: {start}, End: {end}, Ref: {ref_allele}, Alt: {allele}")
        var_count += 1

        umi_depths, str_score = UMI_Fam_size(bamfile, chrom, start, end, allele, ref_allele)

        # no UMI support found
        # or variant was in repeat region, so coudlnt determine what reads had the variant
        if len(umi_depths) == 0 or (len(umi_depths) >= (0.6 * t_depth) and str_score > 0.6):
            inmaf.at[index, "UMI_mean"] = 0
            inmaf.at[index, "UMI_max"] = 0
            inmaf.at[index, "STR"] = str_score
            continue
        inmaf.at[index, "UMI_mean"] = sum(umi_depths) / len(umi_depths)
        inmaf.at[index, "UMI_max"] = max(umi_depths)
        inmaf.at[index, "STR"] = str_score
        

    return inmaf.copy()

def expand_seq_wcigar(seq: str, cigar: str) -> str:
    """Expand the sequence by the cigar string.

    Args:
        seq (str): Sequence of the read.
        cigar (str): Cigar string of the read.

    Returns:
        str: Expanded sequence with gaps as -.
    """
    expanded_seq = ""
    bases_used = 0
    # split cigar by letters

    cigar_list = re.findall(r'\d+[A-Za-z^]*', cigar)

    for c in cigar_list:

        if c.endswith("M") or c.endswith("I") or c.endswith("S"): # uses query bases
            bases = int(re.sub(r'[A-Za-z]', '', c))
            expanded_seq += seq[bases_used: bases+bases_used]
            bases_used += bases
        elif c.endswith("D"):
            bases = int(c.replace("D", ""))
            expanded_seq += "-"*bases
        elif c.endswith("N"):
            bases = int(c.replace("N", ""))
            expanded_seq += "-"*bases

    return expanded_seq.upper()

def main():
    args = get_args()
    inmaf = pd.read_csv(args.input_maf, sep='\t')
    # recalc AF
    inmaf = recalculate_af(inmaf)
    # add UMI support
    inmaf = add_umi_support(inmaf, args.input_bam)
    inmaf.to_csv(args.output_maf, sep='\t', index=False)

if __name__ == '__main__':
    main()
