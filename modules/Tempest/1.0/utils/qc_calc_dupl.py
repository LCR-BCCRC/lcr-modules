import argparse
import pysam
import fastq_utils as fu
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--collapsed_bam',required=True,type=str,help='')
    parser.add_argument('--all_reads_bam',required=True,type=str,help='')
    parser.add_argument('--samplename',required=True,type=str,help='')
    parser.add_argument('--output',required=True,type=str,help='')

    return parser.parse_args()

def qc_calc_dupl(collapsed_bam: str, all_reads_bam: str, samplename: str , output_txt: str) -> None:
    # This definitely isn't the most efficient way to calculate duplicate rate, but it works
    # We are going to generate a Picard-style output file
    # First, parse the final (collapsed) BAM, and calculate the total number of reads
    col_bam = pysam.AlignmentFile(collapsed_bam)
    consensus_reads = col_bam.count(until_eof=True, read_callback=lambda x: not x.is_duplicate and x.is_mapped and x.is_paired and not x.is_supplementary and not x.is_secondary)
    col_read_pairs = consensus_reads / 2
    col_unpaired_reads = col_bam.count(until_eof=True, read_callback=lambda x: not x.is_paired and not x.is_supplementary and not x.is_secondary)

    # Now, parse the original (non-consensus) BAM, and calculate the total number
    # of read pairs, unmapped reads, and unpaired reads
    orig_bam = pysam.AlignmentFile(all_reads_bam, require_index=False)
    orig_reads = orig_bam.count(until_eof=True, read_callback=lambda x: not x.is_duplicate and x.is_mapped and x.is_paired and not x.is_supplementary and not x.is_secondary)
    orig_pairs = orig_reads / 2
    unmapped_reads = orig_bam.count(until_eof=True, read_callback=lambda x: x.is_unmapped and not x.is_supplementary and not x.is_secondary)
    unpaired_reads = orig_bam.count(until_eof=True, read_callback=lambda x: not x.is_paired and not x.is_supplementary and not x.is_secondary)
    secondary_reads = orig_bam.count(until_eof=True, read_callback=lambda x: x.is_supplementary or x.is_secondary)

    # Now, calculate the output stats
    library = samplename
    read_pairs_examined = int(orig_pairs)
    secondary_or_supplemental = secondary_reads
    unpaired_dups = int(unpaired_reads - col_unpaired_reads)
    read_pair_dups = int(orig_pairs - col_read_pairs)
    optical_dup = 0  # In this consensus approach (via UMIs) we lose info on which are optical duplicates
    per_dupl = read_pair_dups / orig_pairs
    # Custom added by Chris. Calculate the total size of this library, and how much we have sequenced
    estimated_library_size = int(fu.estimateLibrarySize(orig_pairs, col_read_pairs))  # Use MultiQC's function to calculate the library size
    prop_library_seq = col_read_pairs / estimated_library_size

    # Close input
    col_bam.close()
    orig_bam.close()

    # Write output
    with open(f"{output_txt}", "w") as o:
        # To trick MultiQC into thinking this is a Picard output
        o.write("##picard.sam.markduplicates.MarkDuplicates BUT NOT ACTUALLY PICARD")
        o.write(os.linesep)
        header = ["LIBRARY", "UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "SECONDARY_OR_SUPPLEMENTARY_RDS", "UNMAPPED_READS", "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES",
                    "READ_PAIR_OPTICAL_DUPLICATES", "PERCENT_DUPLICATION", "ESTIMATED_LIBRARY_SIZE", "PERCENT_LIBRARY_SEQUENCED"]
        o.write("\t".join(header))
        o.write(os.linesep)
        # Write out stats for this sample
        out_values = [library, str(unpaired_reads), str(read_pairs_examined), str(secondary_or_supplemental), str(unmapped_reads),
                        str(unpaired_dups), str(read_pair_dups), str(optical_dup), str(per_dupl), str(estimated_library_size), str(prop_library_seq)]
        o.write("\t".join(out_values))
        o.write(os.linesep)


def main():
    args = get_args()
    qc_calc_dupl(args.collapsed_bam, args.all_reads_bam, args.samplename, args.output)

if __name__ == '__main__':
    main()