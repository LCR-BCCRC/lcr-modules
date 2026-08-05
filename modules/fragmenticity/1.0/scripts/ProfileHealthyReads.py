"""Get length distributions of reads in a bam file.



"""

import pandas as pd
import argparse
import pysam
import random
from collections import Counter

def argument_parser():
    parser = argparse.ArgumentParser(description="Get length distributions of reads in a bam file.")
    parser.add_argument("--bam", type=str, help="BAM file to analyze")
    parser.add_argument("-l", "--min_length", type=int, default=1, help="Minimum fragment length to consider")
    parser.add_argument("-L", "--max_length", type=int, default=900, help="Maximum fragment length to consider")
    parser.add_argument("-c", "--read_count", type=int, default=100000, help="Number of reads to sample")
    parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("-o", "--output", type=str, default="length_distribution.tsv", help="Output TSV file")
    return parser


class ProfileHealthyReads:
    def __init__(self, bam_file, output_file, read_count, seed, min_length=1, max_length=900):
        # user inputs
        self.bam_file = bam_file
        self.output_file = output_file
        self.min_length = min_length
        self.max_length = max_length
        self.read_count = read_count
        self.seed = seed
        # outputs

    def sample_fragment_lengths(self) -> list:
        """
        Two-pass approach for random sampling without replacement
        """
        random.seed(self.seed)
        
        # First pass: collect all valid fragment lengths
        valid_lengths = []
        total_processed = 0
        processed_UMIs = set() # for tracking read mates

        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            for read in bam.fetch():
                total_processed += 1
                # Skip poorl mapped reads
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < 20:
                    continue

                # UMI-based duplicate avoidance
                try:
                    um = read.get_tag("MI")
                    if um in processed_UMIs:
                        continue
                    processed_UMIs.add(um)
                except KeyError:
                    continue

                # Filter only by read length
                template_length = abs(read.template_length)
                if self.min_length <= template_length <= self.max_length:
                    valid_lengths.append(template_length)

        print(f"Found {len(valid_lengths)} valid reads out of {total_processed} total")
        
        # If we have fewer valid reads than requested, return all
        if len(valid_lengths) <= self.read_count:
            print(f"Using all {len(valid_lengths)} valid reads (less than requested {self.read_count})")
            return valid_lengths
        
        # Randomly sample from valid lengths
        sampled_lengths = random.sample(valid_lengths, self.read_count)
        
        print(f"Randomly sampled {len(sampled_lengths)} reads")

        length_counts = Counter(sampled_lengths)
        print(f"Most common lengths: {length_counts.most_common(20)}")
        return sampled_lengths

    def summarize_read_length_distribution(self):
        """
        Summarize the read length distribution by sampling fragment lengths
        """

        sampled_lengths = self.sample_fragment_lengths()
        
        # Create a DataFrame for the sampled lengths
        df = pd.DataFrame(sampled_lengths, columns=["length"])
        
        # Count occurrences of each length
        length_counts = df["length"].value_counts().sort_index()
        
        # Create a DataFrame for the counts
        summary_df = pd.DataFrame(length_counts).reset_index()
        summary_df.columns = ["length", "count"]
        
        return summary_df

def main():
    args = argument_parser().parse_args()

    readcounter = ProfileHealthyReads(
        bam_file=args.bam,
        output_file=args.output,
        min_length=args.min_length,
        max_length=args.max_length,
        read_count=args.read_count,
        seed=args.seed
    )
    summary_df = readcounter.summarize_read_length_distribution()
    summary_df.to_csv(args.output, sep="\t", index=False)
    print("Done!" + '\033[92m' + " ✔" + '\033[0m')

if __name__ == "__main__":
    main()