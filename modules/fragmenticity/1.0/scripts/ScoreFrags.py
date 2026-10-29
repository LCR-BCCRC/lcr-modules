"""Script for sampling reads from bams and calculating their fragmentation score (FS)

The FS score is calculate by sampling reads from the bam file, then for each read take the length and look
up the score in a reference file, then the mean of all scampled reads for a given sample is reported as the FS.

The reference file should contain one column with no headers, the column is the fragment score (1-900)

The output is a TSV file with the sample name and the fragmentation score.

The script uses pysam to read the BAM file and pandas to handle the data. And for running multiple samples use the snakemake worfklow.

Minimal example:

python /home/kyakimovich/repos/kurt/modules/LymphoFrag/scripts/ScoreFrags.py \
    --bam_file /projects/rmorin_scratch/prospective_rrDLBCL_trial_data/pipeline_outputs/bam_pipeline/99-final/JGH_001_DX.consensus.mapped.annot.bam\
    --fragment_score_reference /home/kyakimovich/repos/kurt/modules/LymphoFrag/FragmentLengthReference_Vessies.tsv \
    --sample_name JGH_001 \
    --output /home/kyakimovich/repos/

Written by: Kurt Yakimovich
"""

import pandas as pd
import glob
import os
import pysam
import argparse
import random
import plotnine as pn
from collections import Counter
import warnings

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate fragmentation scores from BAM files.")
    parser.add_argument("--bam_file", required=True, help="Path to the input BAM file.")
    parser.add_argument("--fragment_score_reference", type=str, required=True, help="Path to the reference file with scores for fragment lengths 1-900, one column no headers")
    parser.add_argument("--sample_name", type=str, default=None, help="Sample name to use in the output file, default is the BAM file name without extension.")
    parser.add_argument("--output", default="./", help="Output directory")
    parser.add_argument("--read_count", type=int, default=1000000, help="Number of reads to sample from the BAM file, default 1 million reads.")
    parser.add_argument("--min_length", type=int, default=50, help="Minimum length of reads to be considered, default 50.")
    parser.add_argument("--max_length", type=int, default=900, help="Maximum length of reads to be considered, default 900.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility, default 42.")

    return parser.parse_args()

class FragScorer(object):
    def __init__(self, bam_file, fragment_score_reference=None, read_count=1000000, min_length=50, max_length=900, seed=42):
        self.bam_file = bam_file
        self.fragment_score_reference = fragment_score_reference
        self.read_count = read_count
        self.min_length = min_length
        self.max_length = max_length
        self.seed = seed

        self.sampled_lengths = []
        self.all_valid_lengths = []  # Store all valid lengths for full histogram
        self.actual_sampled_count = 0  # Track how many reads we actually sampled
        self.FS = self.calculate_fragmentation_score()

    def read_score_reference(self):
        """Reads the fragment score reference file and creates a dictionary."""
        scores_dict = {}
        # read in score reference file as .tsv
        raw_ref  = pd.read_csv(self.fragment_score_reference, sep="\t")
        # check if the file has a single column with no headers
        if raw_ref.shape[1] == 1 and "0" in raw_ref.columns:
            # turn that single column into a dictionary, with index as key and value as the score
            scores_dict = {i + 1: score for i, score in enumerate(raw_ref.iloc[:, 0])}
        # in format from LymphoFrag code
        elif "fragmentation_score" in raw_ref.columns and "length" in raw_ref.columns:
            # if the file has two columns, one for length and one for score, length as key and score as value
            scores_dict = raw_ref.set_index("length")["fragmentation_score"].to_dict()
        else:
            raise ValueError("""Fragment score reference file must have one column with no 
                                headers or two columns with 'length' and 'fragmentation_score' headers.""")
        
        print(f"Loaded scores for fragment lengths 1-{len(scores_dict)} from reference file")
        return scores_dict

    def sample_fragment_lengths(self) -> list:
        """
        Two-pass approach for random sampling without replacement.
        Also collects ALL valid lengths (50-900bp) for histogram display.
        """
        random.seed(self.seed)
        
        # First pass: collect all valid fragment lengths for the specified range AND for full histogram
        valid_lengths_for_sampling = []  # Only lengths in min_length to max_length range
        all_valid_lengths = []  # All lengths 50-900bp for histogram
        total_processed = 0
        processed_UMIs = set() # for tracking read mates

        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            for read in bam.fetch():
                total_processed += 1
                # Skip poorly mapped reads
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

                template_length = abs(read.template_length) # make sure its positive, regardless of orientation
                
                # Collect for full histogram (50-900bp)
                if 50 <= template_length <= 900:
                    all_valid_lengths.append(template_length)
                
                # Collect for sampling (specified range)
                if self.min_length <= template_length <= self.max_length:
                    valid_lengths_for_sampling.append(template_length)
        
        # Store all valid lengths for histogram
        self.all_valid_lengths = all_valid_lengths
        
        print(f"Found {len(valid_lengths_for_sampling)} valid reads in range {self.min_length}-{self.max_length}bp out of {total_processed} total")
        print(f"Found {len(all_valid_lengths)} total valid reads in range 50-900bp for histogram")
        
        # **WARNING CHECK**: If we have fewer reads than requested (any shortfall triggers warning)
        if len(valid_lengths_for_sampling) < self.read_count:
            warning_msg = (f"⚠️  WARNING: Only found {len(valid_lengths_for_sampling)} valid reads "
                          f"in range {self.min_length}-{self.max_length}bp, but {self.read_count} were requested. "
                          f"Continuing with available reads.")
            print(warning_msg)
            warnings.warn(warning_msg, UserWarning)
            
            self.actual_sampled_count = len(valid_lengths_for_sampling)
            self.sampled_lengths = valid_lengths_for_sampling
            return valid_lengths_for_sampling
        
        # Randomly sample from valid lengths
        sampled_lengths = random.sample(valid_lengths_for_sampling, self.read_count)
        self.actual_sampled_count = len(sampled_lengths)
        
        print(f"Successfully sampled {len(sampled_lengths)} reads from range {self.min_length}-{self.max_length}bp")

        self.sampled_lengths = sampled_lengths
        length_counts = Counter(self.sampled_lengths)
        print(f"Most common lengths in sampled reads: {length_counts.most_common(10)}")
        return sampled_lengths

    def calculate_fragmentation_score(self):
        """Calculate the fragmentation score based on sampled reads."""
        read_lengths = self.sample_fragment_lengths()
        
        # Get scores as a dictionary
        scores_dict = self.read_score_reference()
        
        # Calculate scores for each fragment length
        fs_scores = [scores_dict[length] for length in read_lengths]

        # Calculate the mean fragmentation score
        fragmentation_score = sum(fs_scores) / len(fs_scores)
        median_score = sorted(fs_scores)[len(fs_scores)//2]  # Simple median

        print(f"Calculated fragmentation score from {len(fs_scores)} reads: {fragmentation_score}")
        print(f"Median fragmentation score: {median_score}")

        return fragmentation_score

    def make_histogram(self, sample_id, outpath):
        """Create a smooth histogram using plotnine with proper grouping."""
        
        # Prepare data for all valid reads
        all_counts = pd.Series(self.all_valid_lengths).value_counts().reset_index()
        all_counts.columns = ['length', 'count']
        all_counts['relative_frequency'] = all_counts['count'] / all_counts['count'].sum()
        all_counts['read_type'] = 'All Valid Reads'
        
        # Prepare data for sampled reads  
        if len(self.sampled_lengths) > 0:
            sampled_counts = pd.Series(self.sampled_lengths).value_counts().reset_index()
            sampled_counts.columns = ['length', 'count']
            sampled_counts['relative_frequency'] = sampled_counts['count'] / sampled_counts['count'].sum()
            sampled_counts['read_type'] = 'Sampled Reads'
            
            # Combine datasets
            plot_data = pd.concat([all_counts, sampled_counts], ignore_index=True)
        else:
            plot_data = all_counts
        
        # Create the plot
        plot = (pn.ggplot(plot_data, pn.aes(x='length', y='relative_frequency', 
                                        color='read_type', group='read_type')) +
                pn.geom_line(size=1, alpha=0.8) +
                
                # Add vertical lines for sampling range
                pn.geom_vline(xintercept=self.min_length, color='red', linetype='dashed', alpha=0.8) +
                pn.geom_vline(xintercept=self.max_length, color='red', linetype='dashed', alpha=0.8) +
                
                # Add reference lines
                pn.geom_vline(xintercept=167, color='orange', linetype='dotted', alpha=0.6) +
                pn.geom_vline(xintercept=333, color='purple', linetype='dotted', alpha=0.6) +
                
                pn.theme_bw() +
                pn.theme(figure_size=(12, 6), legend_position='top') +
                pn.labs(title=f'Fragment Length Distribution - {sample_id}',
                    subtitle=f'Sampled {self.actual_sampled_count:,} reads from {self.min_length}-{self.max_length}bp range',
                    x='Fragment Length (bp)',
                    y='Relative Frequency',
                    color='Dataset') +
                pn.scale_x_continuous(limits=(50, 500))
        )
        
        # Save plot
        plot_file = os.path.join(outpath, f'{sample_id}_fragment_length_distribution.png')
        plot.save(plot_file, width=12, height=6, dpi=300)
        print(f"Histogram saved as {plot_file}")


    def calculate_mean_fragment_length(self):
        """Calculate the mean fragment length from sampled reads."""
        if not self.sampled_lengths:
            return None
        
        mean_length = sum(self.sampled_lengths) / len(self.sampled_lengths)
        print(f"Calculated mean fragment length: {mean_length:.2f}")
        return mean_length


def main():
    args = parse_args()
    
    print(f"🔬 ScoreFrags.py - Fragmentation Score Calculator")
    print(f"📁 BAM file: {args.bam_file}")
    print(f"📊 Reference: {args.fragment_score_reference}")
    print(f"📏 Length range: {args.min_length}-{args.max_length}bp")
    print(f"🎯 Target reads: {args.read_count:,}")
    print("-" * 60)
    
    scorer = FragScorer(
        bam_file=args.bam_file,
        fragment_score_reference=args.fragment_score_reference,
        read_count=args.read_count,
        min_length=args.min_length,
        max_length=args.max_length,
        seed=args.seed
    )
    
    sample_name = args.sample_name or os.path.basename(args.bam_file).replace('.bam', '')
    mean_frag_length = scorer.calculate_mean_fragment_length()
    
    # Enhanced output with sampling information
    out_df = pd.DataFrame({
        "sample_name": [sample_name],
        "fragmentation_score": [scorer.FS],
        "mean_fragment_length": [mean_frag_length],
        "sampled_reads": [scorer.actual_sampled_count],
        "sampling_range": [f"{args.min_length}-{args.max_length}bp"],
        "requested_reads": [args.read_count]
    })
    
    out_file = os.path.join(args.output, f"{sample_name}_FragmentScore.tsv")
    out_df.to_csv(out_file, sep="\t", index=False)
    scorer.make_histogram(sample_id=sample_name, outpath=args.output)

    print(f"\n✅ Results:")
    print(f"   Fragmentation Score: {scorer.FS:.4f}")
    print(f"   Mean Fragment Length: {mean_frag_length:.2f}bp")
    print(f"   Sampled Reads: {scorer.actual_sampled_count:,}")
    print(f"   Results saved to: {out_file}")

if __name__ == "__main__":
    main()
