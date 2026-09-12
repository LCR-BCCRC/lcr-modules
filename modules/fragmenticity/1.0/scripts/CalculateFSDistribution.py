"""Script to read in length profiles of mutated and healthy cfDNA reads
and calculate their the ratio of the probability density for each fragment
length. That value is log2 transformed, and then a file is output with
900 lines, with each line represnting the fragmentation score for each length
of cfDNA fragment from 1-900.

Pseudocounts are added to handle zero cases. Fragment lengths that never appear in either population get a fragmentation score of 0.
The script also outputs a file with the number of reads for each length.

Both the mutated and healthy reads will be subsampled by the same amount, that amount is set
to a minimum of 10000 reads from both, but can be increased if more reads are available.

This whole process is then bootstrapped 1000 times and the mean is taken for the FS of each length.

This script is designed to take in the outputs from ProfileMutatedReads.py and ProfileHealthyReads.py,
those outputs are formatted such that the first column is the length of the read, and the second column
is the number of reads of that length.

"""
import pandas as pd
import argparse
import numpy as np
import plotnine as pn
import os

COLORS = [
"#1D252D",
"#EAAA00",
"#279989",
]

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate fragmentation score distribution from mutated and healthy cfDNA reads.")
    parser.add_argument("--mutated", type=str, nargs = "+",required=True, help="List of paths to the files containing length profile of mutated reads.")
    parser.add_argument("--healthy", type=str, nargs = "+", required=True, help="List of paths to the file containing length profile of healthy reads.")
    parser.add_argument("--output", type=str, required=True, help="Path to the output file for fragmentation score atlas.")
    parser.add_argument("--max_reads", type=int, default=100000, help="Max number of reads to subsample from the mutated and healthy read profiles (default: 100000).")
    return parser.parse_args()


class FragmentAtlasCalculator(object):
    def __init__(self, mutated_files, healthy_files, output_file, max_reads=100000):
        self.mutated_files = mutated_files
        self.healthy_files = healthy_files
        self.output_file = output_file
        self.reads = max_reads
        # outputs
        self.fragmentation_scores = None

    def read_and_aggregate_length_profiles(self, files, type):
        """Read and add together length profiles from provided files."""
        df_list = []
        skipped_files = 0
        for file in files:

            if os.path.exists(file) and os.path.getsize(file) > 1:
                df = pd.read_csv(file, sep="\t")
                df.columns = ["length", "count"]
                # if empty, skip
                if df.empty:
                    print(f"Warning: {file} is empty, skipping.")
                    skipped_files += 1
                    continue
                # if less than 100 reads, skip
                if df["count"].sum() < 100:
                    print(f"Warning: {file} has less than 100 reads, skipping.")
                    skipped_files += 1
                    continue
            else:
                print(f"Warning: {file} does not exist or is empty, skipping.")
                skipped_files += 1
                continue

            df_list.append(df)

            # Handle case where no valid files were found
        if not df_list:
            raise ValueError("No valid input files found. All files were empty or unreadable.")
        print("\n")
        print(f"Processing {type} files:")
        print(f"Total input files processed: {len(files)}")
        print(f"Total valid files processed: {len(df_list)}")
        print(f"Total files skipped due to being empty or unreadable: {skipped_files}")
        print("\n")
        combined_df = pd.concat(df_list).groupby("length").sum().reset_index()
        combined_df["length"] = combined_df["length"].astype(int)
        return combined_df

    def calculate_fragmentation_scores(self, sample_size=10000, bootstrap_iterations=1000, pseudo_count= 20):
        """
        calculate the FS atlas for each fragment size from 0-900
        by subsampling the mutated and healthy reads to the same size.


        """
        # get all read profiles
        mutated_read_profile = self.read_and_aggregate_length_profiles(self.mutated_files, "mutated")
        healthy_read_profile = self.read_and_aggregate_length_profiles(self.healthy_files, "healthy")
        # figure out how many reads to subsample
        subsample_size = self.calc_subsample_size(mutated_read_profile["count"].sum(), healthy_read_profile["count"].sum())
        # subsample the mutated reads, bootstrapping the process
        mutated_bootstrapped_reads = self.bootstrap_sample_reads(mutated_read_profile, bootstrap_iterations, subsample_size, pseudo_count)
        # subsample the healthy reads, bootstrapping the process
        healthy_bootstrapped_reads = self.bootstrap_sample_reads(healthy_read_profile, bootstrap_iterations, subsample_size, pseudo_count)
        # calculate the FS scores atlas
        fs_scores = self.calc_log_FS(mutated_bootstrapped_reads, healthy_bootstrapped_reads, pseudo_count=pseudo_count)

        self.frag_plot = self.plot_sampled_read_profiles(mutated_bootstrapped_reads, healthy_bootstrapped_reads)
        self.fragmentation_scores = fs_scores.copy()

    def plot_sampled_read_profiles(self, mutated_read_profile, healthy_read_profile):

        mut_df = mutated_read_profile.copy()
        mut_df["type"] = "mutated"
        healthy_df = healthy_read_profile.copy()
        healthy_df["type"] = "healthy"

        plot = (pn.ggplot(pd.concat([mut_df, healthy_df]), pn.aes(x="length", y="probability_density", fill="type", color="type") ) +
                pn.geom_line() +
                pn.geom_area(alpha=0.5, position = "identity") +
                pn.theme_bw() +
                pn.theme(legend_position="top") +
                pn.scale_fill_manual(values=COLORS) +
                pn.scale_color_manual(values=COLORS) +
                pn.theme(figure_size=(10, 6),
                        strip_background= pn.element_rect(fill="white")
                ) +
                pn.labs(
                    title="Subsampled Read Length Profiles",
                    x="Fragment Length (bp)",
                    y="Fragment Length Frequency",
                    fill="Read Type",
                    color="Read Type"
                )

        )
        return plot


    def calc_log_FS(self, mut_bootstrapped_profile, healthy_bootstrapped_profile, pseudo_count=20):
        """
        Calculate FS, setting to 0 for lengths that never appeared in either population.
        """
        FS_scores = {}
        
        for length in range(1, 901):
            mut_row = mut_bootstrapped_profile[mut_bootstrapped_profile["length"] == length]
            healthy_row = healthy_bootstrapped_profile[healthy_bootstrapped_profile["length"] == length]
            
            mut_pd = mut_row["probability_density"].iloc[0]
            healthy_pd = healthy_row["probability_density"].iloc[0]
            mut_count = mut_row["count"].iloc[0]
            healthy_count = healthy_row["count"].iloc[0]
            
            # If both counts equal pseudocount exactly, the length never appeared naturally
            if mut_count == pseudo_count and healthy_count == pseudo_count:
                FS_scores[length] = 0
            else:
                fs_score = np.log2(mut_pd / healthy_pd)
                FS_scores[length] = fs_score
        
        fs_df = pd.DataFrame(list(FS_scores.items()), columns=["length", "fragmentation_score"])
        fs_df.set_index("length", inplace=True)
        return fs_df


    def bootstrap_sample_reads(self, read_profile, bootstrap_size, subsample_size, pseudo_count=20):
        """
        Bootstrap sample the read profile and return average probability densities.
        Ensures ALL fragment lengths 1-900 are represented in each bootstrap sample.

        Args:
            read_profile (pd.DataFrame): DataFrame with columns "length" and "count"
            bootstrap_size (int): Number of bootstrap iterations
            subsample_size (int): Number of reads to sample in each bootstrap
            pseudo_count (int): Pseudocount to add to all lengths
        Returns:
            pd.DataFrame: Average probability densities across all bootstrap iterations
        """
        # Create template with all possible fragment lengths 1-900
        all_lengths_template = pd.DataFrame({'length': range(1, 901), 'count': 0})

        # Expand read lengths from counts to individual reads
        read_length_list = read_profile["length"].repeat(read_profile["count"]).reset_index(drop=True)

        bootstrap_results = []

        for i in range(bootstrap_size):
            # Sample with replacement
            if len(read_length_list) >= subsample_size:
                sampled_lengths = read_length_list.sample(n=subsample_size, replace=True)
            else:
                sampled_lengths = read_length_list.sample(n=len(read_length_list), replace=True)
            
            # Count occurrences of each length
            length_counts = sampled_lengths.value_counts().reset_index()
            length_counts.columns = ['length', 'count']
            
            # Merge with template to ensure all lengths 1-900 are represented
            # This fills missing lengths with count=0
            complete_counts = all_lengths_template.merge(length_counts, on='length', how='left', suffixes=('', '_sampled'))
            complete_counts['count'] = complete_counts['count_sampled'].fillna(0)
            complete_counts = complete_counts[['length', 'count']]
            
            # Add pseudocount to ALL lengths (including those with count=0)
            complete_counts['count'] += pseudo_count
            
            # Calculate probability density
            complete_counts['probability_density'] = complete_counts['count'] / complete_counts['count'].sum()
            
            bootstrap_results.append(complete_counts)

        # Average across all bootstrap iterations
        all_lengths = pd.concat(bootstrap_results, ignore_index=True)
        averaged_results = all_lengths.groupby('length').agg({
            'count': 'mean',
            'probability_density': 'mean'
        }).reset_index()

        return averaged_results

    def calc_subsample_size(self, mutated_read_count, healthy_read_count):

        # round down mutated read count to the nearest 1000
        max_mutated_reads = (mutated_read_count // 1000) * 1000

        # check to make sure there is enough healthy reads to subsample from
        if max_mutated_reads > healthy_read_count:
            print(f"Mutated read count ({mutated_read_count}) is greater than healthy read count ({healthy_read_count})")

        subsample_size = min(max_mutated_reads, healthy_read_count, self.reads)

        print(f"\n\n###########Subsampling to {subsample_size} reads from both mutated and healthy read profiles.##########")
        print(f"Mutated read count: {mutated_read_count}, Healthy read count: {healthy_read_count}\n\n")

        return subsample_size

def main():
    args = parse_args()
    print(f"Calculating fragmentation scores for mutated reads from: {args.mutated}")
    print(f"Calculating fragmentation scores for healthy reads from: {args.healthy}")
    calculator = FragmentAtlasCalculator(args.mutated, args.healthy, args.output, max_reads=args.max_reads)
    calculator.calculate_fragmentation_scores()
    # save the fragmentation scores to a file
    calculator.fragmentation_scores.to_csv(args.output, sep="\t", header=True, index=True)
    # save the plot to a file
    plot_file = args.output.replace(".tsv", "_read_profiles.pdf")
    calculator.frag_plot.save(plot_file, width=10, height=6, dpi=300)
    print(f"Fragmentation scores saved to {args.output}")
    print(f"Read profiles plot saved to {plot_file}")

if __name__ == "__main__":
    main()