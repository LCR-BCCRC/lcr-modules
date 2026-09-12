"""Script for regional fragmentation scoring using bed file targets

This script extends the original ScoreFrags.py to calculate fragmentation scores
per genomic region defined in a bed file. It provides both genome-wide and
adaptive regional scoring approaches for enhanced ctDNA detection sensitivity.

Key features:
- Collapses adjacent bed regions by name (using min/max coordinates)
- Calculates FS per region with configurable minimum coverage
- Provides adaptive scoring (mean of top-K and max regional scores)
- Dual output: main results + detailed regional scores

Example usage:
python ScoreFragsRegional.py \
    --bam_file sample.bam \
    --bed_file targets.bed \
    --fragment_score_reference FragmentLengthReference.tsv \
    --sample_name SAMPLE_001 \
    --top_k_regions 10 \
    --min_reads_per_region 1000

Written by: Kurt Yakimovich (adapted from original ScoreFrags.py)
"""

import pandas as pd
import os
import pysam
import argparse
import random
import plotnine as pn
from collections import Counter, defaultdict
import warnings
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate regional fragmentation scores from BAM files using bed file targets.")
    parser.add_argument("--bam_file", required=True, help="Path to the input BAM file.")
    parser.add_argument("--bed_file", required=True, help="Path to bed file with target regions (chr, start, end, name format).")
    parser.add_argument("--fragment_score_reference", type=str, required=True, help="Path to the reference file with scores for fragment lengths 1-900, one column no headers")
    parser.add_argument("--sample_name", type=str, default=None, help="Sample name to use in the output file, default is the BAM file name without extension.")
    parser.add_argument("--output", default="./", help="Output directory")
    parser.add_argument("--read_count", type=int, default=100000, help="Number of reads to sample from the BAM file for genome-wide score, default 100k reads.")
    parser.add_argument("--min_length", type=int, default=50, help="Minimum length of reads to be considered, default 50.")
    parser.add_argument("--max_length", type=int, default=900, help="Maximum length of reads to be considered, default 900.")
    parser.add_argument("--top_k_regions", type=int, default=10, help="Number of top-scoring regions to use for adaptive FS calculation, default 10.")
    parser.add_argument("--min_reads_per_region", type=int, default=1000, help="Minimum number of reads required per region to calculate FS, default 1000.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility, default 42.")

    return parser.parse_args()

class RegionalFragScorer(object):
    """Enhanced fragmentation scorer with regional analysis capabilities."""
    
    def __init__(self, bam_file, bed_file, fragment_score_reference=None, read_count=100000, 
                 min_length=50, max_length=900, top_k_regions=10, min_reads_per_region=1000, seed=42):
        self.bam_file = bam_file
        self.bed_file = bed_file
        self.fragment_score_reference = fragment_score_reference
        self.read_count = read_count
        self.min_length = min_length
        self.max_length = max_length
        self.top_k_regions = top_k_regions
        self.min_reads_per_region = min_reads_per_region
        self.seed = seed

        # Initialize data storage
        self.collapsed_regions = {}
        self.regional_reads = {}
        self.regional_fs_scores = pd.DataFrame()
        self.sampled_lengths = []
        self.all_valid_lengths = []
        self.actual_sampled_count = 0
        
        # Calculate scores
        self.scores_dict = self.read_score_reference()
        self.collapsed_regions = self.parse_and_collapse_bed()
        self.genome_wide_FS = self.calculate_genome_wide_fragmentation_score()
        self.regional_reads = self.sample_reads_by_region()
        self.regional_fs_scores = self.calculate_regional_fs_scores()
        self.adaptive_FS, self.max_regional_FS = self.calculate_adaptive_scores()

    def read_score_reference(self):
        """Reads the fragment score reference file and creates a dictionary."""
        scores_dict = {}
        raw_ref = pd.read_csv(self.fragment_score_reference, sep="\t")
        
        if raw_ref.shape[1] == 1 and "0" in raw_ref.columns:
            scores_dict = {i + 1: score for i, score in enumerate(raw_ref.iloc[:, 0])}
        elif "fragmentation_score" in raw_ref.columns and "length" in raw_ref.columns:
            scores_dict = raw_ref.set_index("length")["fragmentation_score"].to_dict()
        else:
            raise ValueError("""Fragment score reference file must have one column with no 
                                headers or two columns with 'length' and 'fragmentation_score' headers.""")
        
        print(f"📊 Loaded scores for fragment lengths 1-{len(scores_dict)} from reference file")
        return scores_dict

    def parse_and_collapse_bed(self):
        """Parse bed file and collapse regions by name using min/max coordinates."""
        print(f"📁 Parsing bed file: {self.bed_file}")
        
        bed_df = pd.read_csv(self.bed_file, sep="\t", header=None, 
                           names=["chr", "start", "end", "name"])
        
        if bed_df.shape[1] < 4:
            raise ValueError("Bed file must have at least 4 columns: chr, start, end, name")
        
        # Group by name and collapse coordinates
        collapsed = {}
        region_groups = bed_df.groupby("name")
        
        for name, group in region_groups:
            # Use the first chromosome (assume all regions with same name are on same chr)
            chr_name = group["chr"].iloc[0]
            min_start = group["start"].min()
            max_end = group["end"].max()
            
            collapsed[name] = {
                "chr": chr_name,
                "start": min_start,
                "end": max_end,
                "original_regions": len(group)
            }
        
        print(f"🔗 Collapsed {len(bed_df)} bed regions into {len(collapsed)} named regions")
        for name, region in list(collapsed.items())[:5]:  # Show first 5
            orig_count = region["original_regions"]
            print(f"   {name}: {region['chr']}:{region['start']}-{region['end']} (from {orig_count} regions)")
        
        if len(collapsed) > 5:
            print(f"   ... and {len(collapsed) - 5} more regions")
            
        return collapsed

    def sample_reads_by_region(self):
        """Sample reads from each collapsed region with filtering."""
        print(f"🎯 Sampling reads by region (min {self.min_reads_per_region} reads per region)")
        
        regional_reads = {}
        processed_UMIs = set()
        
        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            for region_name, region_info in self.collapsed_regions.items():
                region_reads = []
                region_UMIs = set()
                
                # Fetch reads in this region
                try:
                    for read in bam.fetch(region_info["chr"], region_info["start"], region_info["end"]):
                        # Apply same filtering as original script
                        if read.is_unmapped or read.is_secondary or read.is_supplementary:
                            continue
                        if read.mapping_quality < 20:
                            continue

                        # UMI-based duplicate avoidance
                        try:
                            um = read.get_tag("MI")
                            if um in region_UMIs:
                                continue
                            region_UMIs.add(um)
                        except KeyError:
                            continue

                        template_length = abs(read.template_length)
                        
                        # Check if read is in our length range
                        if self.min_length <= template_length <= self.max_length:
                            region_reads.append(template_length)
                
                except Exception as e:
                    print(f"⚠️  Warning: Could not fetch reads for region {region_name}: {e}")
                    continue
                
                # Store reads for this region
                regional_reads[region_name] = region_reads
                
                # Progress reporting
                if len(region_reads) >= self.min_reads_per_region:
                    status = "✅"
                else:
                    status = "❌"
                print(f"   {status} {region_name}: {len(region_reads)} reads")
        
        # Summary
        valid_regions = sum(1 for reads in regional_reads.values() if len(reads) >= self.min_reads_per_region)
        print(f"📈 Found {valid_regions}/{len(self.collapsed_regions)} regions with ≥{self.min_reads_per_region} reads")
        
        return regional_reads

    def calculate_regional_fs_scores(self):
        """Calculate fragmentation scores for each region."""
        print(f"🧮 Calculating FS scores for regions with ≥{self.min_reads_per_region} reads")
        
        regional_scores = []
        
        for region_name, read_lengths in self.regional_reads.items():
            if len(read_lengths) < self.min_reads_per_region:
                continue
                
            # Calculate FS for this region
            fs_scores = [self.scores_dict[length] for length in read_lengths]
            mean_fs = np.mean(fs_scores)
            
            regional_scores.append({
                "region_name": region_name,
                "read_count": len(read_lengths),
                "fragmentation_score": mean_fs,
                "mean_fragment_length": np.mean(read_lengths)
            })
        
        # Convert to DataFrame and sort by FS (highest first)
        regional_df = pd.DataFrame(regional_scores)
        if not regional_df.empty:
            regional_df = regional_df.sort_values("fragmentation_score", ascending=False)
        
        print(f"✅ Calculated FS scores for {len(regional_df)} regions")
        if not regional_df.empty:
            print(f"   Top region: {regional_df.iloc[0]['region_name']} (FS: {regional_df.iloc[0]['fragmentation_score']:.4f})")
            print(f"   Bottom region: {regional_df.iloc[-1]['region_name']} (FS: {regional_df.iloc[-1]['fragmentation_score']:.4f})")
        
        return regional_df

    def calculate_adaptive_scores(self):
        """Calculate adaptive FS (mean of top-K) and max regional FS."""
        if self.regional_fs_scores.empty:
            print("⚠️  No valid regional scores - cannot calculate adaptive scores")
            return None, None
        
        # Get top K regions
        top_k = min(self.top_k_regions, len(self.regional_fs_scores))
        top_regions = self.regional_fs_scores.head(top_k)
        
        # Calculate adaptive scores
        adaptive_fs = top_regions["fragmentation_score"].mean()
        max_regional_fs = self.regional_fs_scores["fragmentation_score"].max()
        
        print(f"🎯 Adaptive scoring results:")
        print(f"   Top {top_k} regions used for adaptive_FS")
        print(f"   Adaptive FS (mean of top {top_k}): {adaptive_fs:.4f}")
        print(f"   Max regional FS: {max_regional_fs:.4f}")
        print(f"   Top {top_k} regions: {', '.join(top_regions['region_name'].tolist())}")
        
        return adaptive_fs, max_regional_fs

    def sample_fragment_lengths_genome_wide(self):
        """Sample fragment lengths genome-wide (original approach)."""
        random.seed(self.seed)
        
        valid_lengths_for_sampling = []
        all_valid_lengths = []
        total_processed = 0
        processed_UMIs = set()

        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            for read in bam.fetch():
                total_processed += 1
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < 20:
                    continue

                try:
                    um = read.get_tag("MI")
                    if um in processed_UMIs:
                        continue
                    processed_UMIs.add(um)
                except KeyError:
                    continue

                template_length = abs(read.template_length)
                
                if 50 <= template_length <= 900:
                    all_valid_lengths.append(template_length)
                
                if self.min_length <= template_length <= self.max_length:
                    valid_lengths_for_sampling.append(template_length)
        
        self.all_valid_lengths = all_valid_lengths
        
        print(f"🌍 Genome-wide sampling: Found {len(valid_lengths_for_sampling)} valid reads in range {self.min_length}-{self.max_length}bp")
        
        if len(valid_lengths_for_sampling) < self.read_count:
            warning_msg = (f"⚠️  WARNING: Only found {len(valid_lengths_for_sampling)} valid reads "
                          f"for genome-wide sampling, but {self.read_count} were requested.")
            print(warning_msg)
            warnings.warn(warning_msg, UserWarning)
            
            self.actual_sampled_count = len(valid_lengths_for_sampling)
            self.sampled_lengths = valid_lengths_for_sampling
            return valid_lengths_for_sampling
        
        sampled_lengths = random.sample(valid_lengths_for_sampling, self.read_count)
        self.actual_sampled_count = len(sampled_lengths)
        self.sampled_lengths = sampled_lengths
        
        print(f"✅ Successfully sampled {len(sampled_lengths)} reads for genome-wide FS")
        return sampled_lengths

    def calculate_genome_wide_fragmentation_score(self):
        """Calculate genome-wide fragmentation score (original approach)."""
        read_lengths = self.sample_fragment_lengths_genome_wide()
        
        fs_scores = [self.scores_dict[length] for length in read_lengths]
        fragmentation_score = np.mean(fs_scores)
        
        print(f"🌍 Genome-wide FS: {fragmentation_score:.4f} (from {len(fs_scores)} reads)")
        return fragmentation_score

    def calculate_mean_fragment_length(self):
        """Calculate mean fragment length from genome-wide sampled reads."""
        if not self.sampled_lengths:
            return None
        return np.mean(self.sampled_lengths)

    def make_histogram(self, sample_id, outpath):
        """Create histogram showing fragment length distribution."""
        if not self.all_valid_lengths:
            print("⚠️  No data for histogram generation")
            return
            
        all_counts = pd.Series(self.all_valid_lengths).value_counts().reset_index()
        all_counts.columns = ['length', 'count']
        all_counts['relative_frequency'] = all_counts['count'] / all_counts['count'].sum()
        all_counts['read_type'] = 'All Valid Reads'
        
        if len(self.sampled_lengths) > 0:
            sampled_counts = pd.Series(self.sampled_lengths).value_counts().reset_index()
            sampled_counts.columns = ['length', 'count']
            sampled_counts['relative_frequency'] = sampled_counts['count'] / sampled_counts['count'].sum()
            sampled_counts['read_type'] = 'Genome-wide Sampled'
            
            plot_data = pd.concat([all_counts, sampled_counts], ignore_index=True)
        else:
            plot_data = all_counts
        
        plot = (pn.ggplot(plot_data, pn.aes(x='length', y='relative_frequency', 
                                        color='read_type', group='read_type')) +
                pn.geom_line(size=1, alpha=0.8) +
                pn.geom_vline(xintercept=self.min_length, color='red', linetype='dashed', alpha=0.8) +
                pn.geom_vline(xintercept=self.max_length, color='red', linetype='dashed', alpha=0.8) +
                pn.geom_vline(xintercept=167, color='orange', linetype='dotted', alpha=0.6) +
                pn.geom_vline(xintercept=333, color='purple', linetype='dotted', alpha=0.6) +
                pn.theme_bw() +
                pn.theme(figure_size=(12, 6), legend_position='top') +
                pn.labs(title=f'Fragment Length Distribution - {sample_id} (Regional Analysis)',
                    subtitle=f'Genome-wide: {self.actual_sampled_count:,} reads | Regional: {len(self.regional_fs_scores)} valid regions',
                    x='Fragment Length (bp)',
                    y='Relative Frequency',
                    color='Dataset') +
                pn.scale_x_continuous(limits=(50, 500))
        )
        
        plot_file = os.path.join(outpath, f'{sample_id}_regional_fragment_distribution.png')
        plot.save(plot_file, width=12, height=6, dpi=300)
        print(f"📊 Histogram saved as {plot_file}")

    def generate_outputs(self, sample_name, output_dir):
        """Generate main output file and regional scores file."""
        
        # Main output file (enhanced with regional scores)
        mean_frag_length = self.calculate_mean_fragment_length()
        
        main_output = pd.DataFrame({
            "sample_name": [sample_name],
            "fragmentation_score": [self.genome_wide_FS],
            "adaptive_FS": [self.adaptive_FS],
            "max_regional_FS": [self.max_regional_FS],
            "mean_fragment_length": [mean_frag_length],
            "sampled_reads": [self.actual_sampled_count],
            "valid_regions": [len(self.regional_fs_scores)],
            "top_k_regions": [self.top_k_regions],
            "sampling_range": [f"{self.min_length}-{self.max_length}bp"],
            "requested_reads": [self.read_count]
        })
        
        main_file = os.path.join(output_dir, f"{sample_name}_FragmentScore.tsv")
        main_output.to_csv(main_file, sep="\t", index=False)
        
        # Regional scores output file
        if not self.regional_fs_scores.empty:
            regional_file = os.path.join(output_dir, f"{sample_name}_RegionalFragmentScores.tsv")
            self.regional_fs_scores.to_csv(regional_file, sep="\t", index=False)
            print(f"📋 Regional scores saved to: {regional_file}")
        else:
            print("⚠️  No regional scores to save")
        
        print(f"💾 Main results saved to: {main_file}")
        return main_file

def main():
    args = parse_args()
    
    print(f"🔬 ScoreFragsRegional.py - Regional Fragmentation Score Calculator")
    print(f"📁 BAM file: {args.bam_file}")
    print(f"🎯 Bed file: {args.bed_file}")
    print(f"📊 Reference: {args.fragment_score_reference}")
    print(f"📏 Length range: {args.min_length}-{args.max_length}bp")
    print(f"🎯 Genome-wide reads: {args.read_count:,}")
    print(f"🔝 Top K regions: {args.top_k_regions}")
    print(f"🎛️  Min reads/region: {args.min_reads_per_region}")
    print("-" * 80)
    
    # Initialize scorer
    scorer = RegionalFragScorer(
        bam_file=args.bam_file,
        bed_file=args.bed_file,
        fragment_score_reference=args.fragment_score_reference,
        read_count=args.read_count,
        min_length=args.min_length,
        max_length=args.max_length,
        top_k_regions=args.top_k_regions,
        min_reads_per_region=args.min_reads_per_region,
        seed=args.seed
    )
    
    # Generate outputs
    sample_name = args.sample_name or os.path.basename(args.bam_file).replace('.bam', '')
    main_file = scorer.generate_outputs(sample_name, args.output)
    scorer.make_histogram(sample_id=sample_name, outpath=args.output)

    # Summary output
    print(f"\n✅ Final Results for {sample_name}:")
    print(f"   Genome-wide FS: {scorer.genome_wide_FS:.4f}")
    print(f"   Adaptive FS (top {args.top_k_regions}): {scorer.adaptive_FS:.4f}")
    print(f"   Max Regional FS: {scorer.max_regional_FS:.4f}")
    print(f"   Valid Regions: {len(scorer.regional_fs_scores)}")
    print(f"   Genome-wide Reads: {scorer.actual_sampled_count:,}")

if __name__ == "__main__":
    main()
