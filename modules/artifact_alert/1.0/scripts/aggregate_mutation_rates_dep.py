#!/usr/bin/env python3

import pandas as pd
import numpy as np
from typing import List, Dict, Tuple

class PositionStats:
    """Calculate running mean and std deviation for a position using Welford's algorithm"""
    
    def __init__(self):
        self.n = 0
        self.mean = 0.0
        self.M2 = 0.0
    
    def add_sample_to_stats(self, value):
        """Add a new sample value to the running statistics"""
        self.n += 1
        delta = value - self.mean
        self.mean += delta / self.n
        delta2 = value - self.mean
        self.M2 += delta * delta2
    
    def get_std(self):
        """Return sample standard deviation"""
        if self.n < 2:
            return 0.0
        return np.sqrt(self.M2 / (self.n - 1))
    
    def get_stats(self):
        """Return (mean, std, count)"""
        return self.mean, self.get_std(), self.n


def detect_mutation_types(file_path: str) -> List[str]:
    """Detect mutation type columns from first input file."""
    df = pd.read_csv(file_path, sep='\t', nrows=1)
    return [col for col in df.columns if col not in ['chromosome', 'position', 'ref_base']]


def process_row(row, mutation_types: List[str], all_positions_dict: Dict[Tuple, PositionStats]) -> None:
    """Process a single row and update position statistics for all mutation types."""
    chr_pos_key = (row['chromosome'], row['position'], row['ref_base'])
    
    for mut_type in mutation_types:
        key = (*chr_pos_key, mut_type)

        if key not in all_positions_dict:
            all_positions_dict[key] = PositionStats()

        all_positions_dict[key].add_sample_to_stats(row[mut_type])


def aggregate_samples(input_files: List[str], mutation_types: List[str], log_file) -> Dict[Tuple, PositionStats]:
    """
    Process all sample files and aggregate statistics across samples.
    Uses pandas apply() for fastest performance.
    
    Returns:
        Dictionary mapping (chr, pos, ref, mut_type) to PositionStats object
    """
    all_positions_dict = {}
    
    for i, file_path in enumerate(input_files, 1):
        if i % 10 == 0 or i == 1:
            log_progress(f"Processing file {i}/{len(input_files)}...", log_file, flush=True)
        
        # Read file
        df = pd.read_csv(file_path, sep='\t')
        
        # Use apply, slightly faster than using itertuples(), a lot faster than iterrows()
        df.apply(process_row, axis=1, mutation_types=mutation_types, all_positions_dict=all_positions_dict)
    
    return all_positions_dict


def build_output_dataframe(all_positions_dict: Dict[Tuple, PositionStats]) -> pd.DataFrame:
    """Convert position statistics to output DataFrame."""
    results = []
    
    # Group by position
    position_data = {}
    for (chromosome, position, ref_base, mut_type), stats in all_positions_dict.items():
        pos_key = (chromosome, position, ref_base)
        
        if pos_key not in position_data:
            position_data[pos_key] = {'n_samples': stats.n}
        
        mean, std, n = stats.get_stats()
        position_data[pos_key][f"{mut_type}_mean"] = round(mean, 8)
        position_data[pos_key][f"{mut_type}_std"] = round(std, 8)
    
    # Convert to list of records
    for (chromosome, position, ref_base), data in position_data.items():
        results.append({
            'chromosome': chromosome,
            'position': position,
            'ref_base': ref_base,
            **data
        })
    
    # Create and sort DataFrame
    df_out = pd.DataFrame(results)
    return df_out.sort_values(['chromosome', 'position'])


def log_progress(message: str, log_file, flush: bool = False):
    """Helper to print and optionally flush log."""
    print(message, file=log_file)
    if flush:
        log_file.flush()


def aggregate_mutation_rates(input_files: List[str], output_file: str, log_file) -> None:
    """Process all samples and write aggregated results."""
    log_progress(f"Starting mutation rate aggregation", log_file, flush=True)
    log_progress(f"Number of input files: {len(input_files)}", log_file)
    log_progress(f"Output file: {output_file}\n", log_file)
    
    # Detect mutation types from first file
    mutation_types = detect_mutation_types(input_files[0])
    log_progress(f"Detected mutation types: {mutation_types}\n", log_file, flush=True)
    
    # Process all samples and aggregate
    all_positions_dict = aggregate_samples(input_files, mutation_types, log_file)
    
    log_progress(f"\nCalculating final statistics...", log_file, flush=True)
    
    # Build output DataFrame
    df_out = build_output_dataframe(all_positions_dict)
    df_out.to_csv(output_file, sep='\t', index=False)
    
    # Log summary statistics
    n_sample_dist = df_out['n_samples'].describe()
    log_progress(f"\nResults:", log_file)
    log_progress(f"  Total positions: {len(df_out)}", log_file)
    log_progress(f"\nSample coverage distribution:", log_file)
    log_progress(f"  Min: {int(n_sample_dist['min'])}, Max: {int(n_sample_dist['max'])}", log_file)
    log_progress(f"  Median: {int(n_sample_dist['50%'])}, Mean: {n_sample_dist['mean']:.1f}", log_file)
    log_progress(f"\nMean background rates:", log_file)
    for mut_type in mutation_types:
        mean_col = f"{mut_type}_mean"
        if mean_col in df_out.columns:
            log_progress(f"  {mut_type}: {df_out[mean_col].mean():.8f}", log_file)
    log_progress(f"\nCompleted successfully", log_file)


def main():
    """Main entry point for Snakemake workflow"""
    log_file = open(snakemake.log[0], 'w')
    aggregate_mutation_rates(
        input_files=snakemake.input.mutation_files,
        output_file=snakemake.output.aggregated,
        log_file=log_file
    )
    log_file.close()


if __name__ == "__main__":
    main()
