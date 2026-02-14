#!/usr/bin/env python3

import pandas as pd
import numpy as np

def load_hotspots(hotspot_file):
    """Load hotspot positions from text file (chr:position format)"""
    
    hotspot_positions = set()
    
    with open(hotspot_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):  # Skip empty lines and comments
                try:
                    chromosome, position = line.split(':')
                    hotspot_positions.add((chromosome, int(position)))
                except ValueError:
                    print(f"Warning: Could not parse line: {line}")
    
    print(f"Loaded {len(hotspot_positions)} hotspot positions")
    return hotspot_positions

def identify_artifacts(mutation_rates_df, hotspot_positions, threshold_sd):
    """Identify artifact positions using statistical threshold"""
    
    # Calculate statistics
    rates = mutation_rates_df['mutation_rate']
    mean_rate = rates.mean()
    std_rate = rates.std()
    threshold = mean_rate + (threshold_sd * std_rate)
    
    print(f"Mutation rate statistics:")
    print(f"  Mean: {mean_rate:.6f}")
    print(f"  Std:  {std_rate:.6f}")
    print(f"  Threshold ({threshold_sd}SD): {threshold:.6f}")
    
    # Find high-mutation positions
    high_mutation = mutation_rates_df[rates > threshold].copy()
    print(f"Positions above {threshold_sd}SD: {len(high_mutation)}")
    
    # Remove hotspot positions
    initial_count = len(high_mutation)
    artifact_positions = []
    
    for _, row in high_mutation.iterrows():
        position_key = (row['chromosome'], row['position'])
        if position_key not in hotspot_positions:
            artifact_positions.append(f"{row['chromosome']}:{row['position']}")
    
    final_count = len(artifact_positions)
    hotspot_removed = initial_count - final_count
    
    print(f"Removed {hotspot_removed} hotspot positions")
    print(f"Final blacklist positions: {final_count}")
    
    return artifact_positions

def main():
    """Generate artifact blacklist"""
    
    # Load inputs
    mutation_rates = pd.read_csv(snakemake.input.aggregated, sep='\t')
    hotspot_positions = load_hotspots(snakemake.params.hotspot_manifest)
    threshold_sd = snakemake.params.threshold_sd
    
    print(f"Processing {len(mutation_rates)} positions with {threshold_sd}SD threshold")
    
    # Identify artifacts
    blacklist_positions = identify_artifacts(mutation_rates, hotspot_positions, threshold_sd)
    
    # Write simple chr:position format
    with open(snakemake.output.blacklist, 'w') as f:
        for position in blacklist_positions:
            f.write(f"{position}\n")
    
    print(f"Wrote blacklist to: {snakemake.output.blacklist}")
    if blacklist_positions:
        print(f"Sample positions: {blacklist_positions[:5]}")

if __name__ == "__main__":
    main()
