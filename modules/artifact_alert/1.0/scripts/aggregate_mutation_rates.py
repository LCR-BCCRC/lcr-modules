#!/usr/bin/env python3
"""
Aggregate Mutation Rates across multiple samples.

This script aggregates mutation rate data across multiple samples and chromosomes,
supporting incremental updates to an existing mutation rate index.

Features:
---------
- Parallel processing by chromosome
- Incremental updates: merges new samples with existing aggregated data
- Outputs bgzip-compressed, tabix-indexed files for fast genomic lookups
- Tracks processed samples to avoid duplicates
- Uses Welford's algorithm for aggregating variance

Input:
------
- Per-sample, per-chromosome mutation rate files (.tsv.gz format)
- Optional: existing aggregated mutation rate file to update

Output:
-------
- Aggregated mutation rates with mean, std, and M2 (sum of squares) for each mutation type
- bgzip-compressed TSV with tabix index
- Updated sample tracker file

Columns: chromosome, position, ref_base, A_mean, A_std, A_M2, T_mean, T_std, T_M2, C_mean, 
    C_std, C_M2, G_mean, G_std, G_M2, INS_mean, INS_std, INS_M2, DEL_mean, DEL_std, DEL_M2, n_samples

Author: Kurt Yakimovich
"""

import pandas as pd
import numpy as np
import pysam
import os
import shutil
import gzip
from typing import List, Dict, Tuple, Optional
from multiprocessing import Pool
from dataclasses import dataclass
from pathlib import Path
import argparse
import datetime
import sys

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Aggregate mutation rates across samples with incremental updates")
    
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='Input mutation rate files (per-sample, per-chromosome .tsv.gz files)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output aggregated mutation rates file (.tsv.gz, will be bgzip compressed and tabix indexed)'
    )
    parser.add_argument(
        '-s', '--sample-ids',
        nargs='+',
        required=True,
        help='Sample IDs being processed (used for tracking)'
    )
    parser.add_argument(
        '-c', '--chromosomes',
        nargs='+',
        required=True,
        help='Chromosome names (e.g., chr1 chr2 ... chrX chrY)'
    )
    parser.add_argument(
        '-t', '--sample-tracker',
        required=True,
        help='Path to sample tracker file (tracks processed samples, one per line)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of parallel threads for chromosome processing (default: 1)'
    )
    parser.add_argument(
        '--log',
        default=None,
        help='Path to log file (default: write to stdout)'
    )
    
    return parser.parse_args()

ALL_cols = ['chromosome', 'position', 'ref_base', 'A_mean', 'A_std', 'A_M2', 
            'T_mean', 'T_std', 'T_M2', 'C_mean', 'C_std', 'C_M2',
            'G_mean', 'G_std', 'G_M2', 'INS_mean', 'INS_std', 'INS_M2',
            'DEL_mean', 'DEL_std', 'DEL_M2', 'n_samples']
MUTATION_TYPES = ['A', 'T', 'C', 'G', 'INS', 'DEL']

@dataclass
class MutationStats:
    """
    Statistics for a single mutation type using Welford's algorithm.
    Does not store n - that's shared at the position level.
    
    Attributes:
        mean : float
            Running mean value
        M2 : float
            Sum of squared differences from mean
    """
    
    def __init__(self, mean: float = 0.0, M2: float = 0.0):
        self.mean = mean
        self.M2 = M2
    
    def add_sample(self, value: float, n: int) -> None:
        """
        Add a new sample value using Welford's algorithm.
        
        Arguments:
            value : float
                The new sample value
            n : int
                Current sample count (after incrementing)
        """
        delta = value - self.mean
        self.mean += delta / n
        delta2 = value - self.mean
        self.M2 += delta * delta2
    
    def get_variance(self, n: int) -> float:
        """Return sample variance given n samples."""
        if n < 2:
            return 0.0
        return self.M2 / (n - 1)
    
    def get_std(self, n: int) -> float:
        """Return sample standard deviation given n samples."""
        return np.sqrt(self.get_variance(n))


class PositionData:
    """
    Holds all data for a single genomic position.
    
    Attributes:
        n : int
            Number of samples (shared across all mutation types)
        ref : str
            Reference base
        mutations : dict
            Maps mutation type -> MutationStats
    """
    
    def __init__(self, ref: str, n: int = 0):
        self.n = n
        self.ref = ref
        self.mutations: Dict[str, MutationStats] = {}
    
    def add_sample(self, mutation_values: Dict[str, float]) -> None:
        """
        Add one sample with values for all mutation types.
        
        Arguments:
            mutation_values : dict
                Maps mut_type -> value for this sample
        """
        self.n += 1  # Increment once per sample
        
        for mut_type, value in mutation_values.items():
            if mut_type not in self.mutations:
                self.mutations[mut_type] = MutationStats()
            
            self.mutations[mut_type].add_sample(value, self.n)

    def get_mutation_stats(self, mut_type: str) -> Optional[MutationStats]:
        """Get MutationStats for a specific mutation type."""
        return self.mutations.get(mut_type)

class ChromIndex:
    """
    Holds statistics for all positions in a chromosome.
    
    Attributes:
        chrom : str
            Chromosome name
        all_stats : dict
            Maps position (int) -> PositionData
    """
    
    def __init__(self, chrom: str):
        self.chrom = chrom
        self.all_stats: Dict[int, PositionData] = {}
    
    def add_sample(self, pos: int, ref: str, mutation_values: Dict[str, float]) -> None:
        """
        Add one sample's mutation values for a specific position.
        
        Arguments:
            pos : int
                Genomic position
            ref : str
                Reference base
            mutation_values : dict
                Maps mut_type -> value for this sample
        """
        if pos not in self.all_stats:
            self.all_stats[pos] = PositionData(ref=ref)
        
        self.all_stats[pos].add_sample(mutation_values)
    
    def get_stats(self, pos: int, mut_type: str) -> Optional[MutationStats]:
        """Get MutationStats for a specific position and mutation type."""
        if pos not in self.all_stats:
            return None
        return self.all_stats[pos].get_mutation_stats(mut_type)
    
    def get_all_positions_as_df(self) -> pd.DataFrame:
        """Return DataFrame with all positions and mutation stats for this chromosome."""
        sorted_positions = sorted(self.all_stats.keys())
        
        rows = []
        for pos in sorted_positions:
            pos_data = self.all_stats[pos]
            
            row = {
                'chromosome': self.chrom,
                'position': pos,
                'ref_base': pos_data.ref
            }
            
            # Add stats for each mutation type
            for mut in MUTATION_TYPES:
                mut_stats = pos_data.get_mutation_stats(mut)
                if mut_stats:
                    row[f'{mut}_mean'] = mut_stats.mean
                    row[f'{mut}_std'] = mut_stats.get_std(pos_data.n)
                    row[f'{mut}_M2'] = mut_stats.M2
            
            row['n_samples'] = pos_data.n
            rows.append(row)
        
        return pd.DataFrame(rows)

def natural_sort_key(chrom: str) -> Tuple[int, int, str]:
    """
    Generate sort key for natural chromosome ordering.
    
    Handles both GRCh37 (1, 2, ..., 22, X, Y) and GRCh38 (chr1, chr2, ..., chr22, chrX, chrY).
    
    Arguments:
        chrom : str
            Chromosome name
        
    Returns:
        tuple
            Sort key: (type_order, numeric_value, string_value)
            - Numeric chromosomes: (0, number, "")
            - X chromosome: (1, 0, "X")
            - Y chromosome: (1, 1, "Y")
            - Other: (2, 0, chrom)
    """
    chrom_str = str(chrom)
    # Remove 'chr' prefix if present
    chrom_clean = chrom_str.replace('chr', '').replace('Chr', '').replace('CHR', '')
    
    # Try to convert to int for numeric chromosomes
    try:
        return (0, int(chrom_clean), "")
    except ValueError:
        # Sex chromosomes
        if chrom_clean.upper() == 'X':
            return (1, 0, "X")
        elif chrom_clean.upper() == 'Y':
            return (1, 1, "Y")
        else:
            # Other chromosomes (shouldn't happen with standard genomes)
            return (2, 0, chrom_clean)

def group_files_by_chromosome(file_list: List[str], chromosomes: List[str]) -> Dict[str, List[str]]:
    """
    Group mutation rate files by chromosome.
    
    Arguments:
    file_list : list of str
        List of all mutation rate file paths
    chromosomes : list of str
        List of expected chromosome names
        
    Returns:
    dict
        Dictionary mapping chromosome -> list of files for that chromosome
    """
    grouped = {chrom: [] for chrom in chromosomes}
    
    for filepath in file_list:
        basename = os.path.basename(filepath)
        # Find which chromosome this file belongs to
        for chrom in chromosomes:
            if f"_{chrom}_mutation_rates" in basename:
                grouped[chrom].append(filepath)
                break
    return grouped


def fetch_existing_chr_index(input_file: str, chrom: str) -> ChromIndex:
    """
    Read existing aggregated stats from a tabix-indexed file for a single chromosome.
    
    The file format has columns:
    chromosome, position, ref_base, A_mean, A_std, A_M2, T_mean, T_std, T_M2, 
    C_mean, C_std, C_M2, G_mean, G_std, G_M2, INS_mean, INS_std, INS_M2, 
    DEL_mean, DEL_std, DEL_M2, n_samples
    
    Parameters:
        input_file : str
            Path to bgzipped + tabix indexed aggregated file
        chrom : str
            Chromosome to fetch
            
    Returns:
        ChromIndex with existing statistics loaded
    """
    chrom_index = ChromIndex(chrom)
    tb = pysam.TabixFile(input_file)

    try:
        for row in tb.fetch(chrom, parser=pysam.asTuple()):
            pos = int(row[1])
            ref_base = row[2]
            n_samples = int(row[-1])
            
            # Create PositionData with the n value
            pos_data = PositionData(ref=ref_base, n=n_samples)
            
            # Parse each mutation type
            for i, mut_type in enumerate(MUTATION_TYPES):
                col_offset = 3 + (i * 3)
                mean = float(row[col_offset])
                std = float(row[col_offset + 1])
                M2 = float(row[col_offset + 2])

                # Store MutationStats (without n)
                pos_data.mutations[mut_type] = MutationStats(mean=mean, M2=M2)
            
            chrom_index.all_stats[pos] = pos_data
    
    except (StopIteration, ValueError) as e:
        # StopIteration: normal end of iteration
        # ValueError: chromosome not in file (e.g., "could not create iterator for region")
        if isinstance(e, ValueError) and "could not create iterator" in str(e):
            print(f"Warning: Chromosome {chrom} not found in existing file. Starting fresh for this chromosome.")
        elif isinstance(e, ValueError):
            # Some other ValueError - re-raise it
            raise
        # If StopIteration, just pass (normal end)
    finally:
        tb.close()

    return chrom_index

def process_chromosome(args: Tuple) -> Optional[pd.DataFrame]:
    """
    Process mutation rates for a single chromosome.
    
    This function aggregates mutation rates across samples for one chromosome,
    and optionally merges with existing aggregated data.
    
    Arguments:
    args : tuple
        (chromosome, new_files, existing_df, mutation_types)
        - chromosome: str, chromosome name
        - new_files: list of str, paths to new mutation rate files for this chromosome
        - update: bool, if true read in chromosome stats from already existing file
        - mutation_types: list of str, mutation type column names
        
    Returns:
    --------
    pd.DataFrame or None
        Aggregated data for this chromosome, or None if no data
    """
    chromosome, new_files, update, existing_file, mutation_types, temp_dir = args
    
    # STEP 1: Load existing index if it exists
    if update and existing_file:
        MasterChrIndex = fetch_existing_chr_index(existing_file, chromosome)
    else:
        MasterChrIndex = ChromIndex(chromosome)
        
    # STEP 2: Process new sample files ONE AT A TIME
    for file in new_files:
        new_sample = pd.read_csv(file, sep="\t", compression="gzip")
        
        for row in new_sample.itertuples(index=False):
            pos = row.position
            ref = row.ref_base
            
            # Collect all mutation values for this position
            mutation_values = {mut_type: getattr(row, mut_type) for mut_type in mutation_types}
            
            # Add to ChromIndex (handles creation if needed)
            MasterChrIndex.add_sample(pos, ref, mutation_values)
    
    # STEP 3: Write out chromosome results as temp file
    chrom_df = MasterChrIndex.get_all_positions_as_df()
    if chrom_df.empty:
        return (chromosome, None)
    
    temp_file = os.path.join(temp_dir, f"{chromosome}_aggregated_temp.tsv")
    chrom_df.to_csv(temp_file, sep="\t", index=False)

    return  (chromosome, temp_file)


def aggregate_all_chromosomes(new_files: List[str], update: bool, existing_file: str,
                              chromosomes: List[str], mutation_types: List[str],
                              threads: int, temp_dir: str,
                              log_file) -> List[str]:
    """
    Aggregate mutation rates across all chromosomes in parallel.

    Processes each chromosome independently, aggregating new sample data and merging
    with existing aggregated data if present. Writes results to temporary files for
    memory-efficient streaming concatenation.

    Arguments:
        new_files : list of str
            New mutation rate files to process ({sample_id}_{chrom}_mutation_rates.tsv.gz)
        existing_file : str
            Path to existing aggregated file (for incremental updates), or empty string
        chromosomes : list of str
            List of chromosome names (e.g., ['chr1', 'chr2', ...] or ['1', '2', ...])
        mutation_types : list of str
            Mutation type column names
        threads : int
            Number of parallel worker processes
        temp_dir : str
            Directory for temporary chromosome files
        
    Returns:
        list of str
            Paths to temporary chromosome files in natural order (headerless TSV format)
    """
    # get list of all new input files by chrom
    grouped_files = group_files_by_chromosome(new_files, chromosomes)

    # Prepare arguments for parallel processing
    process_args = []
    for chrom in chromosomes:
        new_chrom_files = grouped_files.get(chrom, [])
        
        # Skip chromosomes with no new files and not updating
        if len(new_chrom_files) == 0 and not update:
            continue
        # Args: (chromosome, new_files, update, existing_file, mutation_types, temp_dir)
        process_args.append((chrom, new_chrom_files, update, existing_file, mutation_types, temp_dir))

    log(f"Processing {len(process_args)} chromosomes in parallel with {threads} threads...", log_file)
    # Process chromosomes in parallel - returns temp file paths
    temp_pairs = []
    with Pool(processes=threads) as pool:
        for i, result in enumerate(pool.imap_unordered(process_chromosome, process_args), 1):
            chrom, temp_file = result
            if temp_file is not None:
                temp_pairs.append((chrom, temp_file))
            log(f"  Completed: {i}/{len(process_args)} chromosomes", log_file)
    
    # Sort by natural chrom order and return file paths
    temp_pairs.sort(key=lambda x: natural_sort_key(x[0]))

    return [f for _, f in temp_pairs]

def bgzip_and_index(tsv_file: str, output_gz: str, seq_col: int = 0, 
                   start_col: int = 1, end_col: int = 1) -> str:
    """
    Compress TSV file with bgzip and create tabix index.
    
    Arguments:
        tsv_file : str
            Path to uncompressed TSV file
        output_gz : str
            Path for output .tsv.gz file
        seq_col : int
            Column index for chromosome (0-based)
        start_col : int
            Column index for start position (0-based)
        end_col : int
            Column index for end position (0-based)
        
    Returns:
        str
            Path to tabix index file (.tbi)
    """
    # Write to temporary files first
    temp_gz = output_gz + '.tmp'
    temp_tbi = temp_gz + '.tbi'

    # Compress with bgzip
    pysam.tabix_compress(tsv_file, temp_gz, force=True)
    
    # Create tabix index
    index_file = pysam.tabix_index(
        temp_gz,
        preset=None,
        seq_col=seq_col,
        start_col=start_col,
        end_col=end_col,
        line_skip=1,
        zerobased=True,
        force=True
    )
    shutil.move(temp_gz, output_gz)
    shutil.move(temp_tbi, output_gz + '.tbi')
    return index_file


def append_samples_to_tracker(tracker_file: str, sample_ids: List[str]) -> None:
    """
    Append new sample IDs to the sample tracker file.
    
    Creates the file if it doesn't exist. Each sample ID is written on a new line.
    
    Parameters:
    -----------
    tracker_file : str
        Path to sample tracker file
    sample_ids : list of str
        List of sample IDs to append
    """
    if len(sample_ids) == 0:
        return

    # Create parent directory if needed
    os.makedirs(os.path.dirname(tracker_file), exist_ok=True)
    
    # Append samples (one per line)
    with open(tracker_file, 'a') as f:
        for sample_id in sample_ids:
            f.write(f"{sample_id}\n")


def aggregate_mutation_rates(input_files: List[str], output_file: str, 
                             sample_ids: List[str], chromosomes: List[str],
                             sample_tracker: str, threads: int, log_file) -> None:
    """
    Main aggregation pipeline: process samples and write bgzip-compressed, tabix-indexed output.
    """

    log(f"{'='*60}", log_file)
    log(f"MUTATION RATE AGGREGATION", log_file)
    log(f"{'='*60}\n", log_file)
    log(f"Number of new samples: {len(sample_ids)}", log_file)
    log(f"Number of input files: {len(input_files)}", log_file)
    log(f"Output file: {output_file}", log_file)
    log(f"Threads: {threads}\n", log_file)
    
    # Early exit if no files to process
    if len(input_files) == 0:
        log("No input files to process. Exiting.", log_file)
        return

    log(f"Mutation types: {MUTATION_TYPES}\n", log_file, flush=True)

    # Check if we're updating an existing file
    existing_file = output_file if os.path.exists(output_file) else ""
    update_existing = bool(existing_file)
    
    if update_existing:
        log(f"Existing aggregated file found. Will update with new samples.\n", log_file)
    else:
        log(f"No existing aggregated file. Creating new index.\n", log_file)
    
    # Create temp directory
    temp_dir = os.path.join(os.path.dirname(output_file), 'temp_chroms')
    os.makedirs(temp_dir, exist_ok=True)

    # Get list of temp files (in chromosome order)
    temp_files = aggregate_all_chromosomes(
        new_files=input_files,
        update=update_existing,
        existing_file=existing_file,
        chromosomes=chromosomes,
        mutation_types=MUTATION_TYPES,
        threads=threads,
        temp_dir=temp_dir,
        log_file= log_file
    )
    
    # Concatenate temp files with pandas
    log(f"\nConcatenating {len(temp_files)} chromosome files...", log_file)
    
    dfs = []
    for temp_file in temp_files:
        df = pd.read_csv(temp_file, sep="\t")
        dfs.append(df)
    
    combined_df = pd.concat(dfs, ignore_index=True)
    log(f"  Combined {len(combined_df):,} total positions", log_file)
    
    # Write final output and compress
    log(f"\nCompressing and indexing...", log_file)
    final_temp = output_file.replace('.tsv.gz', '_temp.tsv')
    combined_df.to_csv(final_temp, sep="\t", index=False)
    
    index_file = bgzip_and_index(final_temp, output_file, seq_col=0, start_col=1, end_col=1)
    
    # Clean up temporary files
    log(f"\nCleaning up temporary files...", log_file)
    if os.path.exists(final_temp):
        os.remove(final_temp)
    shutil.rmtree(temp_dir)

    # Update sample tracker
    append_samples_to_tracker(sample_tracker, sample_ids)
    
    # Log summary
    log(f"\n{'='*60}", log_file)
    log(f"AGGREGATION SUMMARY", log_file)
    log(f"{'='*60}", log_file)
    log(f"\nTotal positions: {len(combined_df):,}", log_file)
    log(f"Chromosomes processed: {len(temp_files)}", log_file)
    log(f"New samples added: {len(sample_ids)}", log_file)
    log(f"\nOutput: {output_file}", log_file)
    log(f"Index: {index_file}", log_file)
    log(f"\n{'='*60}\n", log_file)
    
    log(f"Completed successfully!", log_file)

def log(msg, log_file, flush=False):
    print(msg, file=log_file)
    # print(msg)
    if flush:
        log_file.flush()

def filter_processed_samples(input_files: List[str], sample_ids: List[str], 
                             sample_tracker: str) -> Tuple[List[str], List[str]]:
    """
    Filter out samples that have already been processed.
    
    Arguments:
        input_files : list of str
            All input mutation rate files
        sample_ids : list of str
            All sample IDs
        sample_tracker : str
            Path to sample tracker file
            
    Returns:
        tuple: (filtered_input_files, filtered_sample_ids)
            Lists with already-processed samples removed
    """
    # Load already-processed samples
    processed = set()
    if os.path.exists(sample_tracker):
        with open(sample_tracker, 'r') as f:
            processed = {line.strip() for line in f if line.strip()}
    
    # Filter sample IDs
    new_sample_ids = [s for s in sample_ids if s not in processed]
    
    # Filter input files (only keep files that contain a new sample ID)
    new_input_files = []
    for filepath in input_files:
        basename = os.path.basename(filepath)
        # Check if this file belongs to any new sample
        for sample_id in new_sample_ids:
            if basename.startswith(f"{sample_id}_"):
                new_input_files.append(filepath)
                break
    
    return new_input_files, new_sample_ids

def main():
    args = parse_args()
    log_file = open(args.log, 'w') if args.log else sys.stdout

    # log date and time start
    start_time = datetime.datetime.now()
    log(f"Starting mutation rate aggregation: {datetime.datetime.now()}\n", log_file)

    # Filter out already-processed samples
    filtered_input_files, filtered_sample_ids = filter_processed_samples(
        args.input, 
        args.sample_ids, 
        args.sample_tracker
    )

    # Early exit if no new samples
    if len(filtered_sample_ids) == 0:
        log("All samples already processed. Nothing to do.", log_file)
        log(f"\nFinished mutation rate aggregation: {datetime.datetime.now()}", log_file)
        if log_file != sys.stdout:
            log_file.close()
        return

    aggregate_mutation_rates(
        input_files=filtered_input_files,
        output_file=args.output,
        sample_ids=filtered_sample_ids,
        chromosomes=args.chromosomes,
        sample_tracker=args.sample_tracker,
        threads=args.threads,
        log_file=log_file
    )
    # log the date and time finishing the script
    log(f"\nFinished mutation rate aggregation: {datetime.datetime.now()}", log_file)
    finish_time = datetime.datetime.now()
    execution_time= finish_time - start_time
    log(f"Total execution time: {execution_time}", log_file)

    if log_file != sys.stdout:
        log_file.close()

if __name__ == "__main__":
    main()
