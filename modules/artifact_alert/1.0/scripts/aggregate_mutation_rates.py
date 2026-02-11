#!/usr/bin/env python3
"""
Aggregate Mutation Rates Accross multiple samples.

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

ALL_cols = ['chromosome', 'position', 'ref_base', 'A_mean', 'A_std', 'A_M2', 
            'T_mean', 'T_std', 'T_M2', 'C_mean', 'C_std', 'C_M2',
            'G_mean', 'G_std', 'G_M2', 'INS_mean', 'INS_std', 'INS_M2',
            'DEL_mean', 'DEL_std', 'DEL_M2', 'n_samples']
MUTATION_COLS = ['A', 'T', 'C', 'G', 'INS', 'DEL']

@dataclass
class PositionStats:
    """
    Calculate running mean and variance for a genomic position using Welford's algorithm.
    This data class holds the info for a single genomic position and mutation column.

    Attributes:
        n : int
            Number of samples aggregated
        mean : float
            Running mean value
        M2 : float
            Sum of squared differences from mean (used to calculate variance)
    """
    
    def __init__(self, n: int = 0, mean: float = 0.0, M2: float = 0.0):
        self.n = n
        self.mean = mean
        self.M2 = M2
    
    def add_sample(self, value: float) -> None:
        """
        Add a new sample value to the running statistics.

        Arguments:
            value : float
                The new sample value to incorporate
        """
        self.n += 1
        delta = value - self.mean
        self.mean += delta / self.n
        delta2 = value - self.mean
        self.M2 += delta * delta2
    
    def get_variance(self) -> float:
        """Return sample variance."""
        if self.n < 2:
            return 0.0
        return self.M2 / (self.n - 1)
    
    def get_std(self) -> float:
        """Return sample standard deviation."""
        return np.sqrt(self.get_variance())
    
    def get_stats(self) -> Tuple[float, float, float, int]:
        """
        Return complete statistics tuple.
        
        Returns:
        tuple : (mean, std, M2, n)
            Mean, standard deviation, sum of squares, and sample count
        """
        return self.mean, self.get_std(), self.M2, self.n
    
    @staticmethod
    def combine(stats1: 'PositionStats', stats2: 'PositionStats') -> 'PositionStats':
        """
        Combine two independent PositionStats objects.
        
        Uses the parallel algorithm for combining variances from two groups:
        - Combined n = n1 + n2
        - Combined mean = (n1*mean1 + n2*mean2) / (n1 + n2)
        - Combined M2 = M2_1 + M2_2 + (mean1 - mean2)^2 * n1*n2/(n1+n2)
        
        Parameters:
        -----------
        stats1 : PositionStats
            First set of statistics
        stats2 : PositionStats
            Second set of statistics
            
        Returns:
        --------
        PositionStats
            Combined statistics
        """
        if stats1.n == 0:
            return stats2
        if stats2.n == 0:
            return stats1
        
        n_combined = stats1.n + stats2.n
        mean_combined = (stats1.n * stats1.mean + stats2.n * stats2.mean) / n_combined
        
        # Parallel variance combination formula
        delta = stats2.mean - stats1.mean
        M2_combined = stats1.M2 + stats2.M2 + delta * delta * stats1.n * stats2.n / n_combined
        
        return PositionStats(n=n_combined, mean=mean_combined, M2=M2_combined)

@dataclass
class ChromIndex:
    """
    Holds statistics for all positions and mutation types in a chromosome.
    
    Attributes:
        chrom : str
            Chromosome name
        all_stats : dict
            Maps position (int) -> dict(mut_type (str) -> PositionStats)
            Example: all_stats[12345]['A'] = PositionStats(n=5, mean=0.001, M2=0.00002)
    """
    
    def __init__(self, chrom: str):
        self.chrom = chrom
        self.all_stats: Dict[int, Dict[str, PositionStats]] = {}

    def add_or_update(self, pos: int, ref_base: str, mut_type: str, stats: PositionStats) -> None:
        """
        Add new stats or combine with existing stats for a position and mutation type.
        
        Parameters:
            pos : int
                Genomic position
            mut_type : str
                Mutation type (e.g., 'A', 'T', 'C', 'G', 'INS', 'DEL')
            stats : PositionStats
                Statistics to add or combine
        """
        # if position not in index, instantiate
        if pos not in self.all_stats:
            self.all_stats[pos] = {"ref": ref_base}
        
        if mut_type not in self.all_stats[pos]:
            # New position+mut_type: store directly
            self.all_stats[pos][mut_type] = stats
        else:
            # Already exists: combine
            self.all_stats[pos][mut_type] = PositionStats.combine(
                self.all_stats[pos][mut_type], 
                stats
            )
    
    def get_stats(self, pos: int, mut_type: str) -> Optional[PositionStats]:
        """Get PositionStats for a specific position and mutation type."""
        return self.all_stats.get(pos, {}).get(mut_type)
    
    def get_all_positions_as_df(self) -> pd.DataFrame:
        """Return DataFrame with all positions and mutation stats for this chromosome."""
        sorted_positions = sorted(self.all_stats.keys())
        
        rows = []
        for pos in sorted_positions:
            row = {
                'chromosome': self.chrom,
                'position': pos,
                'ref_base': self.all_stats[pos].get('ref', 'N')
            }
            
            # Get stats for each mutation type
            n_samples = 0
            for mut in MUTATION_COLS:
                stats = self.all_stats[pos].get(mut)
                if stats:
                    row[f'{mut}_mean'] = stats.mean
                    row[f'{mut}_std'] = stats.get_std()
                    row[f'{mut}_M2'] = stats.M2
                    if n_samples == 0:
                        n_samples = stats.n
            
            row['n_samples'] = n_samples
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
    
    # Open tabix file
    tb = pysam.TabixFile(input_file)
    
    mutation_types = ['A', 'T', 'C', 'G', 'INS', 'DEL']
    
    for row in tb.fetch(chrom):
        # Columns: chromosome(0), position(1), ref_base(2), then triplets of mean/std/M2 for each mut, n_samples(last)
        pos = int(row[1])
        ref_base = row[2]
        n_samples = int(row[-1])  # last column
        
        # Parse each mutation type (columns 3-20, in groups of 3)
        for i, mut_type in enumerate(mutation_types):
            col_offset = 3 + (i * 3)
            mean = float(row[col_offset])
            std = float(row[col_offset + 1])
            M2 = float(row[col_offset + 2])
            
            # Reconstruct PositionStats from saved values
            stats = PositionStats(n=n_samples, mean=mean, M2=M2)
            chrom_index.add_or_update(pos, ref_base, mut_type, stats)
    
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
    
    # STEP 1 read in existing index, if exists
    if update and existing_file:
        MasterChrIndex = fetch_existing_chr_index(existing_file, chromosome)
    else:
        MasterChrIndex = ChromIndex(chromosome)
    
    print(f"Processing chromosome {chromosome}: {len(new_files)} new samples")
    

    # STEP 2: Process new sample files ONE AT A TIME (memory efficient)
    for file in new_files:
        new_sample = pd.read_csv(file, sep="\t", compression="gzip")
        
        for row in new_sample.itertuples(index=False):
            pos = row.position
            ref = row.ref_base
            
            for mut_type in mutation_types:
                value = getattr(row, mut_type)
                
                # If stats exist, add sample directly; otherwise create new
                existing_stats = MasterChrIndex.get_stats(pos, mut_type)
                if existing_stats:
                    existing_stats.add_sample(value)
                else:
                    MasterChrIndex.add_or_update(pos, ref, mut_type, PositionStats(n=1, mean=value, M2=0.0))


    # STEP 3: Write out chromosome results as temp file
    chrom_df = MasterChrIndex.get_all_positions_as_df()
    if chrom_df.empty:
        return None
    else:
        temp_file = os.path.join(temp_dir, f"{chromosome}_aggregated_temp.tsv")
        chrom_df.to_csv(temp_file, sep="\t", index=False)
        print(f"  Wrote chromosome {chromosome} stats to {temp_file}")
        return temp_file


def aggregate_all_chromosomes(new_files: List[str], update: bool, existing_file: str,
                              chromosomes: List[str], mutation_types: List[str],
                              threads: int, temp_dir: str) -> List[str]:
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
            Mutation type column names (e.g., ['A', 'T', 'C', 'G', 'INS', 'DEL'])
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
    
    print(f"\nProcessing {len(process_args)} chromosomes in parallel with {threads} threads...")
    
    # Process chromosomes in parallel - returns temp file paths
    with Pool(processes=threads) as pool:
        temp_files = pool.map(process_chromosome, process_args)
    
    # Filter out None results and get chromosomes
    chrom_file_pairs = []
    for i, temp_file in enumerate(temp_files):
        if temp_file is not None:
            chrom = process_args[i][0]  # chromosome is first element of args tuple
            chrom_file_pairs.append((chrom, temp_file))
    
    if len(chrom_file_pairs) == 0:
        raise ValueError("No valid data produced from aggregation")
    
    # Sort by chromosome order
    chrom_file_pairs.sort(key=lambda x: natural_sort_key(x[0]))
    
    print(f"\nProcessed {len(chrom_file_pairs)} chromosomes")
    return [f for _, f in chrom_file_pairs]  # Return file paths in order

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
    # Compress with bgzip
    pysam.tabix_compress(tsv_file, output_gz, force=True)
    
    # Create tabix index
    print(f"Creating tabix index...")
    index_file = pysam.tabix_index(
        output_gz,
        preset=None,
        seq_col=seq_col,
        start_col=start_col,
        end_col=end_col,
        meta_char='#',
        zerobased=True,
        force=True
    )
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
    
    print(f"\nUpdating sample tracker: {tracker_file}")
    
    # Create parent directory if needed
    os.makedirs(os.path.dirname(tracker_file), exist_ok=True)
    
    # Append samples (one per line)
    with open(tracker_file, 'a') as f:
        for sample_id in sample_ids:
            f.write(f"{sample_id}\n")
    
    print(f"  Added {len(sample_ids)} samples to tracker")


def aggregate_mutation_rates(input_files: List[str], output_file: str, 
                             sample_ids: List[str], chromosomes: List[str],
                             sample_tracker: str, threads: int, log_file) -> None:
    """
    Main aggregation pipeline: process samples and write bgzip-compressed, tabix-indexed output.
    """
    def log(msg, flush=False):
        print(msg, file=log_file)
        print(msg)
        if flush:
            log_file.flush()
    
    log(f"{'='*60}")
    log(f"MUTATION RATE AGGREGATION")
    log(f"{'='*60}\n")
    log(f"Number of new samples: {len(sample_ids)}")
    log(f"Number of input files: {len(input_files)}")
    log(f"Output file: {output_file}")
    log(f"Threads: {threads}\n")
    
    # Early exit if no files to process
    if len(input_files) == 0:
        log("No input files to process. Exiting.")
        return

    log(f"Mutation types: {MUTATION_COLS}\n", flush=True)

    # Check if we're updating an existing file
    existing_file = output_file if os.path.exists(output_file) else ""
    update_existing = bool(existing_file)
    
    if update_existing:
        log(f"Existing aggregated file found. Will update with new samples.\n")
    else:
        log(f"No existing aggregated file. Creating new index.\n")
    
    # Create temp directory
    temp_dir = os.path.join(os.path.dirname(output_file), 'temp_chroms')
    os.makedirs(temp_dir, exist_ok=True)

    # Get list of temp files (in chromosome order)
    temp_files = aggregate_all_chromosomes(
        new_files=input_files,
        update=update_existing,
        existing_file=existing_file,
        chromosomes=chromosomes,
        mutation_types=MUTATION_COLS,
        threads=threads,
        temp_dir=temp_dir
    )
    
    # Concatenate temp files with pandas
    log(f"\nConcatenating {len(temp_files)} chromosome files...")
    
    dfs = []
    for temp_file in temp_files:
        df = pd.read_csv(temp_file, sep="\t")
        dfs.append(df)
    
    combined_df = pd.concat(dfs, ignore_index=True)
    log(f"  Combined {len(combined_df):,} total positions")
    
    # Write final output and compress
    log(f"\nCompressing and indexing...")
    final_temp = output_file.replace('.tsv.gz', '_temp.tsv')
    combined_df.to_csv(final_temp, sep="\t", index=False)
    
    index_file = bgzip_and_index(final_temp, output_file, seq_col=0, start_col=1, end_col=1)
    
    # Clean up temporary files
    log(f"\nCleaning up temporary files...")
    if os.path.exists(final_temp):
        os.remove(final_temp)
    shutil.rmtree(temp_dir)

    # Update sample tracker
    append_samples_to_tracker(sample_tracker, sample_ids)
    
    # Log summary
    log(f"\n{'='*60}")
    log(f"AGGREGATION SUMMARY")
    log(f"{'='*60}")
    log(f"\nTotal positions: {len(combined_df):,}")
    log(f"Chromosomes processed: {len(temp_files)}")
    log(f"New samples added: {len(sample_ids)}")
    log(f"\nOutput: {output_file}")
    log(f"Index: {index_file}")
    log(f"\n{'='*60}\n")
    
    log(f"Completed successfully!")


def main():
    """
    Main entry point for Snakemake workflow.
    
    Expected Snakemake parameters:
    - input.mutation_files: list of mutation rate files
    - output.aggregated: output .tsv.gz file path
    - output.index: output .tbi index file path
    - params.sample_ids: list of sample IDs being processed
    - params.chromosomes: list of chromosome names
    - params.sample_tracker: path to sample tracker file
    - threads: number of threads for parallel processing
    - log: log file path
    """
    # Open log file
    log_file = open(snakemake.log[0], 'w')
    
    try:
        aggregate_mutation_rates(
            input_files=snakemake.input.mutation_files,
            output_file=snakemake.output.aggregated,
            sample_ids=snakemake.params.sample_ids,
            chromosomes=snakemake.params.chromosomes,
            sample_tracker=snakemake.params.sample_tracker,
            threads=snakemake.threads,
            mutation_types=MUTATION_COLS,
            log_file=log_file
        )
    except Exception as e:
        print(f"ERROR: {str(e)}", file=log_file)
        print(f"ERROR: {str(e)}")
        raise
    finally:
        log_file.close()


if __name__ == "__main__":
    main()
