#!/usr/bin/env python3

import pandas as pd
import re
from dataclasses import dataclass
from typing import Optional, Dict
import gzip

@dataclass
class BaseCounts:
    """Base counts at a genomic position"""
    A: int = 0
    T: int = 0  
    G: int = 0
    C: int = 0
    ref: int = 0
    deletions: int = 0
    insertions: int = 0
    
    def get_total(self) -> int:
        """Return total number of reads (NOT including indels as separate reads)"""
        return self.A + self.T + self.G + self.C + self.ref
    
    def get_base_count(self, base: str) -> int:
        """Get count for a specific base"""
        base_upper = base.upper()
        if base_upper in ['A', 'T', 'G', 'C']:
            return getattr(self, base_upper)
        return 0


class PileupParser:
    """Parses samtools mpileup format"""
    
    def __init__(self):
        self.read_start_pattern = re.compile(r'\^.')
    
    def count_indels(self, bases: str) -> Dict[str, int]:
        """
        Count insertions and deletions in bases string.
        
        - Insertions: Count +N[seq] patterns (insertions occur AFTER this position)
        - Deletions: Count * and # characters (reads with deletion AT this position)
        """
        if pd.isna(bases) or not isinstance(bases, str):
            return {'insertions': 0, 'deletions': 0}
        
        insertions = 0
        deletions = 0
        i = 0
        
        while i < len(bases):
            char = bases[i]
            
            if char == '+':
                # Insertion after this position - count it
                insertions += 1
                # Skip past: parse length, then skip that many bases
                i += 1
                length_match = re.match(r'(\d+)', bases[i:])
                if length_match:
                    indel_len = int(length_match.group(1))
                    i += len(length_match.group(1)) + indel_len
                continue
            
            elif char == '-':
                # Deletion lookahead - do NOT count here, just skip past it
                # The actual deletions are marked with * at subsequent positions
                i += 1
                length_match = re.match(r'(\d+)', bases[i:])
                if length_match:
                    indel_len = int(length_match.group(1))
                    i += len(length_match.group(1)) + indel_len
                continue
            
            elif char in '*#':
                # Deletion AT this position (CIGAR "D")
                deletions += 1
            
            elif char == '^':
                # Read start marker - skip it and the following quality char
                i += 2
                continue
            
            i += 1
        
        return {'insertions': insertions, 'deletions': deletions}
    
    def clean_bases_string(self, bases: str) -> str:
        """Remove indel sequences, deletions, and read markers from bases string."""
        if pd.isna(bases) or not isinstance(bases, str):
            return ""
        
        cleaned = ""
        i = 0
        
        while i < len(bases):
            char = bases[i]
            
            if char in '+-':
                # Skip indel notation: +/-N[seq]
                i += 1
                length_match = re.match(r'(\d+)', bases[i:])
                if length_match:
                    indel_len = int(length_match.group(1))
                    i += len(length_match.group(1)) + indel_len
                continue
            
            elif char == '^':
                # Skip read start marker and quality char
                i += 2
                continue
            
            elif char == '$':
                # Skip read end marker
                i += 1
                continue
            
            elif char in '*#':
                # Skip deletion markers (not a base call)
                i += 1
                continue
            
            else:
                cleaned += char
                i += 1
        
        return cleaned

    def parse_position(self, bases: str, ref_base: str) -> BaseCounts:
        """Parse pileup bases string into base counts"""
        # Count indels
        indel_counts = self.count_indels(bases)
        
        # Clean the bases string (removes indel notation, read markers, deletion markers)
        cleaned_bases = self.clean_bases_string(bases)
        
        # Initialize counts with indels
        counts = BaseCounts(
            deletions=indel_counts['deletions'],
            insertions=indel_counts['insertions']
        )
        
        # Count each base type
        for base in cleaned_bases:
            if base in '.,':  # Reference matches
                counts.ref += 1
            elif base.upper() == 'A':
                counts.A += 1
            elif base.upper() == 'T':
                counts.T += 1
            elif base.upper() == 'G':
                counts.G += 1
            elif base.upper() == 'C':
                counts.C += 1
            # Skip unrecognized characters (like > < for reference skips)
        
        return counts


class MutationRateCalculator:
    """Calculates mutation rates by type from pileup data"""
    
    def __init__(self, min_depth: int):
        self.min_depth = min_depth
        self.parser = PileupParser()
    
    def passes_filters(self, depth: int, total_bases: int) -> bool:
        """Check if position passes quality filters"""
        return (depth >= self.min_depth and total_bases >= self.min_depth)
    
    def analyze_position(self, chromosome: str, position: int, ref_base: str, 
                        depth: int, bases: str) -> Optional[Dict]:
        """Analyze a single genomic position for mutations by type"""
        
        # Parse base counts
        base_counts = self.parser.parse_position(bases, ref_base)
        total_bases = base_counts.get_total()
        
        # Apply filters
        if not self.passes_filters(depth, total_bases):
            return None
        
        ref_upper = ref_base.upper()
        
        # Calculate rates for all bases (ref will be 0)
        result = {
            'chromosome': chromosome,
            'position': position,
            'ref_base': ref_upper,
            'A': base_counts.A / total_bases if ref_upper != 'A' else 0.0,
            'T': base_counts.T / total_bases if ref_upper != 'T' else 0.0,
            'C': base_counts.C / total_bases if ref_upper != 'C' else 0.0,
            'G': base_counts.G / total_bases if ref_upper != 'G' else 0.0,
            'INS': base_counts.insertions / total_bases,
            'DEL': base_counts.deletions / total_bases
        }
        
        return result


def process_pileup_file(input_file: str, output_file: str, min_depth: int):
    """Process entire pileup file and write mutation rate results"""
    
    print(f"Reading pileup file: {input_file}")
    
    # Define column names
    columns = ['chromosome', 'position', 'ref_base', 'depth', 'bases', 'qualities']
    
    # Read pileup data
    df = pd.read_csv(input_file, sep='\t', names=columns, header=None)
    
    print(f"Processing {len(df)} positions...")
    
    # Initialize calculator
    calculator = MutationRateCalculator(min_depth)
    
    # Process each position
    results = []
    for _, row in df.iterrows():
        result = calculator.analyze_position(
            row['chromosome'], row['position'], row['ref_base'],
            row['depth'], row['bases']
        )
        if result:
            results.append(result)
    
    # Write results
    write_results(results, output_file)


def write_results(results: list, output_file: str):
    """Write mutation analysis results to file"""
    
    if results:
        # Create and sort results dataframe
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values(['chromosome', 'position'])
        
        # Calculate summary statistics
        rate_cols = ['A', 'T', 'C', 'G', 'INS', 'DEL']
        
        print(f"Writing {len(results_df)} positions to: {output_file}")
        print(f"\nMean rates by type:")
        for col in rate_cols:
            print(f"  {col}: {results_df[col].mean():.8f}")
        
        # Write to file
        results_df.to_csv(output_file, sep='\t', index=False, compression="gzip")
        
    else:
        print("No positions met criteria - creating empty output file")
        empty_df = pd.DataFrame(columns=['chromosome', 'position', 'ref_base', 
                                         'A', 'T', 'C', 'G', 'INS', 'DEL'])
        empty_df.to_csv(output_file, sep='\t', index=False, compression="gzip")


def main():
    """Main entry point for Snakemake workflow"""
    # Extract parameters
    min_depth = snakemake.params.min_depth
    
    # Process the file
    process_pileup_file(
        snakemake.input.pileup,
        snakemake.output.mutation_rates,
        min_depth
    )

if __name__ == "__main__":
    main()
