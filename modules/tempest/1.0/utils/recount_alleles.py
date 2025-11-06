#!/usr/bin/env python
"""
Recompute tumor allele counts for variants in an input MAF by piling up a tumor BAM/CRAM with pysam.
Writes a MAF with the same rows and coordinates, updating only t_depth, t_ref_count, and t_alt_count.

CLI:
  -m/--maf <MAF>    -b/--bam <BAM|CRAM>    -o/--output <MAF>
  -g/--genome <build>    -t/--threads <int>

Behavior:
- Parallelizes by chromosome inferred from the genome build (chr1–22,X,Y or 1–22,X,Y),
  restricted to chromosomes present in the MAF (skips empty ones).
- For each variant, runs a pysam pileup in a ±100 bp window and computes counts:
  - SNP: collapses paired reads by query_name (1 fragment counted once).
  - DNP: fully matches all bases of the multi-base substitution (case-insensitive).
  - INS: counts '+' tokens at the variant start; may double-count overlapping pairs (kept as-is).
  - DEL: counts '*' within [Start_Position, End_Position]; alt_count = max across the span.
- Coordinates and all non-count columns are preserved unchanged.

Assumptions/requirements:
- MAF has columns: Chromosome, Start_Position, End_Position, Reference_Allele,
  Tumor_Seq_Allele2, Variant_Type, t_ref_count, t_alt_count, t_depth.
- Variant_Type ∈ {SNP, DNP, INS, DEL}. CRAM reading may require reference availability/config.
- Count fields are written as strings to avoid float coercion.

Notes:
- No normal BAM support and no merging of “missing” variants; this is a recount-only step.
"""
import argparse
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import pysam
import multiprocessing

verbose: bool = False
padding: int = 100


def get_genome_chromosomes(genome_build: str) -> List[str]:
    """
    Return a list of chromosome names for the specified genome build.

    Parameters:
    - genome_build: Genome identifier (e.g., 'hg38', 'grch38', 'hg19-reddy').

    Returns:
    - List of chromosome names in the form:
      ['chr1'..'chr22','chrX','chrY'] for hg38/grch38-like builds,
      otherwise ['1'..'22','X','Y'].
    """
    if genome_build.lower() in ["hg38", "hg19-reddy", "grch38"]:
        return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    else:
        return [str(i) for i in range(1, 23)] + ["X", "Y"]


def get_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
    - argparse.Namespace with:
      maf (str): Input MAF path.
      bam (str): Input BAM/CRAM path.
      output (str): Output MAF path.
      genome (str): Genome build identifier.
      threads (int): Number of worker processes to use.
    """
    parser = argparse.ArgumentParser(description="Annotates allele counts")
    parser.add_argument("-m", "--maf", required=True, help="Input MAF file")
    parser.add_argument("-b", "--bam", metavar="BAM", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", metavar="MAF", required=True, help="Output MAF file")
    parser.add_argument("-g", "--genome", metavar="GENOME", required=True, help="Genome build")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use")
    return parser.parse_args()


def get_read_support_with_samfile(
    samfile: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    ref_base: str,
    alt_base: str,
    mut_class: str,
    min_base_qual: int = 10,
    min_mapping_qual: int = 1,
) -> Dict[str, str]:
    """
    Compute tumor allele counts for a single variant using an already-opened AlignmentFile.

    Parameters:
    - samfile: Open pysam AlignmentFile (BAM/CRAM) to read from.
    - chrom: Chromosome name of the variant.
    - start: 1-based Start_Position from the MAF (inclusive).
    - end: 1-based End_Position from the MAF (inclusive).
    - ref_base: Reference allele (string, may be multi-base for DNP/DEL).
    - alt_base: Tumor_Seq_Allele2 allele (string, may be multi-base for DNP/INS).
    - mut_class: Variant type ('SNP', 'DNP', 'INS', 'DEL').
    - min_base_qual: Minimum base quality to include in pileup.
    - min_mapping_qual: Minimum mapping quality to include in pileup.

    Returns:
    - Dict[str, str] with keys:
      'tumor_depth', 'tumor_ref_count', 'tumor_alt_count' (string-typed).
    """
    # Normalize case for comparisons where sequences are used
    ref_base_cmp = ref_base.upper() if isinstance(ref_base, str) else ref_base
    alt_base_cmp = alt_base.upper() if isinstance(alt_base, str) else alt_base

    actual_start = start - 1  # 0-based for pysam pileup comparison
    del_positions: List[int] = []

    if mut_class == "DEL":
        del_positions = list(range(actual_start, end))

    padded_start = int(start) - padding
    padded_end = int(end) + padding
    chrom = str(chrom)

    ref_count = 0
    alt_count = 0
    total_depth = 0

    read_bases: Dict[str, List[str]] = {}
    del_support: Dict[int, int] = {}
    max_depth = 0

    for pos in del_positions:
        del_support[pos] = 0

    for pileupcolumn in samfile.pileup(
        chrom,
        padded_start,
        padded_end,
        ignore_overlaps=True,  # keep as original behavior
        min_base_quality=min_base_qual,
        min_mapping_quality=min_mapping_qual,
        truncate=True,
    ):
        if mut_class == "INS":
            # Count '+' tokens at the variant start; may double-count overlapping reads
            if pileupcolumn.pos == actual_start:
                seqs = pileupcolumn.get_query_sequences(add_indels=True)
                for base in seqs:
                    if "+" in base:
                        alt_count += 1
                    else:
                        ref_count += 1
                    total_depth += 1

        elif mut_class == "DEL":
            # Track deletion support across the deletion span; use max across sites
            if pileupcolumn.pos in del_positions:
                seqs = pileupcolumn.get_query_sequences(add_indels=True)
                ref_count = 0
                if len(seqs) > max_depth:
                    max_depth = len(seqs)
                for base in seqs:
                    if base == "*":
                        del_support[pileupcolumn.pos] += 1
                    else:
                        ref_count += 1

        elif pileupcolumn.pos == actual_start:
            # SNP/DNP: collect per-fragment observations to collapse read pairs
            for pileupread in pileupcolumn.pileups:
                if mut_class == "SNP" or mut_class == "DNP":
                    if pileupread.is_del or pileupread.is_refskip:
                        continue
                    aln = pileupread.alignment
                    qpos = pileupread.query_position
                    if qpos is None:
                        continue
                    qseq = aln.query_sequence

                    if mut_class == "SNP":
                        if not (0 <= qpos < len(qseq)):
                            continue
                        this_seq = qseq[qpos].upper()
                    else:  # DNP: match the full length of the change
                        var_len = len(ref_base_cmp) if isinstance(ref_base_cmp, str) else 1
                        if qpos + var_len > len(qseq):
                            continue
                        this_seq = qseq[qpos : qpos + var_len].upper()

                    thisname = aln.query_name
                    if thisname in read_bases:
                        read_bases[thisname].append(this_seq)
                    else:
                        read_bases[thisname] = [this_seq]

    # Summarize counts post pileup
    if mut_class == "DEL":
        for position in del_support:
            if del_support[position] > alt_count:
                alt_count = del_support[position]
        ref_count = max_depth - alt_count
        total_depth = max_depth

    if mut_class == "SNP" or mut_class == "DNP":
        for readname in read_bases:
            bases = read_bases[readname]
            # At most 2 observations from the same fragment/pair
            if len(bases) == 1:
                total_depth += 1
                if bases[0] == ref_base_cmp:
                    ref_count += 1
                elif bases[0] == alt_base_cmp:
                    alt_count += 1
            elif bases[0] == bases[1]:
                total_depth += 1
                if bases[0] == ref_base_cmp:
                    ref_count += 1
                elif bases[0] == alt_base_cmp:
                    alt_count += 1
            else:
                total_depth += 1
                if bases[0] == ref_base_cmp or bases[1] == ref_base_cmp:
                    ref_count += 1
                elif bases[0] == alt_base_cmp and bases[1] == alt_base_cmp:
                    alt_count += 1

    # Return as strings to prevent pandas from converting to floats
    return {
        "tumor_depth": str(total_depth),
        "tumor_ref_count": str(ref_count),
        "tumor_alt_count": str(alt_count),
    }


def subset_and_run(
    full_index_maf: pd.DataFrame,
    bamfile: str,
    this_chromosome: str,
    _unused_maf_path: Optional[str],
) -> pd.DataFrame:
    """
    Recount allele support for all variants on a single chromosome.

    Parameters:
    - full_index_maf: Full input MAF as a DataFrame.
    - bamfile: Path to the input BAM/CRAM file (tumor).
    - this_chromosome: Chromosome name to process.
    - _unused_maf_path: Unused positional argument retained for interface compatibility.

    Returns:
    - DataFrame with the same rows as the chromosome subset of the input MAF,
      but with t_depth, t_ref_count, and t_alt_count recomputed and overwritten.
    """
    # Subset to rows for the requested chromosome
    to_process_maf = full_index_maf[full_index_maf.Chromosome == this_chromosome].copy()

    total_to_process = to_process_maf.shape[0]
    if verbose:
        print(f"will work on {total_to_process} MAF rows using {bamfile}")

    # Open BAM/CRAM once per worker
    realname = str(Path(bamfile).resolve())
    samfile = pysam.AlignmentFile(realname, "rc" if realname.endswith("cram") else "rb")

    try:
        ref_counts: List[str] = []
        alt_counts: List[str] = []
        depths: List[str] = []

        # Iterate over variants
        for idx, row in enumerate(to_process_maf.itertuples(index=False)):
            read_counts = get_read_support_with_samfile(
                samfile=samfile,
                chrom=str(row.Chromosome),
                start=int(row.Start_Position),
                end=int(row.End_Position),
                ref_base=str(row.Reference_Allele),
                alt_base=str(row.Tumor_Seq_Allele2),
                mut_class=str(row.Variant_Type),
            )
            ref_counts.append(read_counts["tumor_ref_count"])
            alt_counts.append(read_counts["tumor_alt_count"])
            depths.append(read_counts["tumor_depth"])

            if not idx % 500 and verbose and idx > 0:
                percent_complete = int(100 * idx / total_to_process) if total_to_process else 100
                print(f"{idx} of {total_to_process} rows ({percent_complete}%) complete (chromosome {this_chromosome})")
    finally:
        samfile.close()

    # Overwrite count fields
    to_process_maf.t_depth = depths
    to_process_maf.t_ref_count = ref_counts
    to_process_maf.t_alt_count = alt_counts

    return to_process_maf


def main(args: Optional[argparse.Namespace] = None) -> None:
    """
    Program entry point: parse args, load MAF, dispatch per-chromosome workers, and write output.

    Parameters:
    - args: Optional argparse.Namespace; if None, CLI arguments are parsed via get_args().

    Side effects:
    - Writes the output MAF to the path specified by --output.
    """
    if args is None:
        args = get_args()

    # Load index MAF
    index_maf = pd.read_csv(
        args.maf,
        sep="\t",
        comment="#",
        dtype={"Chromosome": "string", "t_ref_count": "string", "t_alt_count": "string", "t_depth": "string"},
        low_memory=False,
    )

    # Build chromosome list from genome build, but restrict to those present in the MAF
    build_chroms = get_genome_chromosomes(args.genome)
    maf_chroms = (
        index_maf["Chromosome"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    chromosomes = [c for c in build_chroms if c in maf_chroms]

    if not chromosomes and verbose:
        print("No matching chromosomes between genome build and MAF; nothing to do.")

    pool_args = [[index_maf, args.bam, chrom, args.maf] for chrom in chromosomes]

    cool_pool = multiprocessing.Pool(processes=args.threads)
    results = cool_pool.starmap(subset_and_run, pool_args)
    cool_pool.close()
    cool_pool.join()  # Keeping Chris' variable naming

    full_augmented_maf_merged = pd.concat(results) if results else index_maf
    full_augmented_maf_merged.to_csv(args.output, sep="\t", index=False)


main()
