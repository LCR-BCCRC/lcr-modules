"""
Drop intervals where the median read count across all case samples is below a
threshold.

Background: V5 (Agilent SureSelect All Exon V5) has ~2.3% of baits in regions
that don't capture cleanly in our Staudt cell-line cohort. The V6+UTR PoN does
capture there, so a PoN-side filter (--minimum-interval-median-percentile)
won't catch them. Instead we identify intervals where the case cohort itself
under-covers and drop them upstream of the PoN build.

Inputs
------
--interval_list  GATK PreprocessIntervals output (.interval_list with SAM-style
                 @-headers and CONTIG/START/END/STRAND/NAME data rows).
--counts         One or more case sample .counts.hdf5 files from
                 CollectReadCounts (initial pass at the unfiltered intervals).
--min_median     Drop intervals whose per-interval median across cases is below
                 this threshold. Default 5.

Output
------
--output         Filtered .interval_list. Same SAM header, only the data rows
                 whose intervals passed the filter.
"""
import argparse
import sys
import h5py
import numpy as np


def read_counts_hdf5(path):
    """Return (counts: np.ndarray of length n_intervals, intervals: list of (contig, start, end))."""
    with h5py.File(path, "r") as f:
        counts = np.array(f["counts/values"])[0]
        idx_se = np.array(f["intervals/transposed_index_start_end"])
        contig_names = [
            c.decode() if isinstance(c, bytes) else c
            for c in np.array(f["intervals/indexed_contig_names"])
        ]
    contig_arr = idx_se[0].astype(int)
    start_arr = idx_se[1].astype(int)
    end_arr = idx_se[2].astype(int)
    intervals = [
        (contig_names[contig_arr[i]], int(start_arr[i]), int(end_arr[i]))
        for i in range(counts.shape[0])
    ]
    return counts, intervals


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--interval_list", required=True)
    p.add_argument("--counts", required=True, nargs="+")
    p.add_argument("--output", required=True)
    p.add_argument("--min_median", type=float, default=5.0)
    args = p.parse_args()

    if len(args.counts) == 0:
        sys.exit("filter_intervals_by_case_coverage: no --counts files provided")

    print(f"Loading {len(args.counts)} case-count HDF5 files", file=sys.stderr)
    first_counts, intervals = read_counts_hdf5(args.counts[0])
    n_intervals = first_counts.shape[0]
    matrix = np.empty((len(args.counts), n_intervals), dtype=np.float64)
    matrix[0] = first_counts
    for i, path in enumerate(args.counts[1:], start=1):
        counts, ivs = read_counts_hdf5(path)
        if counts.shape[0] != n_intervals:
            sys.exit(
                f"filter_intervals_by_case_coverage: {path} has {counts.shape[0]} "
                f"intervals; expected {n_intervals}"
            )
        if ivs != intervals:
            sys.exit(
                f"filter_intervals_by_case_coverage: {path} interval coordinates "
                f"don't match the first file. Were these counts collected against "
                f"the same interval list?"
            )
        matrix[i] = counts

    medians = np.median(matrix, axis=0)
    keep_mask = medians >= args.min_median
    n_keep = int(keep_mask.sum())
    n_drop = n_intervals - n_keep
    print(
        f"Per-interval median across {len(args.counts)} cases: "
        f"min={medians.min():.0f} median={np.median(medians):.0f} max={medians.max():.0f}",
        file=sys.stderr,
    )
    print(
        f"Dropping {n_drop}/{n_intervals} intervals ({100*n_drop/n_intervals:.2f}%) "
        f"with median<{args.min_median}; keeping {n_keep}",
        file=sys.stderr,
    )

    keep_set = {iv for iv, k in zip(intervals, keep_mask) if k}

    with open(args.interval_list) as fh, open(args.output, "w") as out:
        n_data = 0
        n_kept_rows = 0
        for line in fh:
            if line.startswith("@"):
                out.write(line)
                continue
            n_data += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            key = (parts[0], int(parts[1]), int(parts[2]))
            if key in keep_set:
                out.write(line)
                n_kept_rows += 1

    print(
        f"interval_list rewrite: {n_data} data rows in input, {n_kept_rows} kept "
        f"in output",
        file=sys.stderr,
    )

    if n_kept_rows == 0:
        sys.exit("filter_intervals_by_case_coverage: zero intervals kept; threshold too high?")


if __name__ == "__main__":
    main()
