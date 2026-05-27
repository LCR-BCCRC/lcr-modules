#!/usr/bin/env python3
"""
Filter a samples table to only the rows needed for pending Snakemake jobs.

Reads `snakemake -n` output from stdin, extracts wildcard values for a
specified rule, and writes a filtered copy of the samples table containing
only those samples.  Pipe the filtered table back as the samples table on
re-run to avoid wasting scheduler resources on samples that are already done.

Usage
-----
    snakemake -n [--snakefile Snakefile.smk] 2>&1 | \\
        python utils/filter_samples_dryrun.py \\
            --rule _run_battenberg \\
            --samples /path/to/samples.tsv \\
            --match tumour_id=sample_id \\
            --match normal_id=sample_id \\
            --output samples_pending.tsv

Arguments
---------
--rule / -r     Rule name whose pending jobs define the sample set.
                Use the internal rule name (e.g. _run_battenberg).
--samples / -s  Path to the samples TSV (or CSV with --sep ,).
--match / -m    Map a Snakemake wildcard key to a samples-table column,
                as WILDCARD=COLUMN.  Repeatable.  A row is included when
                its value in COLUMN matches ANY wildcard value across
                ALL --match pairs.  Typical battenberg usage:
                    --match tumour_id=sample_id --match normal_id=sample_id
--output / -o   Output file path.  Omit to write to stdout.
--sep           Field separator for the samples table (default: tab).
--list-ids      Print the matched sample IDs and exit (no table written).

Example: dry-run, then re-run only pending samples
---------------------------------------------------
    # 1. Capture dry-run output and build filtered table
    ./demo/dry-run.sh Snakefile.smk _run_battenberg_all "" runtime_config.yaml \\
        2>&1 | python utils/filter_samples_dryrun.py \\
            --rule _run_battenberg \\
            --samples samples.tsv \\
            --match tumour_id=sample_id \\
            --match normal_id=sample_id \\
            --output samples_battenberg_pending.tsv

    # 2. Inspect
    wc -l samples_battenberg_pending.tsv
    head samples_battenberg_pending.tsv

    # 3. Re-run pointing at the filtered table (in your Snakefile config or
    #    by symlinking / overriding the samples variable)
"""

import argparse
import re
import sys
from collections import defaultdict

try:
    import pandas as pd
except ImportError:
    print("Error: pandas is required.  Install with: pip install pandas", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_wildcards_for_rule(text: str, rule_name: str) -> list[dict]:
    """
    Return a list of wildcard dicts (one per job) for *rule_name* found in
    the `snakemake -n` output *text*.

    Snakemake prints each job block like:

        [timestamp] rule _foo:
            input:  ...
            output: ...
            jobid:  42
            wildcards: key1=val1, key2=val2
            ...

    The timestamp prefix is optional (absent in some Snakemake versions).
    """
    # Locate every occurrence of this rule's header
    rule_header = re.compile(
        r'(?:^\[.*?\]\s+)?rule\s+' + re.escape(rule_name) + r'\s*:',
        re.MULTILINE,
    )
    # A new rule block starts with the same pattern (any rule name)
    any_rule_header = re.compile(
        r'(?:^\[.*?\]\s+)?rule\s+\w+\s*:',
        re.MULTILINE,
    )
    wildcards_line = re.compile(r'^\s+wildcards:\s+(.+)$', re.MULTILINE)

    jobs = []
    for m in rule_header.finditer(text):
        # Delimit this block: from this header to the start of the next rule
        block_start = m.start()
        next_rule = any_rule_header.search(text, m.end())
        block_end = next_rule.start() if next_rule else len(text)
        block = text[block_start:block_end]

        wc_m = wildcards_line.search(block)
        if not wc_m:
            continue
        wc_str = wc_m.group(1).strip()

        # Parse "key=val, key2=val2" — values may contain hyphens, dots, slashes
        wc_dict: dict[str, str] = {}
        for pair in re.split(r',\s+', wc_str):
            if '=' in pair:
                k, v = pair.split('=', 1)
                wc_dict[k.strip()] = v.strip()
        if wc_dict:
            jobs.append(wc_dict)

    return jobs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--rule', '-r', required=True,
                   help='Snakemake rule name (e.g. _run_battenberg)')
    p.add_argument('--samples', '-s', required=True,
                   help='Path to samples table (TSV by default)')
    p.add_argument('--match', '-m', action='append', default=[],
                   metavar='WILDCARD=COLUMN',
                   help='Map wildcard key → samples column (repeatable)')
    p.add_argument('--output', '-o', default=None,
                   help='Output file (default: stdout)')
    p.add_argument('--sep', default='\t',
                   help='Samples table field separator (default: tab)')
    p.add_argument('--list-ids', action='store_true',
                   help='Print matched sample IDs and exit without writing table')
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    # ---- validate --match pairs -------------------------------------------
    match_map: dict[str, str] = {}   # wildcard_key → column_name
    for spec in args.match:
        if '=' not in spec:
            print(f"Error: --match '{spec}' must be WILDCARD=COLUMN", file=sys.stderr)
            sys.exit(1)
        wk, col = spec.split('=', 1)
        match_map[wk.strip()] = col.strip()

    if not match_map:
        print("Error: at least one --match WILDCARD=COLUMN is required.", file=sys.stderr)
        sys.exit(1)

    # ---- read dry-run from stdin ------------------------------------------
    if sys.stdin.isatty():
        print("Waiting for `snakemake -n` output on stdin "
              "(pipe it in or redirect a saved log)...", file=sys.stderr)
    text = sys.stdin.read()

    # ---- parse wildcards --------------------------------------------------
    jobs = parse_wildcards_for_rule(text, args.rule)

    if not jobs:
        print(f"No pending jobs found for rule '{args.rule}'.\n"
              f"Check that the rule name is correct and that the dry-run "
              f"output was captured (2>&1).", file=sys.stderr)
        sys.exit(0)

    print(f"Found {len(jobs)} pending job(s) for rule '{args.rule}'.",
          file=sys.stderr)

    # ---- collect matched sample IDs per column ----------------------------
    col_ids: dict[str, set] = defaultdict(set)
    for job in jobs:
        for wk, col in match_map.items():
            if wk in job:
                col_ids[col].add(job[wk])
            else:
                print(f"Warning: wildcard '{wk}' not present in job {job}",
                      file=sys.stderr)

    if args.list_ids:
        for col, ids in sorted(col_ids.items()):
            for sid in sorted(ids):
                print(f"{col}\t{sid}")
        return

    # ---- load + filter samples table -------------------------------------
    try:
        samples = pd.read_csv(args.samples, sep=args.sep, dtype=str)
    except FileNotFoundError:
        print(f"Error: samples table not found: {args.samples}", file=sys.stderr)
        sys.exit(1)

    missing_cols = [c for c in col_ids if c not in samples.columns]
    if missing_cols:
        print(f"Error: column(s) not found in samples table: {missing_cols}",
              file=sys.stderr)
        print(f"Available columns: {list(samples.columns)}", file=sys.stderr)
        sys.exit(1)

    mask = pd.Series(False, index=samples.index)
    for col, ids in col_ids.items():
        mask |= samples[col].isin(ids)

    filtered = samples[mask].reset_index(drop=True)
    n_in, n_out = len(samples), len(filtered)
    print(f"Samples table: {n_in} rows → {n_out} rows after filtering.",
          file=sys.stderr)

    if n_out == 0:
        print("Warning: filtered table is empty.  "
              "Double-check --match column names against the samples table.",
              file=sys.stderr)

    # ---- write output ----------------------------------------------------
    out_text = filtered.to_csv(sep=args.sep, index=False)
    if args.output:
        with open(args.output, 'w') as fh:
            fh.write(out_text)
        print(f"Written to: {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(out_text)


if __name__ == '__main__':
    main()
