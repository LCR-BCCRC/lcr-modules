#!/usr/bin/env python3
"""
List sample IDs with pending Snakemake jobs, or filter a samples table to them.

Reads `snakemake -n` output from stdin and extracts wildcard values for the
relevant rules.  By default prints the pending sample IDs one per line.
Pass --samples to instead write a filtered copy of the samples table.

Default usage (print pending IDs)
----------------------------------
    ./demo/dry-run.sh Snakefile.smk all "" runtime_config.yaml 2>&1 | \\
        utils/filter_samples_dryrun.py

Filter to a specific rule
--------------------------
    ... | utils/filter_samples_dryrun.py --rule _run_battenberg

Write a filtered samples table (requires pandas)
-------------------------------------------------
    ... | utils/filter_samples_dryrun.py \\
            --samples /path/to/samples.tsv \\
            --output samples_pending.tsv

Arguments
---------
--rule / -r     Rule name to collect wildcards from.  Default: "all" —
                collect wildcards from every rule block in the output.
--match / -m    Map a wildcard key to a samples-table column as
                WILDCARD=COLUMN (repeatable).
                Default: tumour_id=sample_id and normal_id=sample_id
--samples / -s  Path to samples table.  When provided, writes a filtered
                copy instead of printing IDs.  Requires pandas.
--output / -o   Output file for the filtered table (default: stdout).
                Ignored unless --samples is given.
--sep           Field separator for the samples table (default: tab).
"""

import argparse
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_wildcard_str(wc_str: str) -> dict:
    """Parse a single 'key=val, key2=val2' wildcards string into a dict."""
    wc_dict: dict[str, str] = {}
    for pair in re.split(r',\s+', wc_str.strip()):
        if '=' in pair:
            k, v = pair.split('=', 1)
            wc_dict[k.strip()] = v.strip()
    return wc_dict


def parse_wildcards_for_rule(text: str, rule_name: str) -> list[dict]:
    """
    Return a list of wildcard dicts (one per job) for *rule_name*.
    When *rule_name* is 'all', collect from every rule block.
    """
    any_rule_header = re.compile(
        r'(?:^\[.*?\]\s+)?rule\s+\w+\s*:',
        re.MULTILINE,
    )
    target_header = any_rule_header if rule_name == 'all' else re.compile(
        r'(?:^\[.*?\]\s+)?rule\s+' + re.escape(rule_name) + r'\s*:',
        re.MULTILINE,
    )
    wildcards_line = re.compile(r'^\s+wildcards:\s+(.+)$', re.MULTILINE)

    jobs = []
    for m in target_header.finditer(text):
        block_start = m.start()
        next_rule = any_rule_header.search(text, m.end())
        block_end = next_rule.start() if next_rule else len(text)
        block = text[block_start:block_end]

        wc_m = wildcards_line.search(block)
        if not wc_m:
            continue
        wc_dict = _parse_wildcard_str(wc_m.group(1))
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
    p.add_argument('--rule', '-r', default='all',
                   help='Rule name to collect wildcards from (default: all rules)')
    p.add_argument('--match', '-m', action='append', default=None,
                   metavar='WILDCARD=COLUMN',
                   help='Wildcard→column mapping (default: tumour_id=sample_id '
                        'and normal_id=sample_id)')
    p.add_argument('--samples', '-s', default=None,
                   help='Samples table path.  When given, writes a filtered table '
                        'instead of printing IDs.  Requires pandas.')
    p.add_argument('--output', '-o', default=None,
                   help='Output path for filtered table (default: stdout)')
    p.add_argument('--sep', default='\t',
                   help='Samples table separator (default: tab)')
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    # ---- match defaults ---------------------------------------------------
    if args.match is None:
        args.match = ['tumour_id=sample_id', 'normal_id=sample_id']

    match_map: dict[str, str] = {}
    for spec in args.match:
        if '=' not in spec:
            print(f"Error: --match '{spec}' must be WILDCARD=COLUMN", file=sys.stderr)
            sys.exit(1)
        wk, col = spec.split('=', 1)
        match_map[wk.strip()] = col.strip()

    # ---- read dry-run from stdin ------------------------------------------
    if sys.stdin.isatty():
        print("Waiting for `snakemake -n` output on stdin "
              "(pipe it in or redirect a saved log)...", file=sys.stderr)
    text = sys.stdin.read()

    # ---- parse wildcards --------------------------------------------------
    all_rules_mode = (args.rule == 'all')
    jobs = parse_wildcards_for_rule(text, args.rule)

    if not jobs:
        msg = ("No rule blocks with wildcards found in dry-run output.\n"
               "Check that stderr was captured with 2>&1."
               if all_rules_mode else
               f"No pending jobs found for rule '{args.rule}'.\n"
               f"Check the rule name and that stderr was captured with 2>&1.")
        print(msg, file=sys.stderr)
        sys.exit(0)

    rule_label = 'all rules' if all_rules_mode else f"rule '{args.rule}'"
    print(f"Found {len(jobs)} pending job(s) across {rule_label}.", file=sys.stderr)

    # ---- collect sample IDs ----------------------------------------------
    col_ids: dict[str, set] = defaultdict(set)
    for job in jobs:
        for wk, col in match_map.items():
            if wk in job:
                col_ids[col].add(job[wk])
            elif not all_rules_mode:
                print(f"Warning: wildcard '{wk}' not in job {job}", file=sys.stderr)

    # ---- default: print IDs ----------------------------------------------
    if args.samples is None:
        all_ids = sorted({sid for ids in col_ids.values() for sid in ids})
        print(f"{len(all_ids)} unique sample ID(s):", file=sys.stderr)
        for sid in all_ids:
            print(sid)
        return

    # ---- --samples mode: write filtered table ----------------------------
    try:
        import pandas as pd
    except ImportError:
        print("Error: pandas is required for --samples filtering.\n"
              "Install with: pip install pandas", file=sys.stderr)
        sys.exit(1)

    try:
        samples = pd.read_csv(args.samples, sep=args.sep, dtype=str)
    except FileNotFoundError:
        print(f"Error: samples table not found: {args.samples}", file=sys.stderr)
        sys.exit(1)

    missing_cols = [c for c in col_ids if c not in samples.columns]
    if missing_cols:
        print(f"Error: column(s) not in samples table: {missing_cols}\n"
              f"Available: {list(samples.columns)}", file=sys.stderr)
        sys.exit(1)

    mask = pd.Series(False, index=samples.index)
    for col, ids in col_ids.items():
        mask |= samples[col].isin(ids)

    filtered = samples[mask].reset_index(drop=True)
    n_in, n_out = len(samples), len(filtered)
    print(f"Samples table: {n_in} → {n_out} rows.", file=sys.stderr)

    if n_out == 0:
        print("Warning: filtered table is empty.  "
              "Check --match column names against the samples table.", file=sys.stderr)

    out_text = filtered.to_csv(sep=args.sep, index=False)
    if args.output:
        with open(args.output, 'w') as fh:
            fh.write(out_text)
        print(f"Written to: {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(out_text)


if __name__ == '__main__':
    main()
