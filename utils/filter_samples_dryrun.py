#!/usr/bin/env python3
"""
Summarise pending Snakemake jobs as a tumour/normal pair × module table.

Reads `snakemake -n` output from stdin and writes a TSV with one row per
(tumour, normal, module) combination that has at least one pending job:

    tumour_id   normal_id   module   first_incomplete_rule   target_file

tumour_id             Tumour sample wildcard value.
normal_id             Normal sample wildcard value.
module                Module name extracted from the output path
                      (results/<module>-<version>/...).
first_incomplete_rule Most upstream pending rule for this pair × module,
                      in the order Snakemake would run them.
target_files          Space-separated list of all terminal output files
                      (paths containing '99-outputs') for this pair × module.
                      Directly usable as Snakemake targets.

Usage
-----
    ./demo/dry-run.sh genome_Snakefile.smk all "" runtime_config.yaml 2>&1 \\
        | utils/filter_samples_dryrun.py

Filter to a specific rule
--------------------------
    ... | utils/filter_samples_dryrun.py --rule _run_battenberg

Write to a file
---------------
    ... | utils/filter_samples_dryrun.py --output pending_pairs.tsv

Re-run only the pending pairs
------------------------------
    # Inspect the table first
    column -t pending_pairs.tsv

    # target_files column is space-separated and directly usable as targets
    snakemake $(cut -f5 pending_pairs.tsv | tail -n +2) \\
        -s genome_Snakefile.smk --use-singularity -j 32 ...

Arguments
---------
--rule / -r     Restrict to jobs from a specific rule (default: all rules).
--match / -m    Wildcard→samples-column mapping as WILDCARD=COLUMN.
                Default: tumour_id=tumour_id, normal_id=normal_id.
                Override if your wildcards use different names.
--output / -o   Output file path (default: stdout).
--sep           Output field separator (default: tab).
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_wildcard_str(wc_str: str) -> dict:
    """Parse 'key=val, key2=val2' into a dict."""
    result: dict[str, str] = {}
    for pair in re.split(r',\s+', wc_str.strip()):
        if '=' in pair:
            k, v = pair.split('=', 1)
            result[k.strip()] = v.strip()
    return result


def _parse_output_files(block: str) -> list[str]:
    """
    Extract output file paths from a job block.

    The output section looks like:
        output: path1, path2,
                path3
    terminated by the next field (log:, jobid:, wildcards:, etc.).
    """
    # Match 'output:' followed by everything until the next unindented field
    m = re.search(
        r'^\s+output:\s+(.+?)(?=^\s+\w+:|$)',
        block,
        re.MULTILINE | re.DOTALL,
    )
    if not m:
        return []
    raw = m.group(1)
    # Strip trailing whitespace/newlines, collapse internal whitespace
    raw = re.sub(r'\s+', ' ', raw).strip()
    # Split on ', ' — paths in Snakemake output don't contain ', '
    paths = [p.strip() for p in raw.split(', ') if p.strip()]
    return paths


def _module_from_paths(paths: list[str]) -> str:
    """
    Extract the module name from lcr-modules output paths.
    Expects paths like  results/<module>-<version>/...
    Returns '' if none match.
    """
    for path in paths:
        m = re.search(r'results/([^/\-]+)-[\d.]+/', path)
        if m:
            return m.group(1)
    return ''


def _terminal_output(paths: list[str]) -> str:
    """
    Return the first path that contains '99-outputs', which is the
    lcr-modules convention for final deliverable files.
    Falls back to the first path if none qualify.
    """
    for p in paths:
        if '99-outputs' in p:
            return p
    return paths[0] if paths else ''


def parse_jobs(text: str, rule_filter: str | None) -> list[dict]:
    """
    Return a list of job dicts parsed from `snakemake -n` output.

    Each dict has keys: rule, wildcards (dict), outputs (list), module,
    terminal_output.
    """
    any_rule_header = re.compile(
        r'(?:^\[.*?\]\s+)?rule\s+(\w+)\s*:',
        re.MULTILINE,
    )
    wildcards_pat = re.compile(r'^\s+wildcards:\s+(.+)$', re.MULTILINE)

    jobs = []
    for m in any_rule_header.finditer(text):
        rule_name = m.group(1)
        if rule_filter and rule_name != rule_filter:
            continue

        block_start = m.start()
        next_m = any_rule_header.search(text, m.end())
        block = text[block_start: next_m.start() if next_m else len(text)]

        wc_m = wildcards_pat.search(block)
        if not wc_m:
            continue
        wildcards = _parse_wildcard_str(wc_m.group(1))

        outputs = _parse_output_files(block)
        module = _module_from_paths(outputs)
        terminal = _terminal_output(outputs)

        jobs.append({
            'rule': rule_name,
            'wildcards': wildcards,
            'outputs': outputs,
            'module': module,
            'terminal_output': terminal,
        })

    return jobs


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def build_pair_table(
    jobs: list[dict],
    tumour_wc: str = 'tumour_id',
    normal_wc: str = 'normal_id',
) -> list[dict]:
    """
    Group jobs by (tumour_id, normal_id, module).

    Returns rows sorted by (tumour_id, normal_id, module), each with:
        tumour_id, normal_id, module, first_incomplete_rule, target_file
    """
    # Key: (tumour_id, normal_id, module)
    # Value: list of jobs in dry-run order (already topological)
    groups: dict[tuple, list[dict]] = defaultdict(list)

    for job in jobs:
        wc = job['wildcards']
        tumour = wc.get(tumour_wc, '')
        normal = wc.get(normal_wc, '')
        module = job['module']

        if not tumour or not normal:
            continue  # skip rules without both wildcards

        groups[(tumour, normal, module)].append(job)

    rows = []
    for (tumour, normal, module), grp_jobs in sorted(groups.items()):
        first_rule = grp_jobs[0]['rule']

        # Collect all unique 99-outputs paths across every job in this group.
        # Using a list to preserve dry-run order while deduplicating.
        seen: set[str] = set()
        targets: list[str] = []
        for j in grp_jobs:
            for p in j['outputs']:
                if '99-outputs' in p and p not in seen:
                    seen.add(p)
                    targets.append(p)

        # Fall back to the last job's first output if nothing in 99-outputs
        if not targets:
            fb = grp_jobs[-1]['terminal_output']
            if fb:
                targets = [fb]

        rows.append({
            'tumour_id': tumour,
            'normal_id': normal,
            'module': module,
            'first_incomplete_rule': first_rule,
            'target_files': ' '.join(targets),
        })

    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--rule', '-r', default=None,
                   help='Restrict output to jobs from a specific rule name')
    p.add_argument('--tumour-wc', default='tumour_id',
                   help='Wildcard name for tumour ID (default: tumour_id)')
    p.add_argument('--normal-wc', default='normal_id',
                   help='Wildcard name for normal ID (default: normal_id)')
    p.add_argument('--output', '-o', default=None,
                   help='Output file path (default: stdout)')
    p.add_argument('--sep', default='\t',
                   help='Output field separator (default: tab)')
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    if sys.stdin.isatty():
        print("Waiting for `snakemake -n` output on stdin "
              "(pipe it in or redirect a saved log)...", file=sys.stderr)
    text = sys.stdin.read()

    jobs = parse_jobs(text, rule_filter=args.rule)

    if not jobs:
        msg = (f"No pending jobs found for rule '{args.rule}'."
               if args.rule else
               "No rule blocks with wildcards found in dry-run output.")
        print(msg + "\nCheck that stderr was captured with 2>&1.", file=sys.stderr)
        sys.exit(0)

    print(f"Found {len(jobs)} pending job(s).", file=sys.stderr)

    rows = build_pair_table(jobs, args.tumour_wc, args.normal_wc)

    if not rows:
        print("No tumour/normal pairs identified.\n"
              "If your wildcards use different names, set --tumour-wc / --normal-wc.",
              file=sys.stderr)
        sys.exit(0)

    print(f"{len(rows)} pending pair × module combination(s).", file=sys.stderr)

    sep = args.sep
    header = sep.join(['tumour_id', 'normal_id', 'module',
                       'first_incomplete_rule', 'target_files'])
    lines = [header]
    for row in rows:
        lines.append(sep.join([
            row['tumour_id'],
            row['normal_id'],
            row['module'],
            row['first_incomplete_rule'],
            row['target_files'],
        ]))
    out_text = '\n'.join(lines) + '\n'

    if args.output:
        with open(args.output, 'w') as fh:
            fh.write(out_text)
        print(f"Written to: {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(out_text)


if __name__ == '__main__':
    main()
