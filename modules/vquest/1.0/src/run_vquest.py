#!/usr/bin/env python3
"""
Submit a per-chain FASTA to IMGT V-QUEST and write the AIRR TSV output.

The vquest package handles batching (50 sequences per request) and merges
results automatically. Requires outbound HTTPS access to www.imgt.org.
"""

import argparse
import re
import sys
import tempfile

import requests as _requests

from vquest.vq import DEFAULTS, layer_configs, vquest
from vquest.util import VquestError

_VALID_NT = re.compile(r"[^ACGTN]", re.IGNORECASE)


def _sanitize_fasta(src_path):
    """
    Write a sanitized copy of src_path to a temp file and return its path.
    All nucleotide sequences are uppercased and any character that is not
    A/C/G/T/N is replaced with N. V-QUEST rejects the entire batch of 50
    when even one sequence contains an invalid character, so cleaning here
    prevents unnecessary whole-batch failures.
    """
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    )
    with open(src_path) as fh:
        for line in fh:
            if line.startswith(">"):
                tmp.write(line)
            else:
                cleaned = _VALID_NT.sub("N", line.rstrip().upper())
                tmp.write(cleaned + "\n")
    tmp.close()
    return tmp.name


def _check_html_errors(text):
    """
    Raise RuntimeError if text looks like a V-QUEST HTML error page.
    V-QUEST returns errors inside <ul class="errorMessage"> elements,
    which the vquest package does not detect (it only checks div.form_error).
    """
    if not text.lstrip().startswith("<!DOCTYPE"):
        return
    msgs = re.findall(
        r'<ul[^>]+class=["\']errorMessage["\'][^>]*>.*?<span>(.*?)</span>',
        text, re.DOTALL
    )
    if msgs:
        raise RuntimeError(
            "V-QUEST returned HTML error page(s):\n"
            + "\n".join(f"  {m.strip()}" for m in msgs[:5])
            + (f"\n  ... ({len(msgs) - 5} more)" if len(msgs) > 5 else "")
        )
    raise RuntimeError(
        "V-QUEST returned an HTML page instead of AIRR data.\n"
        f"First 500 chars: {text[:500]}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta",           required=True,  help="input FASTA file")
    parser.add_argument("--species",         default="human")
    parser.add_argument("--receptor_type",   default="IG",
                        help="receptorOrLocusType: IG or TR")
    parser.add_argument("--molecule_type",   default="cDNA",
                        help="V-QUEST moleculeType: cDNA (default), gDNA, or Unknown")
    parser.add_argument("--output",          required=True,  help="output AIRR TSV path")
    parser.add_argument("--request_timeout", type=int, default=120,
                        help="per-request timeout in seconds for IMGT V-QUEST HTTP calls "
                             "(default: 120). Prevents hung connections from stalling the job.")
    args = parser.parse_args()

    # Patch requests.post so every call to www.imgt.org has a bounded timeout.
    # The vquest library calls requests.post with no timeout, which can hang
    # indefinitely when the IMGT server is unresponsive.
    _orig_post = _requests.post

    def _post_with_timeout(*pargs, **kwargs):
        kwargs.setdefault("timeout", args.request_timeout)
        return _orig_post(*pargs, **kwargs)

    _requests.post = _post_with_timeout

    clean_fasta = _sanitize_fasta(args.fasta)

    config = layer_configs(DEFAULTS, {
        "species":              args.species,
        "receptorOrLocusType":  args.receptor_type,
        "moleculeType":         args.molecule_type,
        "fileSequences":        clean_fasta,
    })

    try:
        result = vquest(config)
    except VquestError as e:
        # e.args[0] is "; ".join(errors), e.args[1] is the raw errors list.
        # Either may be empty if vquest failed to parse the error HTML.
        print(
            f"ERROR: V-QUEST returned an error for {args.fasta}:\n"
            f"  message : {e.args[0]!r}\n"
            f"  errors  : {e.args[1] if len(e.args) > 1 else '(not available)'}",
            file=sys.stderr,
        )
        sys.exit(1)

    if "vquest_airr.tsv" not in result:
        print(
            f"ERROR: vquest_airr.tsv not found in V-QUEST response.\n"
            f"Available keys: {list(result.keys())}",
            file=sys.stderr,
        )
        sys.exit(1)

    airr_text = result["vquest_airr.tsv"]

    try:
        _check_html_errors(airr_text)
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    with open(args.output, "w") as fh:
        fh.write(airr_text)


if __name__ == "__main__":
    main()
