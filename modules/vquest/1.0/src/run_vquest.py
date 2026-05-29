#!/usr/bin/env python3
"""
Submit a per-chain FASTA to IMGT V-QUEST and write the AIRR TSV output.

The vquest package handles batching (50 sequences per request) and merges
results automatically. Requires outbound HTTPS access to www.imgt.org.
"""

import argparse
import sys

from vquest.vq import DEFAULTS, layer_configs, vquest
from vquest.util import VquestError


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta",         required=True,  help="input FASTA file")
    parser.add_argument("--species",       default="human")
    parser.add_argument("--receptor_type", default="IG",
                        help="receptorOrLocusType: IG or TR")
    parser.add_argument("--output",        required=True,  help="output AIRR TSV path")
    args = parser.parse_args()

    config = layer_configs(DEFAULTS, {
        "species":              args.species,
        "receptorOrLocusType":  args.receptor_type,
        "fileSequences":        args.fasta,
    })

    try:
        result = vquest(config)
    except VquestError as e:
        print(
            f"ERROR: V-QUEST returned an error for {args.fasta}:\n  {e}",
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

    with open(args.output, "w") as fh:
        fh.write(result["vquest_airr.tsv"])


if __name__ == "__main__":
    main()
