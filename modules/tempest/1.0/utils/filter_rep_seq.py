import pyfaidx
import argparse


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
    '--reference_genome',
    required=True,
    type=str,
    help='')

    parser.add_argument(
    '--in_maf',
    required=True,
    type=str,
    help='')

    parser.add_argument(
    '--out_maf',
    required=True,
    type=str,
    help='')

    parser.add_argument(
    '--max_repeat_len',
    required=True,
    type=str,
    help='')

    return parser.parse_args()

def filter_repeat(input, output, refseq, max_repeats):
    with open(input) as f, open(output, "w") as o:
        for line in f:
            if line.startswith("#") or line.startswith("Hugo_Symbol"):
                # Skip header lines.
                o.write(line)
                continue
            cols = line.split("\t")
            # Get alternate allele
            alt = cols[12]
            if alt == "-":
                alt = cols[10]
            if len(alt) > 1:  # Handle large indels
                alt = alt[0]
            # Get position of variant.
            
            pos = int(cols[6]) - 1  # Offset by 1 for MAF vs pyfaidx sequence.
            chrom = cols[4]

            # Determine length of repeat via reference sequence.
            repeat_len = 0
            while True:
                pos += 1
                base = refseq[chrom][pos]
                if base != alt or repeat_len > max_repeats:
                    break
                repeat_len += 1

            if repeat_len < max_repeats:
                # Not a repeat (or not sufficiently long)
                o.write(line)

def main():
    args = get_args()
    refseq = pyfaidx.Fasta(args.reference_genome)
    filter_repeat(args.in_maf, args.out_maf, refseq, int(args.max_repeat_len))

if __name__ == '__main__':
    main()



