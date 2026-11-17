import argparse
import pandas as pd
import numpy

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--in_maf',required=True,type=str,help='')
    parser.add_argument('--out_maf',required=True,type=str,help='')
    parser.add_argument('--blacklist',required=True,type=str,help='')
    parser.add_argument('--filter_cols',required=True,type=str,nargs='+',help='')
    parser.add_argument('--exac_freq',required=True,type=float,help='')

    return parser.parse_args()

def filter_maf(input, output, blacklist, filter_cols, exac_freq):
        # Load blacklist
        blacklist_pos = []
        with open(blacklist) as f:
            for line in f:
                if line.startswith("#"):
                    continue  # Ignore comment lines
                # Assuming the first column of this file is chrom:pos
                line = line.rstrip("\n").rstrip("r")
                pos = line.split("\t")[0]
                blacklist_pos.append(pos)
        blacklist_pos = set(blacklist_pos)

        # Load variants
        in_maf = pd.read_csv(input, sep ="\t", comment = "#")
        if in_maf.shape[0] == 0:  # Empty input MAF file
            in_maf.to_csv(output, sep="\t", header=True, index = False)
        else:
            # Filter based on allele frequencies
            # This is slow and inefficient, but is probably fine for this implementation
            filtered_maf = in_maf.copy()
            for f_col in filter_cols:
                filtered_maf = filtered_maf[(numpy.isnan(filtered_maf[f_col])) | (filtered_maf[f_col] < exac_freq)]

            # Filter out any variants at positions overlapping the blacklist
            blacklist_maf = filtered_maf
            filtered_maf["key"] = filtered_maf["Chromosome"].astype(str) + ":" + filtered_maf["Start_Position"].astype(str)
            blacklist_maf = blacklist_maf[filtered_maf["key"].isin(blacklist_pos) == False]

            blacklist_maf.to_csv(output, sep="\t", header=True, index = False)



def main():
    args = get_args()
    filter_maf(args.in_maf, args.out_maf, args.blacklist, args.filter_cols, args.exac_freq)

if __name__ == '__main__':
    main()