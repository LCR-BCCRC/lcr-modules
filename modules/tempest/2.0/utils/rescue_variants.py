"""A script that can rescue variants from VCF files
by looking at the QUAL, filters and AD

Designed to rescue variants labelled as germline but
probably aint.

Author KY
"""
import pandas as pd
import argparse
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', required=True, type=str, help='path to vcf to be parsed. Designed on SAGE v4.1 outputs')
    # filtering params
    parser.add_argument('--lower_af', required=False, type=float, default=0.05, help='lower bound of normal AF to consider for rescue')
    parser.add_argument('--upper_af', required=False, type=float, default=0.3, help='upper bound of normal AF to consider for rescue')
    parser.add_argument('--min_qual', required=False, type=int, default=20, help='minimum QUAL score to consider for rescue')
    parser.add_argument('--hotspots_only', required=False, type=bool, action='store_true', default=False, help='only rescue hotspot variants')

    parser.add_argument('--output', required=False, type=str, default='rescued_variants.tsv', help='output path for rescued variants tsv')
    return parser.parse_args()

class parse_vcf(object):
    """Parse a vcf file just how we like it"""

    def __init__(self, vcf_path):
        super(parse_vcf, self).__init__()
        self.vcf_path = vcf_path
        self.sample_id = None
        self.vcfDF = self.load_vcf()

    def load_vcf(self):
        vcf_names = self.get_vcf_names(self.vcf_path)
        # remove # symbol from names
        vcf_names = [name.lstrip('#') for name in vcf_names]
        vcf_df = pd.read_csv(self.vcf_path, sep='\t', comment='#', names=vcf_names, dtype=str)

        # move sample_id into a new matched_nromal and sample_id column, and rename the two sample columns to matched_normal, and tumour
        matched_normal= vcf_df.columns[-2]
        sample_id = vcf_df.columns[-1]

        self.sample_id = sample_id # add to object as its parsed

        # get column values
        normal_format_values= vcf_df[matched_normal].values
        sample_format_values= vcf_df[sample_id].values

        # extract tthe AD values, with are the third, delim by :
        normal_ad_values = [x.split(':')[2] if pd.notna(x) else './.' for x in normal_format_values]
        sample_ad_values = [x.split(':')[2] if pd.notna(x) else './.' for x in sample_format_values]

        # add as new columns
        vcf_df['n_ref_depth'] = [int(x.split(",")[0]) for x in normal_ad_values]
        vcf_df['n_alt_depth'] = [int(x.split(",")[1]) for x in normal_ad_values]
        vcf_df['t_ref_depth'] = [int(x.split(",")[0]) for x in sample_ad_values] 
        vcf_df['t_alt_depth'] = [int(x.split(",")[1]) for x in sample_ad_values]

        vcf_df['t_depth'] = vcf_df['t_ref_depth'] + vcf_df['t_alt_depth']
        vcf_df['n_depth'] = vcf_df['n_ref_depth'] + vcf_df['n_alt_depth']

        vcf_df['t_AF'] = vcf_df['t_alt_depth'] / vcf_df['t_depth']
        vcf_df['n_AF'] = vcf_df['n_alt_depth'] / vcf_df['n_depth']

        # rename columns
        vcf_df.rename(columns={sample_id: "sample_id", matched_normal: "matched_normal"}, inplace=True)
        vcf_df["matched_normal"] = matched_normal
        vcf_df["sample_id"] = sample_id.strip()

        # label hotspots
        vcf_df['is_hotspot'] = vcf_df.apply(self.find_hotspots, axis=1)
        # force QUAL col into numeric
        vcf_df["QUAL"] = pd.to_numeric(vcf_df["QUAL"], errors='coerce')

        return vcf_df

    def get_vcf_names(self, vcf_path):
        with open(vcf_path, "rt") as ifile:
            for line in ifile:
                if line.startswith("#CHROM"):
                    vcf_names = [x for x in line.split('\t')]
                    break
        ifile.close()
        return vcf_names

    def find_hotspots(self, row):
        info_fields = row['INFO'].split(';')
        for field in info_fields:
            if field.startswith('TIER='):
                tier_value = field.split('=')[1]
                if tier_value == 'HOTSPOT':
                    return True
        return False

def rescue_variants(indf: pd.DataFrame, lower_af: float, upper_af: float, min_qual: int = 20) -> pd.DataFrame:

    # get candidate variants based on VAF bounds and min qual
    candidates = indf.loc[(indf["QUAL"] > min_qual) & (indf["n_AF"] > lower_af)  & (indf["n_AF"] < upper_af)]
    # if a variants FILTER column only contains keeper filts, keep variant
    rescued_vars = candidates[candidates.apply(check_filters, axis=1)]

    return rescued_vars

def check_filters(row):
    keeper_filters = ["maxGermlineVAF", "maxGermlineRelQual"]
    filters = row["FILTER"].split(";")

    if any(filt not in keeper_filters for filt in filters):
        return False
    else:
        return True

    for filt in filters:
        if filt not in keeper_filters:
            return False
    return True


def main():
    args = get_args()
    print("Parsing vcf...")
    vcf_parser = parse_vcf(args.vcf)
    # filter for vars to rescue
    print("Rescuging variants...")
    rescued_variants = rescue_variants(vcf_parser.vcfDF, args.lower_af, args.upper_af, args.min_qual)
    print(f"""Found {str(rescued_variants.shape[0])} variants to rescue""")
    # output rescued variants
    print("Writing output...")
    if args.output:
        rescued_variants.to_csv(os.path.join(args.output, f"{vcf_parser.sample_id}_rescued_variants.tsv"), 
                        sep='\t', index=False)
    else:
        basename = os.path.basename(args.vcf)
        outpath = args.vcf.replace(basename , vcf_parser.sample_id ,"_rescued_variants.tsv")
        rescued_variants.to_csv(outpath, sep='\t', index=False)



if __name__ == "__main__":
    main()