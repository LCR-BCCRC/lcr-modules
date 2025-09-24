"""
This script is used to add variants complete with read depths found in additional mafs to an index (original) maf.

The code looks for variants in using start position, chromosome, and alt allele that are found in a list of
additional maf files but not in the index maf file. If found, the script will update the index maf
file with the new variants with allele counts > 0 from the index bam file. If no new variants are found
it will return the index maf file as is.

This script is a derivitive of the "augment_ssm" by the Morin Lab. 
Key features added here are CLI support, and packaing of the code into a class that
has been refactored for readability and maintainability.

The core logic on what the script does is the same, just
in some prettier syntax and more documentation.

Refactor notes:
- UMI collection is optional (enabled with --compute_umi_metrics).
- Variant-type helpers remain, and insertions/deletions are single-pass loops that collect UMI tags inline.
- UMI family size tag is cD; UMI/group identifier tag is MI (fallback to read name if MI missing).
"""

import pandas as pd
from pathlib import Path 
import pysam
import os
import multiprocessing
import argparse
import shutil

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_id',required=True,type=str,help='tumour sample id')
    parser.add_argument('--threads',required=False,type=int,default=24,help='number of threads to use, default is 24')
    parser.add_argument('--index_maf',required=True,type=str,help='index maf file')
    parser.add_argument('--index_bam',required=True,type=str,help='index bam file')
    parser.add_argument('--add_maf_files',required=False,nargs='*',type=str,help='additional maf files')
    parser.add_argument('--genome_build',required=True,type=str,help='genome build')
    parser.add_argument('--alt_count_min',required=True,type=int,help='minimum alt count reads needed to keep variant')
    parser.add_argument('--output',required=True,type=str,help='output file')

    # Optional: compute UMI metrics (MI and cD tags)
    parser.add_argument('--compute_umi_metrics', action='store_true',
                        help='If set, compute UMI_mean, UMI_max, UMI_3_count for alt-supporting reads (uses MI and cD tags)')
    parser.add_argument("--min_UMI_3_count", type=int, default=1,
                        help="Filters UMI_3_count metric. Which is a count of the number of reads that have UMI family sizes of at least 3.")
    return parser.parse_args()


class AugmentMAF(object):
    """Add variants complete with read depths found in additional mafs to an index maf.

    Useful if you have multiple samples from the same patient and want to see if there is supporting
    reads for a variant in one sample that is not called in another sample.
    """

    def __init__(self, sample_id: str, index_maf:str,
                        index_bam: str, add_maf_files: list, genome_build: str,
                        output: str , threads: int = 6, min_alt_count: int= 3,
                        compute_umi_metrics: bool = False,
                        min_UMI_3_count: int = None):
        super(AugmentMAF, self).__init__()
        # user provided inputs
        self.sample_id = sample_id
        self.threads = threads
        self.index_maf = self.read_maf(index_maf)
        self.index_bam = index_bam
        self.add_maf_files = add_maf_files or []
        self.genome_build = genome_build
        self.output = output
        self.min_alt_count = min_alt_count
        self.compute_umi_metrics = compute_umi_metrics  # optional UMI collection
        self.min_UMI_3_count = min_UMI_3_count 

        # computed variables
        self.chromosomes = self.get_genome_chromosomes()
        self.master_maf = self.run_or_not()
        self.augmented_vaf = self.augment_maf() # the magix

    def read_maf(self, maf_file: str) -> pd.DataFrame:
        return pd.read_csv(maf_file, sep="\t", comment='#',
        dtype={'Chromosome':'string','t_ref_count':'string','t_alt_count':'string','t_depth':'string', "Tumor_Sample_Barcode":"string"},
        low_memory=False)

    def get_genome_chromosomes(self)-> list:
        """Return a list of chromosomes for the genome build.
        """
        if self.genome_build.lower() in ["hg38","hg19-reddy", "grch38"]:
            return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        else:
            return [str(i) for i in range(1, 23)] + ["X", "Y"]

    def write_output(self, output_file: str = None) -> None:
        """Write the augmented maf to a file.

        Args:
            output_file (str): output file name
        """
        if output_file:
            out_name = output_file
        elif self.output:
            out_name = self.output
        else:
            raise print("No output file name provided ...")

        if self.min_UMI_3_count and self.augmented_vaf.shape[0] > 0:
            self.augmented_vaf = self.augmented_vaf[self.augmented_vaf["UMI_3_count"].astype(int) >= self.min_UMI_3_count].copy()

        self.augmented_vaf.to_csv(out_name, sep="\t", index=False)

    def augment_maf(self):
        """Loop through each chromosome in parallel
        and check for variants in the index bam from
        the add maf files. If found, update the index
        maf file with the new variant information.
        """
        # test if master_maf is a df or not
        if self.master_maf is None:
            print("No additional maf files provided ... returning index maf")
            # add standard column anyways
            self.index_maf["origin_samples"] = ""
            return self.index_maf

        # if there are no new variants to add, return the index maf
        new_vars = self.master_maf.loc[self.master_maf["variant_source"] == "additional_maf"].shape[0]
        if new_vars == 0:
            print("No new variants to add to index maf ... returning index maf")
            # add standard column anyways
            self.index_maf["origin_samples"] = ""
            return self.index_maf
        elif new_vars > 0:
            print(f"Found {new_vars} new variants to add to index maf ...")

        # run in parallel for all chromosomes
        cool_pool = multiprocessing.Pool(processes = self.threads)
        results = cool_pool.starmap(self._subset_and_run, [[chrm] for chrm in self.chromosomes])
        cool_pool.close()
        cool_pool.join()

        # concatenate mafs and write them out to a file
        augmented_maf_merged = pd.concat(results)

        # add autmented vaf to index maf, reset sample_id in all rows, and remove duplicates
        full_augmented_maf_merged = pd.concat([self.index_maf, augmented_maf_merged])
        full_augmented_maf_merged["Tumor_Sample_Barcode"] = self.sample_id
        full_augmented_maf_merged = full_augmented_maf_merged.drop_duplicates().reset_index(drop=True)

        # drop key column
        if "key" in full_augmented_maf_merged.columns:
            full_augmented_maf_merged.drop(columns=["key"], inplace=True)
        # recalculate AF
        full_augmented_maf_merged = full_augmented_maf_merged.dropna(subset=["Hugo_Symbol","t_alt_count","t_depth"]).copy()
        full_augmented_maf_merged["AF"] = full_augmented_maf_merged["t_alt_count"].astype(int) / full_augmented_maf_merged["t_depth"].astype(int)
        # insert AF column at 43 column
        full_augmented_maf_merged.insert(43, "AF", full_augmented_maf_merged.pop("AF"))

        return full_augmented_maf_merged

    def _get_missing_maf_rows(self):
        """Compare other maf to index maf and 
        return rows that are missing from the index maf
        based on chromosome and start position of variants.
        """
        additional_mafs = []
        for maf in self.add_maf_files:
            other_maf = self.read_maf(maf)
            other_maf["variant_source"] = "additional_maf"
            additional_mafs.append(other_maf)
        # add column saying what samples vars came from
        additional_vars = self.mark_origin(pd.concat(additional_mafs))

        # concatenate all mafs with index maf and drop duplicate keys, but keep index maf version
        all_mafs = pd.concat([self.index_maf, additional_vars]).reset_index(drop=True).copy()
        # add variant key
        all_mafs["variant_key"] = all_mafs["Chromosome"] + "_" + all_mafs["Start_Position"].astype(str) + "_" + all_mafs["Tumor_Seq_Allele2"].astype(str)
        # keep only the first instance of a variant, so if its in the index maf, keep that one
        all_mafs = all_mafs.drop_duplicates(subset="variant_key", keep="first").copy()
        all_mafs.drop(columns=["variant_key"], inplace=True)

        return all_mafs

    def run_or_not(self):
        """Check if there are additional maf files to work with."""
        if len(self.add_maf_files) > 0:
            return self._get_missing_maf_rows()
        else:
            return None

    def mark_origin(self, add_mafs: pd.DataFrame) -> pd.DataFrame:
        """Mark the origin sample(s) of each variant if it is from an additional maf file."""
        add_mafs["variant_key"] = add_mafs["Chromosome"] + "_" + add_mafs["Start_Position"].astype(str) + "_" + add_mafs["Tumor_Seq_Allele2"].astype(str)
        variant_key_dict = add_mafs.groupby("variant_key")["Tumor_Sample_Barcode"].agg(lambda x: ",".join(list(x))).to_dict()
        add_mafs["origin_samples"] = add_mafs["variant_key"].map(variant_key_dict)
        return add_mafs

    def _subset_and_run(self, this_chromosome: str) -> pd.DataFrame:
        """Subset variants per chr and find support in index sample bam file."""
        print(f"Working on chromosome {this_chromosome} ...")
        missing_maf = self.master_maf[self.master_maf["Chromosome"] == this_chromosome].copy()
        # find only rows with the source "additional_maf" aka were not called in the index maf
        missing_maf = missing_maf[missing_maf["variant_source"] == "additional_maf"].copy()

        print(f"Found these variants to augment: {missing_maf}")

        for var in missing_maf.itertuples():
            read_counts = self._get_read_support(chrom=var.Chromosome,
                            start= int(var.Start_Position),
                            end=int(var.End_Position),
                            ref_base=var.Reference_Allele,
                            alt_base=var.Tumor_Seq_Allele2,
                            bamfile=self.index_bam,
                            mut_class=var.Variant_Type)

            # write counts (as strings)
            missing_maf.at[var.Index, "t_depth"] = read_counts["tumor_depth"]
            missing_maf.at[var.Index, "t_ref_count"] = read_counts["tumor_ref_count"]
            missing_maf.at[var.Index, "t_alt_count"] = read_counts["tumor_alt_count"]

            # optional UMI columns
            if self.compute_umi_metrics:
                missing_maf.at[var.Index, "UMI_mean"] = read_counts.get("UMI_mean", "")
                missing_maf.at[var.Index, "UMI_max"] = read_counts.get("UMI_max", "")
                missing_maf.at[var.Index, "UMI_3_count"] = read_counts.get("UMI_3_count", "")

        print(f"Finished working on chromosome {this_chromosome} ...")
        # only return variants with a t_alt_count > min_alt_count and depth > min_alt_count
        missing_maf = missing_maf[(missing_maf["t_alt_count"].astype(int) > self.min_alt_count) & (missing_maf["t_depth"].astype(int) > self.min_alt_count)].copy()
        return missing_maf

    def _get_read_support(self, chrom: str, start: int, end: int, ref_base: str, alt_base: str, 
                        bamfile: str, mut_class: str, min_base_qual=10, min_mapping_qual=1,
                        padding: int = 100) -> dict:
        """Main coordinator function for getting read support for a variant."""
        
        # Setup common parameters
        params = self._setup_variant_params(start, end, ref_base, mut_class, padding)
        
        # Open BAM file
        samfile = self._open_bam_file(bamfile)
        
        # Process variant based on type
        if mut_class == "INS":
            ref_count, alt_count, total_depth, umi_sizes = self._process_insertion(
                samfile, chrom, params, min_base_qual, min_mapping_qual, self.compute_umi_metrics)
        elif mut_class == "DEL":
            ref_count, alt_count, total_depth, umi_sizes = self._process_deletion(
                samfile, chrom, params, ref_base, min_base_qual, min_mapping_qual, self.compute_umi_metrics)
        elif mut_class in ["SNP", "DNP"]:
            ref_count, alt_count, total_depth, umi_sizes = self._process_snp_dnp(
                samfile, chrom, params, ref_base, alt_base, mut_class, 
                min_base_qual, min_mapping_qual, self.compute_umi_metrics)
        else:
            samfile.close()
            raise ValueError(f"Unsupported mutation class: {mut_class}")
        
        samfile.close()
        
        result = self._format_results(ref_count, alt_count, total_depth)

        # Optional UMI summaries (as strings to match your existing string outputs)
        if self.compute_umi_metrics:
            umi_mean, umi_max, umi_ge3 = self._compute_umi_metrics(umi_sizes)
            result["UMI_mean"] = "" if umi_mean is None else f"{umi_mean:.3f}"
            result["UMI_max"] = "" if umi_max is None else str(umi_max)
            result["UMI_3_count"] = "" if umi_ge3 is None else str(umi_ge3)

        return result

    def _setup_variant_params(self, start: int, end: int, ref_base: str, mut_class: str, padding: int) -> dict:
        """Setup common parameters for variant processing."""
        actual_start = start - 1  # Convert to 0-based indexing
        actual_end = end
        
        # Calculate padded region
        padded_start = start - padding
        padded_end = end + padding
        
        params = {
            'original_start': start,
            'actual_start': actual_start,
            'actual_end': actual_end,
            'padded_start': padded_start,
            'padded_end': padded_end,
            'variant_length': len(ref_base)
        }
        
        if mut_class == "DEL":
            params['del_positions'] = list(range(actual_start, end))
        
        return params

    def _open_bam_file(self, bamfile: str):
        """Open BAM or CRAM file."""
        realname = str(Path(bamfile).resolve())
        if realname.endswith("cram"):
            return pysam.AlignmentFile(realname, 'rc')
        else:
            return pysam.AlignmentFile(bamfile, 'rb')

    # ---------- Refactored helpers with optional UMI collection ----------

    def _process_insertion(self, samfile, chrom: str, params: dict, min_base_qual: int, min_mapping_qual: int,
                           collect_umi: bool) -> tuple:
        """Process insertion variants at the anchor column in a single pass.
        Counting semantics identical to your prior approach; UMI collection is optional.
        """
        ref_count = alt_count = total_depth = 0
        umi_sizes = []
        seen = set()  # MI/readname dedup

        for col in samfile.pileup(
            chrom, params['padded_start'], params['padded_end'],
            ignore_overlaps=False, min_base_quality=min_base_qual,
            min_mapping_quality=min_mapping_qual, truncate=True
        ):
            if col.pos != params['actual_start']:
                continue

            for pu in col.pileups:
                if pu.is_refskip:
                    continue
                total_depth += 1
                if pu.indel and pu.indel > 0:  # insertion starts here
                    alt_count += 1
                    if collect_umi:
                        key, size = self._umi_key_and_size(pu.alignment)
                        if size is not None and key not in seen:
                            seen.add(key)
                            umi_sizes.append(size)
                else:
                    ref_count += 1

        return ref_count, alt_count, total_depth, umi_sizes

    def _process_deletion(self, samfile, chrom: str, params: dict, ref_base: str, 
                          min_base_qual: int, min_mapping_qual: int, collect_umi: bool) -> tuple:
        """Process deletion variants across the span in a single pass.
        - alt_count is the maximum number of deletion-supporting reads ('is_del') across the span.
        - total_depth is the maximum column depth across the span.
        - UMI sizes are taken from the column that determines alt_count (to mirror your counting).
        """
        max_depth = 0
        max_alt = 0
        umi_at_max_alt = []

        del_start = params['del_positions'][0] if params.get('del_positions') else params['actual_start']
        del_end = params['del_positions'][-1] + 1 if params.get('del_positions') else params['actual_end']

        for col in samfile.pileup(
            chrom, del_start, del_end,
            ignore_overlaps=False, min_base_quality=min_base_qual,
            min_mapping_quality=min_mapping_qual, truncate=True
        ):
            if not (del_start <= col.pos < del_end):
                continue

            col_depth = 0
            col_alt = 0
            col_seen = set()
            col_umis = []

            for pu in col.pileups:
                if pu.is_refskip:
                    continue
                col_depth += 1
                if pu.is_del:
                    col_alt += 1
                    if collect_umi:
                        key, size = self._umi_key_and_size(pu.alignment)
                        if size is not None and key not in col_seen:
                            col_seen.add(key)
                            col_umis.append(size)

            if col_depth > max_depth:
                max_depth = col_depth
            if col_alt > max_alt:
                max_alt = col_alt
                umi_at_max_alt = col_umis  # align UMIs with the column that set alt_count

        ref_count = max_depth - max_alt
        alt_count = max_alt
        total_depth = max_depth

        return ref_count, alt_count, total_depth, umi_at_max_alt

    def _process_snp_dnp(self, samfile, chrom: str, params: dict, ref_base: str, alt_base: str, 
                        mut_class: str, min_base_qual: int, min_mapping_qual: int,
                        collect_umi: bool) -> tuple:
        """
         SNP/DNP processing:
        - One pileup pass to collect per-fragment observations (up to 2) and remember an alt-bearing alignment.
        - One small finalize pass to apply fragment-level counting semantics and collect UMIs for alt fragments.
        """
        per_read = {}  # read_name -> {"bases": [base1, base2], "alt_aln": AlignedSegment | None}

        for pileupcolumn in samfile.pileup(
            chrom, params['padded_start'], params['padded_end'],
            ignore_overlaps=False, min_base_quality=min_base_qual,
            min_mapping_quality=min_mapping_qual, truncate=True
        ):
            if pileupcolumn.pos != params['actual_start']:
                continue

            for pu in pileupcolumn.pileups:
                if pu.is_del or pu.is_refskip:
                    continue

                base = self._extract_variant_sequence(pu, mut_class, params['variant_length'])
                if not base:
                    continue

                read_name = pu.alignment.query_name
                rec = per_read.get(read_name)
                if rec is None:
                    rec = {"bases": [], "alt_aln": None}
                    per_read[read_name] = rec

                if len(rec["bases"]) < 2:
                    rec["bases"].append(base)

                if collect_umi and base == alt_base and rec["alt_aln"] is None:
                    rec["alt_aln"] = pu.alignment

        # Finalize counts and collect UMI sizes only for alt fragments
        ref_count = alt_count = total_depth = 0
        umi_sizes = []
        seen_keys = set()

        for read_name, rec in per_read.items():
            bases = rec["bases"]
            total_depth += 1

            if len(bases) == 1:
                b = bases[0]
                if b == ref_base:
                    ref_count += 1
                elif b == alt_base:
                    alt_count += 1
                    if collect_umi and rec["alt_aln"] is not None:
                        key, size = self._umi_key_and_size(rec["alt_aln"])
                        if size is not None and key not in seen_keys:
                            seen_keys.add(key)
                            umi_sizes.append(size)

            else:
                # Two observations
                if bases[0] == bases[1]:
                    b = bases[0]
                    if b == ref_base:
                        ref_count += 1
                    elif b == alt_base:
                        alt_count += 1
                        if collect_umi and rec["alt_aln"] is not None:
                            key, size = self._umi_key_and_size(rec["alt_aln"])
                            if size is not None and key not in seen_keys:
                                seen_keys.add(key)
                                umi_sizes.append(size)
                else:
                    # Disagree: count as ref if any base is ref; never alt (matches your logic)
                    if ref_base in bases:
                        ref_count += 1
                    # else: neither

        return ref_count, alt_count, total_depth, umi_sizes

    def _umi_key_and_size(self, aln):
        """Return (dedup_key, family_size) where key is MI tag if present else read name; family size from cD."""
        try:
            key = aln.get_tag("MI")
        except KeyError:
            key = aln.query_name
        try:
            size = int(aln.get_tag("cD"))
        except KeyError:
            size = None
        return key, size

    def _compute_umi_metrics(self, sizes):
        if not sizes:
            return None, None, None
        mean_v = sum(sizes) / len(sizes)
        max_v = max(sizes)
        ge3 = sum(1 for x in sizes if x >= 3)
        return mean_v, max_v, ge3

    def _extract_variant_sequence(self, pileupread, mut_class: str, variant_length: int) -> str:
        """Extract the appropriate sequence from a read based on variant type."""
        query_pos = pileupread.query_position
        query_sequence = pileupread.alignment.query_sequence
        
        if query_pos is None:
            return None

        if mut_class == "SNP":
            if 0 <= query_pos < len(query_sequence):
                return query_sequence[query_pos]
            return None
        elif mut_class == "DNP":
            if query_pos + variant_length <= len(query_sequence):
                return query_sequence[query_pos:query_pos + variant_length]
            else:
                return None
        else:
            return None

    def _count_read_support(self, read_bases: dict, ref_base: str, alt_base: str) -> tuple:
        """Count read support from collected read bases (fragment-level aggregation)."""
        ref_count = alt_count = total_depth = 0
        
        for read_name in read_bases:
            bases = read_bases[read_name]
            
            # Handle cases where we have 1 or 2 observations from the same fragment
            if len(bases) == 1:
                total_depth += 1
                if bases[0] == ref_base:
                    ref_count += 1
                elif bases[0] == alt_base:
                    alt_count += 1
            elif bases[0] == bases[1]:
                # Both reads from pair agree
                total_depth += 1
                if bases[0] == ref_base:
                    ref_count += 1
                elif bases[0] == alt_base:
                    alt_count += 1
            else:
                # Reads from pair disagree
                total_depth += 1
                if bases[0] == ref_base or bases[1] == ref_base:
                    ref_count += 1
                elif bases[0] == alt_base and bases[1] == alt_base:
                    alt_count += 1
        
        return ref_count, alt_count, total_depth

    def _format_results(self, ref_count: int, alt_count: int, total_depth: int) -> dict:
        """Format results as strings to prevent pandas float conversion."""
        return {
            'tumor_depth': str(total_depth),
            'tumor_ref_count': str(ref_count),
            'tumor_alt_count': str(alt_count)
        }


def main():
    args = get_args()
    augment_maf = AugmentMAF(
        args.sample_id, args.index_maf, args.index_bam, 
        args.add_maf_files, args.genome_build, args.output, 
        threads=args.threads, min_alt_count=args.alt_count_min,
        compute_umi_metrics=args.compute_umi_metrics,
        min_UMI_3_count=args.min_UMI_3_count
    )
    augment_maf.write_output()

if __name__ == '__main__':
    main()
