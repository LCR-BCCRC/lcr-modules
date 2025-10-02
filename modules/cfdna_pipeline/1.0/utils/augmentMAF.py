import pandas as pd
from pathlib import Path 
import pysam
import os
import multiprocessing
import argparse
import shutil
from dataclasses import dataclass
import math
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_id',required=True,type=str,help='tumour sample id')
    parser.add_argument('--threads',required=False,type=int,default=24,help='number of threads to use, default is 24')
    parser.add_argument('--index_maf',required=True,type=str,help='index maf file')
    parser.add_argument('--index_bam',required=True,type=str,help='index bam file')
    parser.add_argument('--add_maf_files',required=False,nargs='*',type=str,help='additional maf files')
    parser.add_argument('--genome_build',required=True,type=str,help='genome build')
    parser.add_argument('--min_alt_count',required=True,type=int,help='minimum alt count reads needed to keep variant')
    parser.add_argument('--min_t_depth',required=False,type=int,default=50,help='minimum depth needed to keep variant, default is 50')
    parser.add_argument('--min_VAF',required=False,type=float, help='minimum VAF needed to keep variant, default is 0.01')
    parser.add_argument('--output',required=True,type=str,help='output file')

    # Optional: compute UMI metrics (MI and cD tags)
    parser.add_argument('--compute_umi_metrics', action='store_true',
                        help='If set, compute UMI_mean, UMI_max, UMI_3_count for alt-supporting reads (uses MI and cD tags)')
    parser.add_argument("--min_UMI_3_count", type=int, default=1,
                        help="Filters UMI_3_count metric. Which is a count of the number of reads that have UMI family sizes of at least 3.")

    # phased vars arguments 
    parser.add_argument("--phase_ID_col",required=False, type=str, help="""Column from input mafs that contains the IDs for the phase sets if a given variant is apart of any,
                                                                            comma delimited if multiple. Empty if none. If not specified the script wont filter by phase sets.""")
    parser.add_argument("--phased_min_t_alt_count",required=False, type=int, default=3, help="Minimum alt count for a variant to be considered for phasing. Default is 3")
    return parser.parse_args()

@dataclass(frozen=True)
class AugmnentMAFArgs:
    sample_id: str
    threads: int
    index_maf: str
    index_bam: str
    add_maf_files: list
    genome_build: str
    output: str
    min_alt_count: int
    min_t_depth: int = 50
    compute_umi_metrics: bool = False
    min_UMI_3_count: int = None
    phase_ID_col: str = None
    phased_min_t_alt_count: int = 3
    min_VAF: float = None

class AugmentMAF(object):
    """Add variants complete with read depths found in additional mafs to an index maf."""

    def __init__(self, config: AugmnentMAFArgs) -> None:
        super(AugmentMAF, self).__init__()
        self.cfg = config

        # Load index maf once as DataFrame
        self.index_maf_df = self.read_maf(self.cfg.index_maf)

        # Precompute
        self.chromosomes = self.get_genome_chromosomes()
        self.master_maf = self.run_or_not()
        self.augmented_vaf = self.augment_maf()  # core processing

    def read_maf(self, maf_file: str) -> pd.DataFrame:
        return pd.read_csv(
            maf_file,
            sep="\t",
            comment='#',
            dtype={'Chromosome':'string','t_ref_count':'string','t_alt_count':'string','t_depth':'string', "Tumor_Sample_Barcode":"string"},
            low_memory=False
        )

    def get_genome_chromosomes(self)-> list:
        if self.cfg.genome_build.lower() in ["hg38","hg19-reddy", "grch38"]:
            return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        else:
            return [str(i) for i in range(1, 23)] + ["X", "Y"]

    def write_output(self, output_file: str = None) -> None:
        if output_file:
            out_name = output_file
        elif self.cfg.output:
            out_name = self.cfg.output
        else:
            raise RuntimeError("No output file name provided.")
        self.augmented_vaf.to_csv(out_name, sep="\t", index=False)

    def augment_maf(self):
        if self.master_maf is None:
            print("No additional maf files provided ... returning index maf")
            self.index_maf_df["origin_samples"] = ""
            return self.index_maf_df

        new_vars = self.master_maf.loc[self.master_maf["variant_source"] == "additional_maf"].shape[0]
        if new_vars == 0:
            print("No new variants to add to index maf ... returning index maf")
            self.index_maf_df["origin_samples"] = ""
            return self.index_maf_df
        elif new_vars > 0:
            print(f"Found {new_vars} new variants to add to index maf ...")

        pool = multiprocessing.Pool(processes=self.cfg.threads)
        results = pool.starmap(self._subset_and_run, [[chrm] for chrm in self.chromosomes])
        pool.close()
        pool.join()

        augmented_maf_merged = pd.concat(results)
        full_augmented_maf_merged = pd.concat([self.index_maf_df, augmented_maf_merged])
        full_augmented_maf_merged["Tumor_Sample_Barcode"] = self.cfg.sample_id
        full_augmented_maf_merged = full_augmented_maf_merged.drop_duplicates().reset_index(drop=True)

        if "key" in full_augmented_maf_merged.columns:
            full_augmented_maf_merged.drop(columns=["key"], inplace=True)

        full_augmented_maf_merged = full_augmented_maf_merged.dropna(subset=["Hugo_Symbol","t_alt_count","t_depth"]).copy()
        full_augmented_maf_merged["AF"] = (
            full_augmented_maf_merged["t_alt_count"].astype(int) /
            full_augmented_maf_merged["t_depth"].astype(int)
        )
        full_augmented_maf_merged.insert(43, "AF", full_augmented_maf_merged.pop("AF"))

        return self._filter_augmented_variants(full_augmented_maf_merged)

    def _merge_overlapping_sets(self, sets):
        """
        Merge list of input sets. Sets with at least one
        element are combined.
        """
        work = [s.copy() for s in sets if s]  # drop empties early
        i = 0
        while i < len(work):
            j = i + 1
            while j < len(work):
                if work[i] & work[j]:  # they overlap
                    work[i] |= work[j]  # in-place union
                    work.pop(j)         # remove merged set
                    # do not advance j; new work[i] might intersect the next one now at index j
                else:
                    j += 1
            i += 1
        return work

    def _fetch_phase_var_sets(self) -> list:
        """Take all input mafs and returns a final list with sets of phased variant IDs found.

        Returns:
            List of sets, where each set contains variant keys (chr_pos_alt) that are phased together.
        """

        phase_col = getattr(self.cfg, "phase_ID_col", None)
        add_mafs = getattr(self.cfg, "add_maf_files", None)
        index_maf = getattr(self.cfg, "index_maf", None)

        print("\nAttempting to fetch phased variant sets ...")
        if not phase_col or not add_mafs or not index_maf:
            return []

        phased_vars = []

        for inmaf_path in ([index_maf] + add_mafs):
            # read in maf
            inmaf = self.read_maf(inmaf_path)
            if inmaf.empty:
                continue
            # drop rows with empty phase ID col
            if inmaf[phase_col].isnull().all():
                continue
            # only keep rows with IDs  in phase col
            inmaf = inmaf[~inmaf[phase_col].isnull()].copy()

            # create variant key
            inmaf["variant_key"] = inmaf["Chromosome"].astype(str) + "_" + inmaf["Start_Position"].astype(str) + "_" + inmaf["Tumor_Seq_Allele2"].astype(str)
            # grab phase IDs
            all_phase_ids = inmaf[phase_col].fillna("").astype(str).str.split(",").explode().str.strip()
            all_phase_ids = all_phase_ids[all_phase_ids != ""].unique().tolist()
            # skip if none found
            if not all_phase_ids:
                continue
                
            # populate pid dict
            local_sets = {pid: set() for pid in all_phase_ids}
            # fill sets
            for row in inmaf.itertuples(index=False):
                row_pids = getattr(row, phase_col, "")
                for pid in all_phase_ids:
                    if any(p.strip() == pid for p in str(row_pids).split(",")):
                        local_sets[pid].add(row.variant_key)

            phased_vars.extend([s for s in local_sets.values() if len(s) > 1])
        # merge overlapping sets
        final_phase_sets = self._merge_overlapping_sets(phased_vars)
        print(f"Found {len(final_phase_sets)} phased variant sets.")
        print(final_phase_sets)

        return final_phase_sets

    def _mark_phased_vars(self, maf: pd.DataFrame) -> pd.DataFrame:
        """take input maf, and mark variants that are present with more
        than one of the variants in the same phase set.
        """
        phased_sets = self._fetch_phase_var_sets()
        if not phased_sets:
            # has empty phase cols, with NaN values
            maf["phase_set_index"] = pd.NA
            maf["phase_set_size"] = 0
            return maf

        # create variant key
        maf["variant_key"] = (
            maf["Chromosome"].astype(str) + "_" +
            maf["Start_Position"].astype(str) + "_" +
            maf["Tumor_Seq_Allele2"].astype(str)
        )

        # create map of variant_key to phase set index
        var_to_set = {}
        for i, pset in enumerate(phased_sets):
            for var in pset:
                var_to_set[var] = i

        maf["phase_set_index"] = maf["variant_key"].map(var_to_set)
        index_counts = maf["phase_set_index"].value_counts().to_dict()
        maf["phase_set_size"] = maf["phase_set_index"].map(index_counts).fillna(0).astype(int)
        # mark as phased if part of a set with more than one variant
        maf["is_phased"] = maf["phase_set_size"] > 1

        return maf

    def _filter_augmented_variants(self, inmaf: pd.DataFrame) -> pd.DataFrame:
        """Apply filtering criteria to augmented maf.
        
        Filters based on min read support and UMI family support.
        If any variant in a phased set (phase_set_size > 1) meets the relaxed
        phased threshold, all variants in that set are kept (regardless of
        their individual t_alt_count). Non-phased variants (and phased sets
        with no member passing) are filtered by the standard thresholds.
        """
        maf = inmaf.copy()
        maf = self._mark_phased_vars(maf)

        # Split phased vs non-phased
        phased_vars = maf[maf["phase_set_size"] > 1].copy()
        maf_nonphased = maf[maf["phase_set_size"] <= 1].copy()

        # Keep whole phase sets if ANY member meets relaxed threshold
        if not phased_vars.empty:
            keeper_phase_indices = phased_vars.loc[
                phased_vars["t_alt_count"].astype(int) >= self.cfg.phased_min_t_alt_count,
                "phase_set_index"
            ].unique()
            phased_vars = phased_vars[phased_vars["phase_set_index"].isin(keeper_phase_indices)]

        # Standard filtering for non-phased (and phased sets that didn't trigger)
        maf_nonphased = maf_nonphased[
            (maf_nonphased["t_alt_count"].astype(int) >= self.cfg.min_alt_count) &
            (maf_nonphased["t_depth"].astype(int) > self.cfg.min_t_depth)
        ].copy()

        if self.cfg.min_UMI_3_count and maf_nonphased.shape[0] > 0:
            maf_nonphased = maf_nonphased[
                maf_nonphased["UMI_3_count"].astype(int) >= self.cfg.min_UMI_3_count
            ].copy()

        outmaf = pd.concat([maf_nonphased, phased_vars], ignore_index=True)

        # if VAF threshold set, apply it now
        if self.cfg.min_VAF is not None and outmaf.shape[0] > 0:
            outmaf = outmaf[outmaf["AF"] >= self.cfg.min_VAF].copy()

        return outmaf


    def _get_missing_maf_rows(self):
        additional_mafs = []
        for maf in (self.cfg.add_maf_files or []):
            other_maf = self.read_maf(maf)
            other_maf["variant_source"] = "additional_maf"
            additional_mafs.append(other_maf)

        additional_vars = self.mark_origin(pd.concat(additional_mafs))
        all_mafs = pd.concat([self.index_maf_df, additional_vars]).reset_index(drop=True).copy()
        all_mafs["variant_key"] = (
            all_mafs["Chromosome"] + "_" +
            all_mafs["Start_Position"].astype(str) + "_" +
            all_mafs["Tumor_Seq_Allele2"].astype(str)
        )
        all_mafs = all_mafs.drop_duplicates(subset="variant_key", keep="first").copy()
        all_mafs.drop(columns=["variant_key"], inplace=True)
        return all_mafs

    def run_or_not(self):
        if self.cfg.add_maf_files and len(self.cfg.add_maf_files) > 0:
            return self._get_missing_maf_rows()
        else:
            return None

    def mark_origin(self, add_mafs: pd.DataFrame) -> pd.DataFrame:
        add_mafs["variant_key"] = (
            add_mafs["Chromosome"] + "_" +
            add_mafs["Start_Position"].astype(str) + "_" +
            add_mafs["Tumor_Seq_Allele2"].astype(str)
        )
        variant_key_dict = add_mafs.groupby("variant_key")["Tumor_Sample_Barcode"].agg(
            lambda x: ",".join(list(x))
        ).to_dict()
        add_mafs["origin_samples"] = add_mafs["variant_key"].map(variant_key_dict)
        return add_mafs

    def _subset_and_run(self, this_chromosome: str) -> pd.DataFrame:
        print(f"Working on chromosome {this_chromosome} ...")
        missing_maf = self.master_maf[self.master_maf["Chromosome"] == this_chromosome].copy()
        missing_maf = missing_maf[missing_maf["variant_source"] == "additional_maf"].copy()
        print(f"Found these variants to augment: {missing_maf}")

        for var in missing_maf.itertuples():
            read_counts = self._get_read_support(
                chrom=var.Chromosome,
                start=int(var.Start_Position),
                end=int(var.End_Position),
                ref_base=var.Reference_Allele,
                alt_base=var.Tumor_Seq_Allele2,
                bamfile=self.cfg.index_bam,
                mut_class=var.Variant_Type
            )

            missing_maf.at[var.Index, "t_depth"] = read_counts["tumor_depth"]
            missing_maf.at[var.Index, "t_ref_count"] = read_counts["tumor_ref_count"]
            missing_maf.at[var.Index, "t_alt_count"] = read_counts["tumor_alt_count"]

            if self.cfg.compute_umi_metrics:
                umi_mean = read_counts.get("UMI_mean")
                umi_max = read_counts.get("UMI_max")
                umi_3 = read_counts.get("UMI_3_count")
                missing_maf.at[var.Index, "UMI_mean"] = float(umi_mean) if umi_mean is not None else math.nan
                missing_maf.at[var.Index, "UMI_max"] = int(umi_max) if umi_max is not None else math.nan
                missing_maf.at[var.Index, "UMI_3_count"] = int(umi_3) if umi_3 is not None else math.nan

        print(f"Finished working on chromosome {this_chromosome} ...")
        return missing_maf

    def _get_read_support(self, chrom: str, start: int, end: int, ref_base: str, alt_base: str, 
                          bamfile: str, mut_class: str, min_base_qual=10, min_mapping_qual=1,
                          padding: int = 100) -> dict:
        params = self._setup_variant_params(start, end, ref_base, mut_class, padding)
        samfile = self._open_bam_file(bamfile)

        if mut_class == "INS":
            ref_count, alt_count, total_depth, umi_sizes = self._process_insertion(
                samfile, chrom, params, min_base_qual, min_mapping_qual, self.cfg.compute_umi_metrics)
        elif mut_class == "DEL":
            ref_count, alt_count, total_depth, umi_sizes = self._process_deletion(
                samfile, chrom, params, ref_base, min_base_qual, min_mapping_qual, self.cfg.compute_umi_metrics)
        elif mut_class in ["SNP", "DNP"]:
            ref_count, alt_count, total_depth, umi_sizes = self._process_snp_dnp(
                samfile, chrom, params, ref_base, alt_base, mut_class, 
                min_base_qual, min_mapping_qual, self.cfg.compute_umi_metrics)
        else:
            samfile.close()
            raise ValueError(f"Unsupported mutation class: {mut_class}")

        samfile.close()
        result = self._format_results(ref_count, alt_count, total_depth)

        if self.cfg.compute_umi_metrics:
            umi_mean, umi_max, umi_ge3 = self._compute_umi_metrics(umi_sizes)
            result["UMI_mean"] = umi_mean
            result["UMI_max"] = umi_max
            result["UMI_3_count"] = umi_ge3

        return result

    # (All helper methods below unchanged except for using self.cfg where needed) 
    def _setup_variant_params(self, start: int, end: int, ref_base: str, mut_class: str, padding: int) -> dict:
        actual_start = start - 1
        actual_end = end
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
        realname = str(Path(bamfile).resolve())
        if realname.endswith("cram"):
            return pysam.AlignmentFile(realname, 'rc')
        else:
            return pysam.AlignmentFile(bamfile, 'rb')

    def _process_insertion(self, samfile, chrom: str, params: dict, min_base_qual: int, min_mapping_qual: int,
                           collect_umi: bool) -> tuple:
        ref_count = alt_count = total_depth = 0
        umi_sizes = []
        seen = set()
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
                if pu.indel and pu.indel > 0:
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
                umi_at_max_alt = col_umis
        ref_count = max_depth - max_alt
        alt_count = max_alt
        total_depth = max_depth
        return ref_count, alt_count, total_depth, umi_at_max_alt

    def _process_snp_dnp(self, samfile, chrom: str, params: dict, ref_base: str, alt_base: str, 
                         mut_class: str, min_base_qual: int, min_mapping_qual: int,
                         collect_umi: bool) -> tuple:
        per_read = {}
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
                    if ref_base in bases:
                        ref_count += 1
        return ref_count, alt_count, total_depth, umi_sizes

    def _umi_key_and_size(self, aln):
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
        ref_count = alt_count = total_depth = 0
        for read_name in read_bases:
            bases = read_bases[read_name]
            if len(bases) == 1:
                total_depth += 1
                if bases[0] == ref_base:
                    ref_count += 1
                elif bases[0] == alt_base:
                    alt_count += 1
            elif bases[0] == bases[1]:
                total_depth += 1
                if bases[0] == ref_base:
                    ref_count += 1
                elif bases[0] == alt_base:
                    alt_count += 1
            else:
                total_depth += 1
                if bases[0] == ref_base or bases[1] == ref_base:
                    ref_count += 1
                elif bases[0] == alt_base and bases[1] == alt_base:
                    alt_count += 1
        return ref_count, alt_count, total_depth

    def _format_results(self, ref_count: int, alt_count: int, total_depth: int) -> dict:
        return {
            'tumor_depth': str(total_depth),
            'tumor_ref_count': str(ref_count),
            'tumor_alt_count': str(alt_count)
        }

def main():
    args = get_args()
    cfg = AugmnentMAFArgs(**vars(args))
    augment_maf_runner = AugmentMAF(cfg)
    augment_maf_runner.write_output()

if __name__ == '__main__':
    main()
