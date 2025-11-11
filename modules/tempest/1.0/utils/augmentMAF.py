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
from collections import defaultdict
from typing import Dict, Set, Tuple, List


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_id', required=True, type=str, help='tumour sample id')
    parser.add_argument('--threads', required=False, type=int, default=24, help='number of threads to use, default is 24')
    parser.add_argument('--index_maf', required=True, type=str, help='index maf file')
    parser.add_argument('--index_bam', required=True, type=str, help='index bam file')
    parser.add_argument('--add_maf_files', required=False, nargs='*', type=str, help='additional maf files')
    parser.add_argument('--genome_build', required=True, type=str, help='genome build')
    parser.add_argument('--min_alt_count', required=True, type=int, help='minimum alt count reads needed to keep variant')
    parser.add_argument('--min_t_depth', required=False, type=int, default=50, help='minimum depth needed to keep variant, default is 50')
    parser.add_argument('--min_VAF', required=False, type=float, help='minimum VAF needed to keep variant, default is 0.01')
    parser.add_argument('--output', required=True, type=str, help='output file')

    # Optional: compute UMI metrics (MI and cD tags)
    parser.add_argument('--compute_umi_metrics', action='store_true',
                        help='If set, compute UMI_mean, UMI_max, UMI_3_count for alt-supporting reads (uses MI and cD tags)')
    parser.add_argument("--min_UMI_3_count", type=int, default=1,
                        help="Filters UMI_3_count metric. Which is a count of the number of reads that have UMI family sizes of at least 3.")

    # phased vars arguments
    parser.add_argument("--phase_ID_col", required=False, type=str, help="""Column from input mafs that contains the IDs for the phase sets if a given variant is apart of any,
                                                                            comma delimited if multiple. Empty if none. If not specified the script wont filter by phase sets.""")
    parser.add_argument("--phased_min_t_alt_count", required=False, type=int, default=3, help="Minimum alt count for a variant to be considered for phasing. Default is 3")
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

        # Phase-sets and helpers (computed once)
        self._phase_sets: List[Set[str]] = self._fetch_phase_var_sets()  # list[set[variant_key]]
        self._phase_group_variant_keys: Set[str] = set().union(*self._phase_sets) if self._phase_sets else set()
        self._var_to_set_index: Dict[str, int] = {}
        for i, s in enumerate(self._phase_sets or []):
            for vk in s:
                self._var_to_set_index[vk] = i

        # This will be populated after workers finish (variant_key -> set[qnames])
        self._alt_support_qnames: Dict[str, Set[str]] = defaultdict(set)
        # Fixed: phased evidence requires >=1 fragment
        self._min_phased_fragments_for_phased_label: int = 1

        # Precompute
        self.chromosomes = self.get_genome_chromosomes()
        self.master_maf = self.run_or_not()
        self.augmented_vaf = self.augment_maf()  # core processing

    def read_maf(self, maf_file: str) -> pd.DataFrame:
        return pd.read_csv(
            maf_file,
            sep="\t",
            comment='#',
            dtype={'Chromosome': 'string', 't_ref_count': 'string', 't_alt_count': 'string', 't_depth': 'string',
                   "Tumor_Sample_Barcode": "string"},
            low_memory=False
        )

    def get_genome_chromosomes(self) -> list:
        if self.cfg.genome_build.lower() in ["hg38", "hg19-reddy", "grch38"]:
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

    def augment_maf(self) -> pd.DataFrame:
        if self.master_maf is None:
            print("No additional maf files provided ... returning index maf")
            out_df = self.index_maf_df.copy()
            out_df["origin_samples"] = ""
            # Ensure variant_source is set for consistency
            out_df["variant_source"] = "index_maf"
            return out_df

        new_vars = self.master_maf.loc[self.master_maf["variant_source"] == "additional_maf"].shape[0]
        if new_vars == 0:
            print("No new variants to add to index maf ... returning index maf")
            out_df = self.index_maf_df.copy()
            out_df["origin_samples"] = ""
            out_df["variant_source"] = "index_maf"
            return out_df
        elif new_vars > 0:
            print(f"Found {new_vars} new variants to add to index maf ...")

        # Multiprocessing across chromosomes.
        # Each worker returns (augmented_maf_subset_df, qname_dict_for_that_chrom)
        with multiprocessing.Pool(processes=self.cfg.threads) as pool:
            results = pool.starmap(self._subset_and_run, [[chrm] for chrm in self.chromosomes])

        # Aggregate augmented subsets
        augmented_maf_merged = pd.concat([r[0] for r in results]) if results else pd.DataFrame()

        # Aggregate qname dicts from all workers
        global_qname_dict: Dict[str, Set[str]] = defaultdict(set)
        for _, qd in results:
            for vk, names in qd.items():
                if not names:
                    continue
                if isinstance(names, set):
                    global_qname_dict[vk].update(names)
                else:
                    global_qname_dict[vk].update(set(names))
        # Store globally for phasing computation
        self._alt_support_qnames = global_qname_dict

        # Combine with index MAF, keep all index rows, filter only augmented later
        dfs = [self.index_maf_df.copy()]
        if not augmented_maf_merged.empty:
            # Ensure we don’t introduce all-NA columns that confuse dtype inference
            dfs.append(augmented_maf_merged)
        full_augmented_maf_merged = pd.concat(dfs, ignore_index=True)
        # Label variant_source
        if "variant_source" not in full_augmented_maf_merged.columns:
            full_augmented_maf_merged["variant_source"] = pd.NA
        full_augmented_maf_merged["variant_source"] = full_augmented_maf_merged["variant_source"].fillna("index_maf")
        full_augmented_maf_merged.loc[
            full_augmented_maf_merged.index[-len(augmented_maf_merged):] if len(augmented_maf_merged) else [],
            "variant_source"
        ] = "additional_maf"

        # Standardize sample id
        full_augmented_maf_merged["Tumor_Sample_Barcode"] = self.cfg.sample_id

        # Drop duplicate variants based on variant_key
        full_augmented_maf_merged["variant_key"] = (
            full_augmented_maf_merged["Chromosome"].astype(str) + "_" +
            full_augmented_maf_merged["Start_Position"].astype(str) + "_" +
            full_augmented_maf_merged["Tumor_Seq_Allele2"].astype(str)
        )
        full_augmented_maf_merged = full_augmented_maf_merged.drop_duplicates(subset="variant_key", keep="first").copy()

        # Compute AF for all rows where counts exist
        # Do NOT drop index variants; filtering will only apply to augmented subset.
        numeric_ok = (
            full_augmented_maf_merged["t_alt_count"].notna() &
            full_augmented_maf_merged["t_depth"].notna()
        )
        full_augmented_maf_merged.loc[numeric_ok, "AF"] = (
            full_augmented_maf_merged.loc[numeric_ok, "t_alt_count"].astype(int) /
            full_augmented_maf_merged.loc[numeric_ok, "t_depth"].astype(int)
        )

        # Apply filtering (phasing computed inside)
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
            # only keep rows with IDs in phase col
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

    def _subset_and_run(self, this_chromosome: str) -> Tuple[pd.DataFrame, Dict[str, Set[str]]]:
        """
        Process a single chromosome:
          - Compute counts for augmented (additional_maf) variants on this chromosome, collecting qnames.
          - Collect qnames-only for relevant index variants in phase sets overlapping augmented variants on this chromosome.
          - Return (augmented_subset_with_counts, qname_dict)
        """
        print(f"Working on chromosome {this_chromosome} ...")

        # Subset additional variants for this chromosome
        missing_maf = self.master_maf[self.master_maf["Chromosome"] == this_chromosome].copy()
        missing_maf = missing_maf[missing_maf["variant_source"] == "additional_maf"].copy()

        # Local qname sink for this worker
        qname_sink: Dict[str, Set[str]] = defaultdict(set)

        # 1) Process augmented variants: compute strict counts and collect qnames
        if not missing_maf.empty:
            for var in missing_maf.itertuples():
                variant_key = f"{var.Chromosome}_{int(var.Start_Position)}_{var.Tumor_Seq_Allele2}"
                read_counts = self._get_read_support(
                    chrom=var.Chromosome,
                    start=int(var.Start_Position),
                    end=int(var.End_Position),
                    ref_base=var.Reference_Allele,
                    alt_base=var.Tumor_Seq_Allele2,
                    bamfile=self.cfg.index_bam,
                    mut_class=var.Variant_Type,
                    variant_key=variant_key,
                    alt_qname_sink=qname_sink
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

        # 2) Collect qnames for relevant index variants in the same phase sets as augmented variants (this chromosome)
        relevant_index_rows = pd.DataFrame()
        if not missing_maf.empty and self._phase_sets:
            # Build augmented variant_key set for this chromosome
            missing_maf["variant_key"] = (
                missing_maf["Chromosome"].astype(str) + "_" +
                missing_maf["Start_Position"].astype(int).astype(str) + "_" +
                missing_maf["Tumor_Seq_Allele2"].astype(str)
            )
            augmented_keys = set(missing_maf["variant_key"].tolist())
            augmented_set_idxs = {self._var_to_set_index[vk] for vk in augmented_keys if vk in self._var_to_set_index}

            if augmented_set_idxs:
                # From index_maf_df, select rows on this chromosome that are members of any of these set indices
                idx_df = self.index_maf_df[self.index_maf_df["Chromosome"] == this_chromosome].copy()
                if not idx_df.empty:
                    idx_df["variant_key"] = (
                        idx_df["Chromosome"].astype(str) + "_" +
                        idx_df["Start_Position"].astype(int).astype(str) + "_" +
                        idx_df["Tumor_Seq_Allele2"].astype(str)
                    )
                    idx_df["set_idx"] = idx_df["variant_key"].map(self._var_to_set_index)
                    relevant_index_rows = idx_df[idx_df["set_idx"].isin(augmented_set_idxs)].copy()

        if not relevant_index_rows.empty:
            # For each relevant index variant, collect alt-supporting qnames ONLY (do not recompute counts)
            for row in relevant_index_rows.itertuples():
                variant_key = f"{row.Chromosome}_{int(row.Start_Position)}_{row.Tumor_Seq_Allele2}"
                # Only collect if this variant is part of a phase set
                track_alt = (variant_key in self._phase_group_variant_keys)
                if not track_alt:
                    continue
                # Invoke scanning using the same strict allele logic, but we ignore counts and only capture qnames
                _ = self._get_read_support(
                    chrom=row.Chromosome,
                    start=int(row.Start_Position),
                    end=int(row.End_Position),
                    ref_base=row.Reference_Allele,
                    alt_base=row.Tumor_Seq_Allele2,
                    bamfile=self.cfg.index_bam,
                    mut_class=row.Variant_Type,
                    variant_key=variant_key,
                    alt_qname_sink=qname_sink
                )

        print(f"Finished working on chromosome {this_chromosome} ...")
        return missing_maf, qname_sink

    def _filter_augmented_variants(self, inmaf: pd.DataFrame) -> pd.DataFrame:
        """
        Apply filtering criteria to augmented variants only.
        Index variants pass through unchanged (but may carry phase metadata).
        """
        maf = inmaf.copy()

        # Compute phasing annotations using observed co-occurrence
        maf = self._mark_phased_vars(maf)

        # Split augmented vs index
        is_aug = maf["variant_source"] == "additional_maf"
        maf_aug = maf[is_aug].copy()
        maf_idx = maf[~is_aug].copy()

        # Drop augmented rows with missing counts
        if not maf_aug.empty:
            maf_aug = maf_aug.dropna(subset=["t_alt_count", "t_depth"]).copy()
            # Drop augmented rows with zero alt (do NOT touch index rows)
            maf_aug = maf_aug[maf_aug["t_alt_count"].astype(int) > 0].copy()

            # Recompute AF for augmented safely
            maf_aug["AF"] = (maf_aug["t_alt_count"].astype(int) / maf_aug["t_depth"].astype(int))

        # Apply phased/unphased thresholds to augmented only
        phased_aug = maf_aug[maf_aug["is_phased"]].copy()
        unphased_aug = maf_aug[~maf_aug["is_phased"]].copy()

        if not phased_aug.empty:
            phased_aug = phased_aug[
                phased_aug["t_alt_count"].astype(int) >= self.cfg.phased_min_t_alt_count
            ].copy()

        if not unphased_aug.empty:
            unphased_aug = unphased_aug[
                (unphased_aug["t_alt_count"].astype(int) >= self.cfg.min_alt_count) &
                (unphased_aug["t_depth"].astype(int) > self.cfg.min_t_depth)
            ].copy()

        # UMI filter (augmented only, if metrics exist)
        if self.cfg.compute_umi_metrics and self.cfg.min_UMI_3_count:
            if "UMI_3_count" in phased_aug.columns and not phased_aug.empty:
                phased_aug = phased_aug[phased_aug["UMI_3_count"].astype(int) >= self.cfg.min_UMI_3_count].copy()
            if "UMI_3_count" in unphased_aug.columns and not unphased_aug.empty:
                unphased_aug = unphased_aug[unphased_aug["UMI_3_count"].astype(int) >= self.cfg.min_UMI_3_count].copy()

        # Apply VAF threshold only to augmented rows if requested
        if self.cfg.min_VAF is not None:
            if not phased_aug.empty:
                phased_aug = phased_aug[phased_aug["AF"] >= self.cfg.min_VAF].copy()
            if not unphased_aug.empty:
                unphased_aug = unphased_aug[unphased_aug["AF"] >= self.cfg.min_VAF].copy()

        # Combine back: keep all index rows unfiltered
        outmaf = pd.concat([maf_idx, unphased_aug, phased_aug], ignore_index=True)

        # Clean up helper columns if desired (keep for QC)
        # Ensure Tumor_Sample_Barcode is set
        outmaf["Tumor_Sample_Barcode"] = self.cfg.sample_id

        return outmaf

    def _mark_phased_vars(self, maf: pd.DataFrame) -> pd.DataFrame:
        """
        Adds:
          - phase_set_index, phase_set_size (from pre-established sets)
          - phased_fragment_count (observed co-occurrence fragments via qname)
          - is_phased (phased_fragment_count >= 1)
        """
        maf = maf.copy()
        maf["variant_key"] = (
            maf["Chromosome"].astype(str) + "_" +
            maf["Start_Position"].astype(str) + "_" +
            maf["Tumor_Seq_Allele2"].astype(str)
        )

        # phase_set_index and size
        maf["phase_set_index"] = maf["variant_key"].map(self._var_to_set_index)
        index_counts = maf["phase_set_index"].value_counts(dropna=True).to_dict()
        maf["phase_set_size"] = maf["phase_set_index"].map(index_counts).fillna(0).astype(int)

        # observed phased fragments
        phased_fragment_counts = self._compute_phased_fragment_counts()
        maf["phased_fragment_count"] = maf["variant_key"].map(phased_fragment_counts).fillna(0).astype(int)

        # label is_phased based on observed evidence (>= 1 fragment)
        maf["is_phased"] = maf["phased_fragment_count"] >= self._min_phased_fragments_for_phased_label

        return maf

    def _compute_phased_fragment_counts(self) -> Dict[str, int]:
        """
        Using self._phase_sets and self._alt_support_qnames, compute:
          dict[variant_key] -> phased_fragment_count
        A phased fragment is a query_name that contains >=2 members of the same phase set.
        """
        if not self._phase_sets:
            return {}
        counts = {vk: 0 for vk in self._phase_group_variant_keys}
        for pset in self._phase_sets:
            frags_to_vars = defaultdict(list)  # qname -> [variant_keys in this set]
            for vk in pset:
                for qn in self._alt_support_qnames.get(vk, ()):
                    frags_to_vars[qn].append(vk)
            for qn, members in frags_to_vars.items():
                if len(members) >= 2:
                    for vk in members:
                        counts[vk] += 1
        return counts

    # ---------- BAM processing and allele support ----------

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

    def _get_read_support(self, chrom: str, start: int, end: int, ref_base: str, alt_base: str,
                          bamfile: str, mut_class: str, min_base_qual=10, min_mapping_qual=40,
                          padding: int = 100, variant_key: str = None,
                          alt_qname_sink: Dict[str, Set[str]] = None) -> dict:
        """
        Compute strict ref/alt counts and (optionally) collect alt-supporting qnames into alt_qname_sink.
        """
        params = self._setup_variant_params(start, end, ref_base, mut_class, padding)
        samfile = self._open_bam_file(bamfile)

        if variant_key is None:
            variant_key = f"{chrom}_{start}_{alt_base}"
        track_alt = (variant_key in self._phase_group_variant_keys)

        if mut_class == "INS":
            ref_count, alt_count, total_depth, umi_sizes = self._process_insertion(
                samfile, chrom, params, alt_base,
                min_base_qual, min_mapping_qual, self.cfg.compute_umi_metrics,
                track_alt, variant_key, alt_qname_sink
            )
        elif mut_class == "DEL":
            ref_count, alt_count, total_depth, umi_sizes = self._process_deletion(
                samfile, chrom, params, ref_base,
                min_base_qual, min_mapping_qual, self.cfg.compute_umi_metrics,
                track_alt, variant_key, alt_qname_sink
            )
        elif mut_class in ["SNP", "DNP"]:
            ref_count, alt_count, total_depth, umi_sizes = self._process_snp_dnp(
                samfile, chrom, params, ref_base, alt_base, mut_class,
                min_base_qual, min_mapping_qual, self.cfg.compute_umi_metrics,
                track_alt, variant_key, alt_qname_sink
            )
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

    def _process_snp_dnp(self, samfile, chrom: str, params: dict, ref_base: str, alt_base: str,
                         mut_class: str, min_base_qual: int, min_mapping_qual: int,
                         collect_umi: bool, track_alt: bool, variant_key: str,
                         alt_qname_sink: Dict[str, Set[str]] = None) -> tuple:
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

        def record_umi(aln):
            if collect_umi and aln is not None:
                key, size = self._umi_key_and_size(aln)
                if size is not None and key not in seen_keys:
                    seen_keys.add(key)
                    umi_sizes.append(size)

        for read_name, rec in per_read.items():
            bases = rec["bases"]
            total_depth += 1
            if len(bases) == 1:
                b = bases[0]
                if b == ref_base:
                    ref_count += 1
                elif b == alt_base:
                    alt_count += 1
                    record_umi(rec["alt_aln"])
                    if track_alt and alt_qname_sink is not None:
                        alt_qname_sink[variant_key].add(read_name)
            else:
                if bases[0] == bases[1]:
                    b = bases[0]
                    if b == ref_base:
                        ref_count += 1
                    elif b == alt_base:
                        alt_count += 1
                        record_umi(rec["alt_aln"])
                        if track_alt and alt_qname_sink is not None:
                            alt_qname_sink[variant_key].add(read_name)
                else:
                    if ref_base in bases:
                        ref_count += 1
                    # Mixed non-concordant alt bases are not counted as strict alt support.
        return ref_count, alt_count, total_depth, umi_sizes

    def _process_insertion(self, samfile, chrom: str, params: dict,
                           alt_base: str,
                           min_base_qual: int, min_mapping_qual: int,
                           collect_umi: bool, track_alt: bool, variant_key: str,
                           alt_qname_sink: Dict[str, Set[str]] = None) -> tuple:
        """
        Strict insertion support (REF == '-', ALT = inserted sequence).
        Depth counts reads that:
        - Cover anchor (left base) with base quality >= min_base_qual
        - Cover right flank (reference base immediately after anchor)
        ALT-supporting read must additionally:
        - Contain a single CIGAR I op of length == len(ALT) occurring immediately after anchor
        - Inserted sequence exactly matches ALT
        Right flank bridging ensures the read spans across the insertion context.
        """
        start_1 = params['original_start']          # 1-based position where insertion occurs AFTER this base
        anchor_pos = start_1 - 1                    # 0-based anchor base
        inserted_seq = alt_base or ""
        ins_len = len(inserted_seq)
        if ins_len == 0:
            return 0, 0, 0, []

        # Right flank = next reference position after anchor
        right_flank_pos = anchor_pos + 1

        # Fetch small padded region
        fetch_start = max(anchor_pos - 10, 0)
        fetch_end = right_flank_pos + 10

        total_depth = 0
        alt_count = 0
        umi_sizes = []
        seen_umi = set()

        for read in samfile.fetch(chrom, fetch_start, fetch_end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapping_qual:
                continue

            cig = read.cigartuples or []
            # Pre-screen: any insertion of required length?
            has_req_len_ins = any(op == 1 and length == ins_len for op, length in cig)

            ref_positions = read.get_reference_positions()
            # Must bridge anchor and right flank to be counted in depth
            if anchor_pos not in ref_positions or right_flank_pos not in ref_positions:
                continue

            # Anchor base quality
            anchor_qpos = None
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos == anchor_pos:
                    anchor_qpos = qpos
                    break
            if anchor_qpos is None:
                continue
            if read.query_qualities and read.query_qualities[anchor_qpos] < min_base_qual:
                continue

            # Read contributes to depth now
            total_depth += 1

            if not has_req_len_ins:
                # No insertion of correct length => ref-supporting read
                continue

            # Walk CIGAR to confirm exact insertion position & sequence
            ref_cursor = read.reference_start
            read_cursor = 0
            alt_supported = False

            for op, length in cig:
                # 0=M,7==,8=X consume ref & read
                if op in (0, 7, 8):
                    ref_cursor += length
                    read_cursor += length
                elif op == 1:  # insertion consumes read only
                    # Insertion occurs between ref_cursor - 1 and ref_cursor
                    if ref_cursor - 1 == anchor_pos and length == ins_len:
                        ins_seq_read = read.query_sequence[read_cursor: read_cursor + length]
                        if ins_seq_read == inserted_seq:
                            alt_supported = True
                            break
                    read_cursor += length
                elif op == 2:  # deletion consumes ref
                    ref_cursor += length
                elif op == 3:  # N consumes ref
                    ref_cursor += length
                elif op == 4:  # soft clip consumes read
                    read_cursor += length
                elif op == 5:  # hard clip consumes neither
                    pass
                else:
                    # Other ops (P, etc.) ignore
                    pass

            if alt_supported:
                alt_count += 1
                if track_alt and alt_qname_sink is not None:
                    alt_qname_sink[variant_key].add(read.query_name)
                if collect_umi:
                    key, size = self._umi_key_and_size(read)
                    if size is not None and key not in seen_umi:
                        seen_umi.add(key)
                        umi_sizes.append(size)

        ref_count = total_depth - alt_count
        return ref_count, alt_count, total_depth, umi_sizes

    def _process_deletion(self, samfile, chrom: str, params: dict,
                          ref_base: str,
                          min_base_qual: int, min_mapping_qual: int,
                          collect_umi: bool, track_alt: bool, variant_key: str,
                          alt_qname_sink: Dict[str, Set[str]] = None) -> tuple:
        """
        Strict deletion support (REF = exactly deleted bases, ALT='-').

        A read contributes to depth if:
        - It covers the left anchor (base before the deletion) with base quality >= min_base_qual, and
        - It covers the right flank (base immediately after the deleted span).
        A read is ALT if, in addition:
        - It contains a CIGAR 'D' that starts exactly at deletion_start_0 (0-based)
            and has length == len(REF).
        """
        # Variant geometry
        start_1 = params['original_start']          # 1-based first deleted base
        deletion_start_0 = start_1 - 1              # 0-based
        del_len = len(ref_base)                     # strict allele length
        if del_len <= 0:
            return 0, 0, 0, []

        deletion_end_0 = deletion_start_0 + del_len - 1
        anchor_pos = deletion_start_0 - 1
        right_flank_pos = deletion_end_0 + 1

        # Small padding for fetch robustness
        fetch_start = max(anchor_pos - 10, 0)
        fetch_end = right_flank_pos + 10

        total_depth = 0
        alt_count = 0
        umi_sizes = []
        seen_umi = set()

        for read in samfile.fetch(chrom, fetch_start, fetch_end):
            # Primary, mapped, MAPQ gate
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapping_qual:
                continue

            cig = read.cigartuples or []

            # Fast pre-screen: any deletion of the required length?
            if not any(op == 2 and length == del_len for op, length in cig):
                # Still may contribute to depth if it bridges
                ref_positions = read.get_reference_positions()
                if anchor_pos in ref_positions and right_flank_pos in ref_positions:
                    # Anchor base quality
                    anchor_qpos = None
                    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                        if rpos == anchor_pos:
                            anchor_qpos = qpos
                            break
                    if anchor_qpos is not None and (not read.query_qualities or read.query_qualities[anchor_qpos] >= min_base_qual):
                        total_depth += 1
                continue

            # Must bridge anchor and right flank
            ref_positions = read.get_reference_positions()
            if anchor_pos not in ref_positions or right_flank_pos not in ref_positions:
                continue

            # Anchor base quality
            anchor_qpos = None
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos == anchor_pos:
                    anchor_qpos = qpos
                    break
            if anchor_qpos is None:
                continue
            if read.query_qualities and read.query_qualities[anchor_qpos] < min_base_qual:
                continue

            total_depth += 1

            # Confirm exact deletion start and length via CIGAR walk
            ref_cursor = read.reference_start
            exact_d_found = False
            for op, length in cig:
                if op in (0, 7, 8):  # M, =, X consume ref
                    ref_cursor += length
                elif op == 2:        # D consumes ref
                    if ref_cursor == deletion_start_0 and length == del_len:
                        exact_d_found = True
                        break
                    ref_cursor += length
                elif op == 3:        # N consumes ref
                    ref_cursor += length
                # I/S/H do not advance ref
                else:
                    pass

            if not exact_d_found:
                # counts as REF
                continue

            # ALT-supporting read
            alt_count += 1
            if track_alt and alt_qname_sink is not None:
                alt_qname_sink[variant_key].add(read.query_name)
            if collect_umi:
                key, size = self._umi_key_and_size(read)
                if size is not None and key not in seen_umi:
                    seen_umi.add(key)
                    umi_sizes.append(size)

        ref_count = total_depth - alt_count
        return ref_count, alt_count, total_depth, umi_sizes

    # ---------- UMI helpers ----------

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

    # ---------- SNP/DNP helpers ----------

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

    # ---------- Formatting ----------

    def _format_results(self, ref_count: int, alt_count: int, total_depth: int) -> dict:
        return {
            'tumor_depth': str(total_depth),
            'tumor_ref_count': str(ref_count),
            'tumor_alt_count': str(alt_count)
        }

    # ---------- Additional maf helpers ----------

    def _get_missing_maf_rows(self):
        additional_mafs = []
        for maf in (self.cfg.add_maf_files or []):
            other_maf = self.read_maf(maf)
            other_maf["variant_source"] = "additional_maf"
            additional_mafs.append(other_maf)

        additional_vars = self.mark_origin(pd.concat(additional_mafs))
        all_mafs = pd.concat([self.index_maf_df, additional_vars]).reset_index(drop=True).copy()
        all_mafs["variant_key"] = (
            all_mafs["Chromosome"].astype(str) + "_" +
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
            add_mafs["Chromosome"].astype(str) + "_" +
            add_mafs["Start_Position"].astype(str) + "_" +
            add_mafs["Tumor_Seq_Allele2"].astype(str)
        )
        variant_key_dict = add_mafs.groupby("variant_key")["Tumor_Sample_Barcode"].agg(
            lambda x: ",".join(list(x))
        ).to_dict()
        add_mafs["origin_samples"] = add_mafs["variant_key"].map(variant_key_dict)
        return add_mafs


def main():
    args = get_args()
    cfg = AugmnentMAFArgs(**vars(args))
    augment_maf_runner = AugmentMAF(cfg)
    augment_maf_runner.write_output()


if __name__ == '__main__':
    main()
