**Module Overview**
- **Purpose**: Compare MAF files from two variant callers per tumour/normal pair. Produces caller-unique MAFs and a long-format summary TSV with variant counts and driver gene statistics.

**Required config**

| Key | Description |
|-----|-------------|
| `inputs.maf1` | Path pattern for caller 1 MAF (supports `{seq_type}`, `{genome_build}`, `{tumour_id}`, `{normal_id}`, `{pair_status}`) |
| `inputs.maf2` | Path pattern for caller 2 MAF |
| `options.caller1_name` | Short label for caller 1 (used in output filenames and stats) |
| `options.caller2_name` | Short label for caller 2 |

**Outputs** (under `99-outputs/{seq_type}--{genome_build}/`)

| File | Description |
|------|-------------|
| `{tumour}--{normal}--{pair_status}.{caller1}-only.maf` | Variants called only by caller 1 |
| `{tumour}--{normal}--{pair_status}.{caller2}-only.maf` | Variants called only by caller 2 |
| `{tumour}--{normal}--{pair_status}.stats.tsv` | Long-format summary (3 rows per pair) |

**Stats TSV format**

One row per category (`intersect`, `{caller1}_only`, `{caller2}_only`), with columns: `tumour_id`, `normal_id`, `seq_type`, `genome_build`, `pair_status`, `caller1`, `caller2`, `category`, `count`, `driver_variants`, `driver_genes`.

`driver_variants` is the number of coding variants in bona fide driver genes; `driver_genes` is a comma-separated list of the genes involved. Both are populated for all three categories. Intersect stats are derived from caller 1 rows for shared keys.

**Variant matching**

Two variants are considered the same if `Chromosome`, `Start_Position`, `Reference_Allele`, and `Tumor_Seq_Allele2` all match. This is configurable via `options.key_cols`.

**Driver gene analysis**

By default uses the bundled LLMPP curated tier-1 gene list (`etc/any_tier1_BL_FL_DLBCL_MCL_MZL_PMBL.tsv`). Only coding variant classes (missense, nonsense, frameshift, in-frame indels, splice site, translation start, nonstop) count toward driver gene statistics. Set `options.driver_genes: ""` to disable driver gene reporting entirely.

**Notes**
- Caller-unique MAFs retain all original columns from the respective input MAF.
- Variants present in both callers are not written to an intersect MAF (only counted in the stats).
- The `pair_status` wildcard propagates from the samples table; unmatched tumours require `unmatched_normal_ids` to be set in `_shared`.
