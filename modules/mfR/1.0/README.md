**Module Overview**
- **Purpose**: Cluster non-coding somatic mutation positions into "foci" via hierarchical
  clustering, independently per sample_set and per chromosome, selecting the clustering height
  that maximizes mean silhouette width.

**Required config**

| Key | Description |
|-----|-------------|
| `inputs.sample_maf` | Path pattern for a single sample's genome MAF, already bgzipped + tabix-indexed (a `.tbi` must exist alongside each file). Available wildcards: `{seq_type}`, `{genome_build}`, `{tumour_id}`, `{normal_id}`, `{pair_status}` |
| `inputs.sample_sets` | TSV mapping samples to sample sets. Must contain the columns named by `options.sample_id_column` and `options.sample_set_column` |
| `sample_set` | List of sample_set values to run (must match values in the sample_sets TSV) |

**Why per-sample, tabix-indexed MAFs**

Earlier versions of this module (as `mutation_foci`) took one combined, genome-wide MAF
covering every sample in the project. Every sample_set's job read that entire file (once per
sample_set), and every chromosome's job then re-read the entire per-sample_set MAF just to filter
down to one chromosome (once per chromosome, per sample_set) — for a genome-wide MAF spanning
many samples this reads far more data into memory than any single job needs, and does so
redundantly.

mfR instead takes each sample's MAF individually. For a given sample_set × chromosome job, only
that sample_set's own samples' MAFs are ever opened, and `tabix` is used to pull only that
chromosome's rows directly out of each — no sample outside the set and no chromosome outside the
one being processed is ever read. See `src/python/extract_chrom_maf.py`.

**Outputs** (under `99-outputs/tsv/`)

One TSV per sample_set: `{sample_set}.foci.tsv`, plus (not symlinked to outputs, but kept in
`03-foci/{sample_set}/chromosomes/`) a per-chromosome silhouette-vs-height diagnostic PDF.

**Notes**
- A sample_set's samples are assumed to share one genome_build (and therefore one chromosome
  naming convention — `"1"` vs `"chr1"`). Mixed-build sample_sets are not supported.
- `options.coding_variant_classifications` lists the `Variant_Classification` values dropped
  before clustering (the remainder is treated as non-coding). Set to `[]` if the per-sample MAFs
  are already non-coding-only.
- `hclust_method: "centroid"` (the default) can produce non-monotonic merge heights; switch to
  `"average"` or `"ward.D2"` if foci look wrong for your data.
