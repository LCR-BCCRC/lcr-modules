"""
Extract one chromosome's rows from one sample_set's per-sample MAFs.

Run under Snakemake's `script:` directive for the `_mfR_extract_chrom` rule.

For each sample belonging to this sample_set, shells out to `tabix` and
streams its stdout one line at a time -- so a job never reads another
chromosome's rows, never reads another sample's MAF at all, and never holds
more than one line (or one sample's worth of file-handle buffering) in memory
at once. Rows in the configured coding Variant_Classification set are dropped
inline as each line streams past; everything else is written straight
through to the output file, unparsed. This replaces the old prepare_maf.R
(which read the whole genome-wide master MAF per sample_set) and the old
per-chromosome read+filter in cluster_foci.R (which read the whole
per-sample_set MAF per chromosome).

No pandas/polars here on purpose: once tabix has scoped the data down to one
sample_set x one chromosome, the only work left is a single-column row
filter, which a line-by-line stream handles without ever materializing the
chromosome's rows as a table in memory.
"""

import gzip
import subprocess
import sys

log = open(snakemake.log[0], "w")
sys.stderr = log

maf_gz_paths = list(snakemake.input.mafs)
chrom = snakemake.wildcards.chrom
sample_set = snakemake.wildcards.sample_set
coding_classes = set(snakemake.params.coding_classes)

print(f"sample_set '{sample_set}', chromosome '{chrom}': {len(maf_gz_paths)} sample MAF(s)", file=log)

# Header is read directly off the first sample's MAF (not via tabix, whose
# region-query output never includes header/comment lines) so the
# Variant_Classification column position is known before streaming rows.
with gzip.open(maf_gz_paths[0], "rt") as fh:
    header_line = fh.readline().rstrip("\n")
columns = header_line.split("\t")
vc_idx = columns.index("Variant_Classification") if "Variant_Classification" in columns else None
if vc_idx is None and coding_classes:
    print("WARNING: Variant_Classification column not found; coding-class filter skipped", file=log)

total_rows = 0
kept_rows = 0
dropped_rows = 0

with open(snakemake.output.maf, "w") as out:
    out.write(header_line + "\n")
    for maf_gz in maf_gz_paths:
        sample_rows = 0
        sample_dropped = 0
        with subprocess.Popen(["tabix", maf_gz, chrom], stdout=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                sample_rows += 1
                if vc_idx is not None and coding_classes:
                    fields = line.rstrip("\n").split("\t")
                    if fields[vc_idx] in coding_classes:
                        sample_dropped += 1
                        continue
                out.write(line if line.endswith("\n") else line + "\n")
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, ["tabix", maf_gz, chrom])
        print(f"  {maf_gz}: {sample_rows} rows on {chrom}, {sample_dropped} dropped (coding)", file=log)
        total_rows += sample_rows
        dropped_rows += sample_dropped
        kept_rows += sample_rows - sample_dropped

print(
    f"Total: {total_rows} rows on {chrom} across {len(maf_gz_paths)} samples; "
    f"{dropped_rows} dropped (coding); {kept_rows} written",
    file=log,
)
print(f"Wrote {snakemake.output.maf}", file=log)

log.close()
