"""
Extract one chromosome's rows from one sample_set's per-sample MAFs.

Run under Snakemake's `script:` directive for the `_mfr_extract_chrom` rule.

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

Logging notes: the log handle is line-buffered and duplicated onto fd 2, so
(a) diagnostics survive a crash instead of dying in an unflushed buffer, and
(b) tabix's own stderr -- stale-index and unknown-region complaints, which are
the usual cause of failure here -- lands in the same file rather than only in
the scheduler's job stderr.

Binary-safety notes: tabix indexes are BGZF, which is gzip-compatible, so
`gzip.open` on a .tbi succeeds and yields decompressed binary as str rather
than raising. The header read below is therefore bounded (a fixed-size read,
not `for line in fh`, since binary content may contain no newline for
megabytes) and validated before any of it reaches a log or an exception
message. Every interpolation of file content goes through `_preview()` so a
misidentified file produces a short, printable diagnostic instead of dumping
an index across the terminal.
"""

import gzip
import io
import os
import shutil
import subprocess
import sys
import traceback

# Line-buffered: partial output survives an uncaught exception, so the log
# shows how far the job got rather than coming back empty.
log = open(snakemake.log[0], "w", buffering=1)

# Point the OS-level stderr fd at the log too. Rebinding sys.stderr alone is
# not enough -- subprocesses inherit fd 2, not Python's sys.stderr object, so
# tabix's messages would otherwise bypass the log entirely.
sys.stderr = log

# How much of the first non-comment line to read when sniffing the header.
# Comfortably larger than any real MAF header row, small enough that reading
# it from a binary file costs nothing.
_HEADER_SNIFF_BYTES = 65536


def _preview(s, n=80):
    """Printable, length-capped preview of file content for messages.

    Non-printable characters become '.', so a misidentified binary file
    yields a short readable marker rather than flooding the terminal or log
    with escaped index bytes.
    """
    trimmed = s[:n]
    cleaned = "".join(c if c.isprintable() or c == "\t" else "." for c in trimmed)
    suffix = "..." if len(s) > n else ""
    return cleaned + suffix


def _read_header(path):
    """Return (header_line, n_comments) from a bgzipped MAF.

    Read via a bounded chunk rather than line iteration: if `path` is
    accidentally a .tbi (or any other binary), there may be no newline for a
    very long stretch and `for line in fh` would pull all of it into memory
    before anything could inspect it.
    """
    with gzip.open(path, "rt", errors="replace") as fh:
        chunk = fh.read(_HEADER_SNIFF_BYTES)

    if not chunk:
        raise ValueError(f"{path} is empty after decompression.")

    n_comments = 0
    header_line = None
    for line in chunk.split("\n"):
        if line.startswith("#"):
            n_comments += 1
            continue
        header_line = line.rstrip("\r")
        break

    if header_line is None:
        raise ValueError(
            f"No header line found in the first {_HEADER_SNIFF_BYTES} bytes "
            f"of {path} (only comment lines?)."
        )

    # Validate before this string reaches any message or the output file.
    if "\x00" in header_line or not header_line.strip():
        raise ValueError(
            f"{path} does not look like a MAF: first non-comment line "
            f"contains NUL bytes or is blank. Is this a .tbi index rather "
            f"than the .maf.gz? Content begins: {_preview(header_line)!r}"
        )
    if "\t" not in header_line:
        raise ValueError(
            f"{path} does not look like a tab-delimited MAF: first "
            f"non-comment line has no tabs. Content begins: "
            f"{_preview(header_line)!r}"
        )

    return header_line, n_comments


def main():
    maf_gz_paths = list(snakemake.input.mafs)
    chrom = snakemake.wildcards.chrom
    sample_set = snakemake.wildcards.sample_set
    coding_classes = set(snakemake.params.coding_classes)

    print(
        f"sample_set '{sample_set}', chromosome '{chrom}': "
        f"{len(maf_gz_paths)} sample MAF(s)",
        file=log,
    )

    # Echo the resolved inputs: if the maf and tbi lists ever get crossed or
    # concatenated upstream, this is where it becomes obvious.
    for p in maf_gz_paths:
        print(f"  input maf: {p}", file=log)

    if not maf_gz_paths:
        raise ValueError(
            f"No input MAFs for sample_set '{sample_set}'. Check the "
            f"sample_set -> sample mapping in the sample_sets TSV."
        )

    stray_tbis = [p for p in maf_gz_paths if p.endswith(".tbi")]
    if stray_tbis:
        raise ValueError(
            f"input.mafs contains {len(stray_tbis)} .tbi path(s), e.g. "
            f"{stray_tbis[0]}. The index files belong in input.tbis only."
        )

    if shutil.which("tabix") is None:
        raise RuntimeError(
            "tabix not found on PATH. Ensure htslib is present in the "
            "module's conda env / container."
        )

    # Header is read directly off the first sample's MAF (not via tabix, whose
    # region-query output never includes header/comment lines) so the
    # Variant_Classification column position is known before streaming rows.
    # Leading '#' comment lines are skipped -- MAFs commonly start with e.g.
    # '#version 2.4', and treating that as the header would silently disable
    # the coding filter and emit a bogus header downstream.
    header_line, n_comments = _read_header(maf_gz_paths[0])
    columns = header_line.split("\t")
    print(
        f"header from {maf_gz_paths[0]}: {len(columns)} columns "
        f"({n_comments} comment line(s) skipped)",
        file=log,
    )

    vc_idx = columns.index("Variant_Classification") if "Variant_Classification" in columns else None
    if vc_idx is None and coding_classes:
        print(
            "WARNING: Variant_Classification column not found; "
            "coding-class filter skipped",
            file=log,
        )
    elif vc_idx is not None:
        print(
            f"Variant_Classification at column {vc_idx + 1}; "
            f"dropping {len(coding_classes)} coding class(es)",
            file=log,
        )

    total_rows = 0
    kept_rows = 0
    dropped_rows = 0
    malformed_rows = 0

    with open(snakemake.output.maf, "w") as out:
        out.write(header_line + "\n")
        for maf_gz in maf_gz_paths:
            sample_rows = 0
            sample_dropped = 0
            cmd = ["tabix", maf_gz, chrom]
            # stderr=log (not PIPE): nothing here reads a stderr pipe, so a
            # chatty tabix could fill the pipe buffer and deadlock.
            with subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=log, text=True
            ) as proc:
                for line in proc.stdout:
                    sample_rows += 1
                    if vc_idx is not None and coding_classes:
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) <= vc_idx:
                            # Short row -- can't evaluate the filter; keep it
                            # and flag rather than crashing the whole job.
                            malformed_rows += 1
                        elif fields[vc_idx] in coding_classes:
                            sample_dropped += 1
                            continue
                    out.write(line if line.endswith("\n") else line + "\n")

            if proc.returncode != 0:
                print(
                    f"ERROR: tabix exited {proc.returncode} for {maf_gz} "
                    f"region '{chrom}' (see tabix message above)",
                    file=log,
                )
                raise subprocess.CalledProcessError(proc.returncode, cmd)

            if sample_rows == 0:
                # Not fatal on its own, but if every sample reports zero the
                # chromosome naming in CFG['chromosomes'] likely disagrees
                # with the MAFs (e.g. '1' vs 'chr1').
                print(
                    f"  {maf_gz}: 0 rows returned for region '{chrom}' "
                    f"-- check chromosome naming convention",
                    file=log,
                )
            else:
                print(
                    f"  {maf_gz}: {sample_rows} rows on {chrom}, "
                    f"{sample_dropped} dropped (coding)",
                    file=log,
                )

            total_rows += sample_rows
            dropped_rows += sample_dropped
            kept_rows += sample_rows - sample_dropped

    if total_rows == 0:
        print(
            f"WARNING: no rows returned for '{chrom}' across any of the "
            f"{len(maf_gz_paths)} sample MAF(s). Either this chromosome is "
            f"genuinely absent, or the region name does not match the "
            f"contig names in the tabix indexes.",
            file=log,
        )
    if malformed_rows:
        print(
            f"WARNING: {malformed_rows} row(s) had fewer than "
            f"{vc_idx + 1} fields and bypassed the coding filter",
            file=log,
        )

    print(
        f"Total: {total_rows} rows on {chrom} across {len(maf_gz_paths)} samples; "
        f"{dropped_rows} dropped (coding); {kept_rows} written",
        file=log,
    )
    print(f"Wrote {snakemake.output.maf}", file=log)


try:
    main()
except Exception:
    # Snakemake prints the traceback to the job's stderr, but that is not
    # necessarily the same place as {log}; write it here too so the log is
    # self-contained.
    traceback.print_exc(file=log)
    log.flush()
    # Re-raise: swallowing this would mark a failed job as successful and
    # leave a truncated .maf for _mfr_cluster to consume.
    raise
finally:
    log.flush()
    log.close()
