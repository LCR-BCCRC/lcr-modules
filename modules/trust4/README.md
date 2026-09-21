# trust4

# Purpose

The `trust4` module is a Level 2 module that operates on `BAM` files to use TRUST4 ([liulab-dfci/TRUST4](https://github.com/liulab-dfci/TRUST4)) and reconstruct BCR/TCR immune repertoires (V/D/J/C genes, CDR3 sequences) from bulk RNA-seq for `mrna` data. It generates `TSV`/`FASTA` files as outputs.

This is a **bulk RNA-seq v1**: it takes an existing coordinate-sorted, indexed BAM as input (by default, `modules/star/1.4`'s own final output) and does not run TRUST4's own barcode/UMI/10x single-cell machinery, nor its raw-FASTQ input mode. TRUST4 operates on one BAM per `sample_id` -- there is no tumour/normal pairing concept anywhere in this module.

TRUST4 is bioconda-installable (`trust4=1.1.10`) with no licensing restriction, unlike `modules/mhc_hammer/1.0`'s Novoalign/HLA-HD/VEP dependencies -- this module installs it via a normal conda/container environment, no user-supplied external script directory needed.

**Input BAMs are content-sniffed, not trusted by extension.** TRUST4's own `bam-extractor` binary links against a bundled, ancient `samtools-0.1.19` (predates CRAM support entirely) -- a real production run hit this directly: a sample's `inputs.sample_bam` was actually CRAM content stored with a `.bam` name (a real, confirmed storage convention on this cohort's cluster, already documented for DNA BAMs in `modules/mhc_hammer/1.0`), and `bam-extractor` crashed outright (`invalid BAM binary header (this is not a BAM file)`) rather than reading it. `_trust4_input_bam` checks the real first 4 bytes: a genuine BAM is symlinked as before (cheap, no wasted duplicate); genuine CRAM is decoded via `samtools view -b -T <genome_fasta>` first (a real htslib tool, so it handles the decode correctly regardless of what upstream naming convention produced it); anything that's neither fails loudly with a clear message rather than being silently mishandled. Whatever lands in `00-inputs/` (symlink or real converted copy) is `temp()`-marked -- `_trust4_run` is its only consumer, so a real converted duplicate doesn't sit on disk forever, and the original file you point `inputs.sample_bam` at is never touched either way. No `inputs.sample_bai` needed any more -- this rule indexes its own output directly.

## Reference files

TRUST4's own GitHub repo ships small, pre-built, plain-text reference FASTAs directly in-repo (not via bioconda): a genome-coordinate-aware V/D/J/C gene FASTA per genome build (`hg19_bcrtcr.fa`/`hg38_bcrtcr.fa`), and a species-only IMGT reference (`human_IMGT+C.fa`). This module downloads both automatically, pinned to a specific commit (`options.trust4_repo_commit`) for reproducibility -- no manual reference-building step needed for standard human genome builds (`grch37`/`hg19`/`grch38`/`hg38`, mapped via `options.trust4_bcrtcr_map`).

**`options.abnormal_unmap_flag`** (default `False`): a real production crash -- `bam-extractor` assumes (confirmed directly from its own source) that a completely unmapped read pair's two mates appear as *consecutive* records in a coordinate-sorted BAM, and hard-fails (`Two reads from the unaligned fragment are not showing up together`) if they don't. Some aligners/alignment parameter sets don't guarantee this placement. Set `True` if your own BAMs hit this crash -- it maps directly to TRUST4's own `--abnormalUnmapFlag` (-> `bam-extractor`'s `-u`), which disables that adjacency assumption entirely. **This is currently a cohort-wide toggle, not per-sample** -- setting it affects every sample's `_trust4_run` invocation, not just the one that hit the crash.

## Opt-in: BAM-reconstructed-FASTQ comparison path

`_trust4_run_fastq` (and its own `_trust4_bam_to_fastq` upstream step) run TRUST4 in its native FASTQ mode (`-1/-2`, via `fastq-extractor`) instead of BAM mode (`-b`, via `bam-extractor`), on FASTQ reconstructed from the same already-normalized BAM via `samtools collate | samtools fastq`. This sidesteps the entire class of BAM/CRAM-convention issue this module has hit twice (CRAM disguised as `.bam`; an aligner that doesn't place unmapped-pair mates adjacently) -- `fastq-extractor` has no such assumptions at all.

**Real caveat**: these FASTQs are reconstructed from an *already-aligned* BAM, not the sample's true original sequencer output -- they inherit whatever alignment-time decisions already happened (secondary/supplementary alignment filtering, soft-clip handling). This is a comparison against "the same reads, different extraction code path," not a genuinely independent rerun from raw data.

**Deliberately not part of `_trust4_all`'s default target list** -- it would double every sample's compute by default. Request it per-sample directly:
```
snakemake ... results/trust4-1.0/99-outputs/report_reconstructed_fastq/{seq_type}--{genome_build}/{sample_id}.trust4_reconstructed_fastq_report.tsv
```
(and the `cdr3_reconstructed_fastq`/`annot_reconstructed_fastq`/`airr_reconstructed_fastq` siblings, same naming pattern). `options.repseq`/`options.abnormal_unmap_flag` still apply in this mode.

## What's not included in v1

- Barcode/UMI/10x single-cell support (`--barcode`, `--barcodeWhitelist`, `--UMI`, and the resulting `trust-barcoderep.pl`-produced `_barcode_report.tsv`/`_barcode_airr.tsv` outputs) -- bulk RNA-seq only.
- Mouse genome builds (TRUST4 also ships `GRCm38_bcrtcr.fa`/`GRCm39_bcrtcr.fa`/`mouse_IMGT+C.fa`) -- a trivial addition to `options.trust4_bcrtcr_map` if ever needed, not wired up now.
- `--assembleWithRef`/`--stage`/`--outputReadAssignment` and other advanced TRUST4 options -- not exposed as module options; open an issue/PR if you need one.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    trust4:
        inputs:
            # Defaults to modules/star/1.4's own output -- override only if your own
            # deployment uses a different star version/path, or a different RNA-seq aligner
            sample_bam: "star-1.4/99-outputs/bam/{seq_type}--{genome_build}/{sample_id}.bam"
```

The example snakefile:

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

configfile: "../modules/trust4/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/trust4/1.0/trust4.smk"

rule all:
    input:
        rules._trust4_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/trust4/CHANGELOG.md)
