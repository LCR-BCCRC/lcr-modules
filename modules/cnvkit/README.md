# CNVkit module

A lightweight CNV calling module built around **CNVkit** with optional B‑allele fraction (BAF) integration. It builds capture‑aware bins, normalizes with a panel of normals (PoN), segments, calls copy number, and emits plots, gene‑level summaries, and SEG files—optionally projected to hg19/GRCh37 or hg38/GRCh38.

---

## Is this the right tool for you?
Use this module if you need:
- CNV calls from WES/targeted tumour (± matched normal) BAMs.
- PoN‑based normalization with capture‑aware targets/antitargets.
- BAF‑aware CNV **calling** using dbSNP sites.
- Clean plots (scatter + diagram), SEG export, basic QC (MAD, IQR, #segments).
- Optional lift‑over and projection‑normalized SEG outputs with consistent `chr` prefixing.

You may want a different tool if you need whole‑genome read‑depth CNV callers, deep purity/ploidy modeling, or allele‑specific CN across the whole genome without dbSNP reliance.

---

## What it produces
Per tumour sample 

- `*.cns` (BAF‑aware **called** segments)
- `*_scatter.png` and `*_diagram.pdf`
- `*.genebreaks.txt`, `segment.gene_cn.txt`, `trusted_genes.txt`, `inferred_sex.txt`
- `*.seg` (native build) and **projection** SEG(s) if requested
- `*.metrics.txt` (MAD, IQR, segment count)

---

## Overview

Essential inputs:
- BAM/BAI per sample (tumour ± normal)
- Genome FASTA and capture BED (or default exome UTRs)
- dbSNP common sites VCF (for BAF; e.g., build 151)

---

## Minimal config (example)
```yaml
lcr-modules:
  cnvkit:
    options:
      new_normals: true           # rebuild PoN when normals change
      male_ref: "--male-reference"
      cns: { method: "cbs" }
      segmetrics: { add_col: "--ci --t-test" }
      BAF: { rescale: "threshold", min_depth: 20, filter_by: "ci" }
    output:
      requested_projections: ["hg38","grch37"]
```
> If `capture_space` isn’t provided, defaults to an exome‑UTR bed per build. You can override per‑capture `target_bed` in config.

---

## Outline of the workflow
1. Build **accessible regions** and capture‑aware target/antitarget bins.
2. Compute **coverage**; (optionally) build **PoN**.
3. **Fix** → **segment**.
4. Create per‑chrom **mpileups** on dbSNP; **concat** VCF.
5. Add **segmetrics** (t‑test & CIs) and **call** with VCF (BAF‑aware).
6. Generate **plots**, **gene metrics**, **sex**, **breaks**, **metrics**, **SEG**.
7. Optional **lift‑over**, **fill neutral regions**, and **projection** normalization.

---

## Key options to know
- `options.new_normals`: `true` to rebuild PoN; `false` to reuse existing PoN.
- `options.BAF.rescale`: `"threshold"` by default (useful when purity/ploidy unknown).
- `output.requested_projections`: e.g., `["hg38","grch37"]` to emit normalized SEGs.

---

## Changelog (highlights)
- **1.0 (2022‑11‑26):** Add `metrics` rule (MAD, IQR, #segments).
- **1.0 (2022‑02‑17):** Split batch into modular rules; `new_normals` option; default exome fallback; BAF from dbSNP; CI‑based filtering; plots and SEG export.
- **1.0 (2022‑01‑22):** Initial release (Jasper Wong), CNVkit batch conventions (.cnn/.cnr/.cns/.call.cns).

---

## Attribution & licenses
- **CNVkit** and **bcftools** retain their respective licenses.
- This module is part of the LCR-modules ecosystem; module authored by **Jasper Wong**.
