# CNAqc Module

This module runs [CNAqc](https://github.com/caravagnalab/CNAqc) (Copy Number Alteration quality check) 
to assess the consistency between somatic mutations and copy number data.

## Inputs

| Name | Description |
|------|-------------|
| `subclones` | Battenberg 1.2 filled subclones file (`_subclones.txt`) |
| `cellularity_ploidy` | Battenberg cellularity/ploidy file (`_cellularity_ploidy.txt`) |
| `sample_maf` | Somatic SNV MAF file (e.g. from slms-3 via vcf2maf) |

## Outputs

| Type | Description |
|------|-------------|
| PDF | Per-sample QC plots (CNA overview, VAF histograms, peak analysis) |
| TSV | Per-sample CNAqc score and QC metrics |

## Configuration

Key options in `default.yaml`:

- `min_depth`: Minimum read depth for a variant to be included (default: 10)
- `min_muts_per_segment`: Minimum mutations per segment for peak analysis (default: 10)
- `run_peak_analysis`: Whether to run peak analysis (default: True)

## Typical Input Paths

**Battenberg 1.2 outputs:**
```yaml
subclones: "results/battenberg-1.2/99-outputs/txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_subclones.txt"
cellularity_ploidy: "results/battenberg-1.2/99-outputs/txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_cellularity_ploidy.txt"
```

**vcf2maf 1.3 MAF outputs (slms-3 pipeline):**
```yaml
sample_maf: "results/vcf2maf-1.3/99-outputs/deblacklisted/maf/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.slms-3.{genome_build}.maf"
```

## Reference

Caravagna G, et al. (2022). Reliable tumor purity and cellularity estimation from high-depth sequencing data.
https://github.com/caravagnalab/CNAqc
