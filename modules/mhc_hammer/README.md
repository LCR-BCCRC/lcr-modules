# mhc_hammer

# Purpose

The `mhc_hammer` module is a Level 2 module that operates on `BAM` files to use MHC Hammer ([McGranahanLab/mhc-hammer](https://github.com/McGranahanLab/mhc-hammer)) and detect HLA class I gene disruption (somatic mutations, loss of heterozygosity, and allelic imbalance) for `capture, genome` (WES/WGS) data. It generates `CSV` files as outputs.

This is a **DNA-only v1** port of MHC Hammer's DNA analysis arm (HLA typing, personalised HLA reference construction, Novoalign alignment, allele-specific BAM splitting, copy-number/allelic-imbalance/LOH detection, and Mutect2+VEP mutation calling), reimplemented as native Snakemake rules instead of wrapping upstream's Nextflow pipeline. Upstream's RNA allelic expression/imbalance/repression analysis and the alternative-splicing arm are **not** included in this version.

In addition to the HLA class I disruption-detection pipeline above (which is entirely upstream MHC Hammer's own scope), this module also includes a separate, parallel **HLA class II germline typing** path (`_mhc_hammer_hla2_*` rules) -- not part of upstream MHC Hammer at all. It types the classical peptide-presenting genes (`DRA`, `DRB1`, `DRB3`, `DRB4`, `DRB5`, `DQA1`, `DQB1`, `DPA1`, `DPB1`) plus the non-classical peptide-loading/editing genes `DMA`, `DMB`, `DOA`, `DOB` (`options.hla2_genes`) from each patient's germline sample via HLA-HD, using this module's own `src/parse_hlahd_output.R` script rather than upstream's `bin/hlahd_parse_output.R` (which needs a class-I-only GTF this path has no equivalent of -- see the comment at the top of that script). This is **germline genotyping only** -- no personalised reference, Novoalign alignment, copy-number/allelic-imbalance, or mutation calling for class II, unlike the class I path. Output: `99-outputs/hla2_alleles/{seq_type}--{genome_build}/{patient_id}.hla2_alleles.csv`.

Also included: `_mhc_hammer_mutations_to_maf` (`src/mutations_to_maf.R`, module-owned), which reformats each patient's `{patient_id}_mutations.csv` into a MAF-like table for tools that expect MAF input (e.g. `maftools`) but don't need real genomic coordinates. Mutations are called against each patient's own personalised, allele-specific reference (one contig per typed HLA allele, e.g. `A*02:01:01:01`), so there is no genomic coordinate system shared across alleles or samples -- `Chromosome`/`Start_Position`/`End_Position`/`NCBI_Build` in the output are therefore synthetic placeholders (see the header comment in `src/mutations_to_maf.R`), while `Hugo_Symbol`, `Variant_Classification`, `Variant_Type`, allele calls, VEP annotations, and read counts are all derived from the real underlying data. The original mhc_hammer-relative allele/position are preserved in the non-standard `HLA_Allele`/`HLA_Allele_Position` columns. Only PASS-filtered mutations (`Filter == "PASS"`) are included. Output: `99-outputs/mutations_maf/{seq_type}--{genome_build}/{patient_id}.mutations.maf`.

HLA-HD's own consolidated HLA class I result file (`{patient_id}_final.result.txt`, in HLA-HD's native standard IMGT/HLA nomenclature -- unlike the mhc_hammer-internal personalised-reference contig names, e.g. `hla_a_01_01_01_01`, used elsewhere in this module's own outputs including the DNA analysis and mutations tables) is also symlinked to `99-outputs/hla_alleles/{seq_type}--{genome_build}/{patient_id}.hla_alleles.txt` for convenience, since it's otherwise only reachable inside the internal `03-hlahd/` working directory.

## Prerequisites (read before use)

This module requires three things the user must obtain and configure separately -- it cannot install them automatically:

1. **Novoalign** (`options.novoalign_dir`): download from [novocraft.com](https://www.novocraft.com/support/download/) (only V3.x and earlier are free for non-commercial/academic use) and point this at the directory containing the `novoalign`/`novoindex` binaries.
2. **HLA-HD** (`options.hlahd_dir`): request a download from the [HLA-HD website](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/download-request/), install it, and point this at the installed directory (must contain `bin/`, `dictionary/`, `freq_data/`, `HLA_gene.split.txt`). Its database version must match `options.imgt_release` below.
3. **MHC Hammer scripts** (`options.mhc_hammer_scripts_dir`): MHC Hammer is distributed by Cancer Research Horizons under an **academic-use licence that prohibits redistribution and modification** (see the [LICENSE](https://github.com/McGranahanLab/mhc-hammer/blob/main/LICENSE) in the upstream repo). This module therefore never bundles or copies any of upstream's `bin/*.R`/`bin/*.sh` scripts -- it calls them from a path you provide. `git clone https://github.com/McGranahanLab/mhc-hammer.git` yourself, read and accept its licence, and point `mhc_hammer_scripts_dir` at that clone. If you publish results produced with this module, cite the MHC Hammer publication per the licence's attribution requirement.
4. **VEP** (`options.vep_path`, `inputs.vep_cache`): also user-supplied, following this repo's own `vcf2maf` module convention, to avoid the bioconda/Perl dependency conflicts already hit there (see `CHANGELOG.md`).
5. **A reference genome FASTA** (`inputs.reference_genome_fasta`): any complete human genome FASTA, doesn't need to match any sample's `genome_build` -- used only to identify non-specific/repetitive kmers genome-wide.

Every patient in your sample table must have **exactly one** germline WES sample (`tissue_status: normal`) -- unlike upstream, this module does not tolerate multiple germline samples per patient by silently picking one, and it will error instead. Tumour WES samples without a matched germline in the same patient are never processed (see `CFG["paired_runs"]` in the module code) since HLA typing and the personalised reference both require the patient's own germline sample.

## Opt-in: RNA-seq input for HLA-HD typing

HLA typing (Class I and Class II) can optionally also run from RNA-seq, via two independent, default-off options -- both `False` by default, zero behaviour change unless explicitly enabled:

- **`options.rna_hla_typing_fallback`**: type HLA from RNA-seq for patients who have RNA but no germline WES/WGS sample at all -- currently such patients are fully excluded from this module entirely. **Does not extend the DNA-analysis arm** (personalised reference, Novoalign, copy-number/allelic-imbalance, mutation calling) -- those still require a real germline+tumour WES pair; a fallback-typed patient gets an HLA call and nothing else from this module.
- **`options.rna_hla_typing_comparison`**: additionally runs RNA-based typing for patients who *do* have a real DNA pair, alongside (not instead of) the existing DNA-based typing, so results can be cross-checked.

This uses `inputs.sample_rna_bam`/`inputs.sample_rna_bai` (same three wildcards as `sample_bam`/`sample_bai`) as the RNA BAM source -- typically `modules/star/1.4`'s own output, but any coordinate-sorted RNA-seq BAM works. RNA-seq input for HLA-HD is an officially supported use case, not a hack: per [HLA-HD's own site](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/), "RNA-Seq data can also be applied" -- same CLI, same exon/intron dictionary reconciliation, just fed different FASTQs.

**Real caveat driving the sample-selection logic**: this module's whole DNA-analysis arm exists to detect HLA LOH. RNA-based germline typing sourced from the *tumour* being analysed for LOH risks a genuine circularity -- a lost allele wouldn't be transcribed and could be mistyped as homozygous. `_mhc_hammer_get_rna_typing_sample` therefore prefers a `tissue_status: normal` RNA sample; only falls back to tumour/tumor RNA if no normal RNA exists for that patient (picking deterministically, by sample_id, if multiple tumour RNA samples exist -- e.g. from different biopsies/timepoints). `99-outputs/hla_typing_source_rna/{seq_type}--{genome_build}/{patient_id}.hla_typing_source.csv` always records which sample and tissue_status was actually used, and flags `tumour_derived_loh_risk` explicitly -- check this file before trusting an RNA-sourced call for LOH-sensitive analysis.

Other RNA-specific outputs: `99-outputs/hla_alleles_rna/` and `99-outputs/hla2_alleles_rna/` (same format as the DNA-sourced `hla_alleles`/`hla2_alleles`, so anything already consuming those files can consume these too). `options.hlahd_rna_min_read_length` (independent of `options.hlahd_min_read_length`, since RNA read lengths may differ from this cohort's own DNA reads) controls HLA-HD's own `-m` minimum tag size for the RNA-sourced typing rules specifically.

## What's not included in v1

- RNA allelic expression, RNA allelic imbalance, and RNA allelic repression (tumour vs. matched-normal RNA) -- upstream's RNA analysis arm.
- Alternative splicing detection (2-pass STAR + splice junction analysis).
- Upstream's graceful per-patient exclusion on HLA-HD failure -- this module fails loudly instead.
- The `exon_snps`-restricted variant of copy-number/allelic-imbalance detection (upstream itself has this disabled).
- BAM-subsetting bypass / pre-typed-HLA-input / preprocessing-only modes.
- HLA class II somatic disruption detection (personalised reference, Novoalign, copy-number/allelic-imbalance, mutation calling) -- only germline typing is implemented for class II (see above); this is left as a possible future extension of the parallel `_mhc_hammer_hla2_*` path.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    mhc_hammer:
        inputs:
            sample_bam: "data/{sample_id}.bam"
            vep_cache: "/path/to/vep/cache"
            reference_genome_fasta: "/path/to/any/reference/genome.fa"
        options:
            novoalign_dir: "/path/to/novocraft"
            hlahd_dir: "/path/to/hlahd"
            mhc_hammer_scripts_dir: "/path/to/your/mhc-hammer/clone"
            vep_path: "/path/to/vep/bin"
```

The example snakefile:

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

configfile: "../modules/mhc_hammer/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/mhc_hammer/1.0/mhc_hammer.smk"

rule all:
    input:
        rules._mhc_hammer_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/mhc_hammer/CHANGELOG.md)
