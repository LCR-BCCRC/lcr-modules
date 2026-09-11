# pvactools

# Purpose

The `pvactools` module is a Level 2 module wrapping the Griffith Lab's pVACtools (`pvacseq`) ([griffithlab/pVACtools](https://github.com/griffithlab/pVACtools)) for DNA-based neoantigen prediction from somatic SNVs/indels, for `capture, genome` (WES/WGS) data. It generates per-pair filtered/all-epitopes `TSV` reports for MHC Class I and Class II as outputs.

This is an alternative to `modules/neo/1.0` (HMF Tools' NEO), not a replacement -- both modules exist side by side. pVACtools supports many more binding-prediction algorithms and both MHC classes natively, at the cost of a heavier per-pair pipeline (its own dedicated VEP annotation pass, on top of whatever annotation `vcf2maf` already did).

**pVACseq only** for v1 -- not pVACfuse (fusion-derived neoantigens), pVACbind (arbitrary peptide binding), or pVACvector (vaccine design). **DNA-only** for v1 -- no RNA-seq expression/coverage integration (pVACtools treats this as first-class; deferred here, matching this repo's own `lilac`/`mhc_hammer`/`neo` DNA-first precedent).

## Why this module has no variant source of its own

Unlike `neo` (which reads `modules/hmftools/1.1`'s own PURPLE somatic VCF directly), this module is deliberately **caller-agnostic**: it reads whatever VCF `modules/vcf2maf/1.3` is already configured with (`inputs.sample_vcf_gz`, itself fully generic -- "full path to your compressed vcfs from your favourite variant caller"). This means `pvactools` never needs to know or care whether the underlying calls came from SAGE, Strelka, Mutect2, or anything else -- that question is already solved once, upstream, by whatever `vcf2maf` deployment you point this module at.

This also solves a real, separate efficiency problem for free: a WGS sample carries far more somatic variants than WES/capture, and the large majority are synonymous/non-coding -- structurally incapable of producing a neoantigen (a neoantigen requires an altered protein sequence). `vcf2maf`'s own final MAF already has a `Variant_Classification` column from its own earlier VEP run; `_pvactools_extract_coding_regions`/`_pvactools_subset_to_coding_variants` reuse that classification to subset the raw VCF down to only protein-altering positions **before** this module's own VEP+Frameshift+Wildtype annotation and pvacseq's own multi-algorithm binding prediction -- avoiding a full, wasteful re-annotation/re-scoring of variants already known to be silent.

## Prerequisites (read before use)

`vcf2maf` is a **required** upstream dependency -- this module has no variant source without it. Point `inputs.vcf2maf_raw_vcf`/`inputs.vcf2maf_maf` at the same `vcf2maf-1.3` deployment (same `_parent`/`vcf_base_name`/`filter` choices) your project already uses.

`inputs.vcf2maf_maf` must be vcf2maf's own **original, non-CrossMap-projected** MAF (`_vcf2maf_output_original`) -- not its final `_vcf2maf_output_maf` (which is projected onto a single cohort-wide `target_build`). This matters for correctness: `_pvactools_subset_to_coding_variants` matches this MAF's positions directly against `vcf2maf_raw_vcf`'s own coordinates, which are in the sample's native `genome_build`. A CrossMap-projected MAF would silently mismatch for any sample whose native build differs from the projection target.

`modules/mhc_hammer/1.0` is the sole HLA typing source for **both** MHC classes -- `inputs.hla_alleles` (Class I, `_mhc_hammer_output_hla_final_result`) is required, there is no fallback typing source built into this module. `inputs.hla2_alleles` (Class II, `_mhc_hammer_output_hla2_alleles`) is optional -- a pair with none available (typing failed for this patient) still gets a Class-I-only pvacseq run rather than failing. (An earlier draft of this module sourced Class I from `modules/lilac/1.0` instead, for its better-validated WGS typing accuracy -- switched to mhc_hammer for both classes to keep this module to a single upstream typing dependency; worth revisiting if mhc_hammer's own WGS depth-tuning caveat, noted in its own README, turns out to matter in practice.)

VEP is user-supplied (`options.vep_path`, `inputs.vep_cache`), following the same pattern already used by `modules/mhc_hammer/1.0`/`modules/vcf2maf/1.3`, to avoid the bioconda/Perl dependency conflicts already documented for those modules. This module additionally needs two VEP plugins (`Frameshift`, `Wildtype`) that pVACtools bundles and installs itself (`_pvactools_install_vep_plugins`) -- nothing to supply manually for those.

`_pvactools_download_mhcflurry_models` fetches MHCflurry's own trained models (release `2.0.0` by default) once, cohort-wide. MHCnuggets needs no equivalent step -- its trained models ship bundled in the package.

By default, only fully ungated (no license/registration) prediction algorithms are enabled: `MHCflurry`, `MHCflurryEL`, `MHCnuggetsI` (Class I), `MHCnuggetsII` (Class II). Anything IEDB-mediated (`NetMHCpan`, `NetMHCIIpan`, etc.) needs a separately downloaded, license-gated local IEDB installation (`options.iedb_install_directory`) -- the same shape as this repo's own `options.novoalign_dir`/`options.hlahd_dir` gated-tool pattern. Requesting a gated algorithm without setting this fails fast at config-parse time with a clear message, rather than a confusing pvacseq runtime error.

## What's not included in v1

- RNA-seq expression/coverage integration -- pVACtools treats this as first-class; this module doesn't yet wire it in (matches this repo's own DNA-first precedent for `lilac`/`mhc_hammer`/`neo`).
- VAF/depth-based filtering inside pvacseq (`--normal-vaf`/`--tdna-vaf`/etc.) -- pVACtools' own docs describe this as optional, requiring a separate bam-readcount + `vatools` re-annotation pipeline this module doesn't build. A future addition, not required to run pvacseq at all.
- Phased proximal-variant handling (`--phased-proximal-variants-vcf`) -- pVACtools' own docs note the tool this needs (`GATK ReadBackedPhasing`) is no longer available in current GATK versions.
- pVACfuse/pVACbind/pVACvector -- see "Purpose" above.

## Neoantigen MAF output (linking pvactools predictions back to your MAF)

`_pvactools_join_neoantigen_maf`/`_pvactools_crossmap_neoantigen_maf` produce one small, MAF-shaped file per pair (`99-outputs/neoantigen_maf/.../{tumour_id}--{normal_id}--{pair_status}.neoantigen.maf`) -- this MAF's own first `options.minimal_maf_columns` columns (default 45) plus a handful of pvactools annotation columns (`options.neoantigen_maf_columns`: `Tier`, `Best Peptide`, `Allele`, `IC50 MT`, `Pres %ile MT`, `DNA VAF`, and the RNA-related columns kept for forward-compatibility even though always `NA` today), two derived booleans (`Pvacseq_Pass_Relaxed`/`Pvacseq_Pass_VeryRelaxed`, from the `additional_filters` presets above), and any columns defined in `options.neoantigen_pass_filters` (see below). **One row per pvactools-predicted mutation** -- this is not a left join against every row of your cohort MAF, and it does **not** concatenate across a cohort; combining these per-pair files is left to your own downstream MAF-merging tooling.

`options.neoantigen_pass_filters` adds your own combined pass/fail columns on top of the above -- each config key becomes its own new column (named exactly by that key), holding `"True"`/`"False"` per row against that preset's own fixed criteria: `binding_threshold` (max `IC50 MT`), `presentation_percentile_threshold` (max `Pres %ile MT`), `min_dna_vaf` (min `DNA VAF` -- pvacseq's own `binding_filter`, which `Pvacseq_Pass_Relaxed`/`VeryRelaxed` reflect, never looks at DNA VAF at all), and `excluded_tiers` (a list of `Tier` values to reject, empty by default -- `"Subclonal"` is never excluded unless you explicitly add it, see CHANGELOG for why). Omit any of these per preset to skip that criterion; a criterion only auto-passes on a genuinely missing (`"NA"`) value, matching pvactools' own `Filter.execute()` convention. Add, rename, or remove as many presets as you want -- cheap per-row evaluation, no re-prediction, same convention as `options.additional_filters`.

No new upstream dependency: the join (Stage 1) reads `inputs.vcf2maf_maf` -- the same native-build MAF this module already requires -- since `_pvactools_extract_coding_regions` already derives its own coding-region VCF subset from that exact file, every pvactools-predicted mutation is guaranteed a corresponding row in it. Two earlier designs were tried and rejected before landing here: joining directly onto your cohort's final, CrossMap-*projected* MAF (by genomic position, and separately by a `Hugo_Symbol`+`HGVSp_Short` protein-level key) -- see this module's CHANGELOG for the full account, including a real, non-hypothetical case where a cohort's native builds differ from the target build at the majority of positions, which is what actually motivated Stage 2 below.

`options.neoantigen_target_build` (default `"grch37"`) is the single, cohort-wide harmonized build Stage 2 lifts every pair's Stage-1 output into -- lowercase `"grch37"`/`"grch38"`, matching `modules/vcf2maf/1.3`'s own internal convention (not VEP's own `"GRCh37"`/`"GRCh38"` `--assembly` convention used by `options.vep_assembly_map` above). Only `grch37`<->`grch38` lifting is supported, mirroring that module's own CrossMap chain-file limitation (hg19/hg38 chains only) -- a pair whose native build already matches needs no lift at all. The lift itself calls `lcr-scripts/crossmap/1.1/convert_maf_coords.sh` directly (same script/chain files/conda env `modules/vcf2maf/1.3`'s own `_vcf2maf_crossmap` rule uses) but deliberately skips that module's own reannotation step (`_vcf2maf_reannotate`/`maf2maf.pl`) -- gene/protein identity doesn't depend on which build's coordinates it's expressed in, and reannotating would reintroduce a transcript-instability risk (see CHANGELOG). That script's own multi-stage shell pipeline was found to be intermittently flaky (real, observed ~1-in-10 failure rate, see CHANGELOG) -- `_pvactools_crossmap_neoantigen_maf` retries it up to 3 times to work around this.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    pvactools:
        inputs:
            vcf2maf_raw_vcf: "results/vcf2maf-sage-1.1_vcf2maf-1.3/00-inputs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/combined.passed.vcf.gz"
            vcf2maf_maf: "results/vcf2maf-sage-1.1_vcf2maf-1.3/99-outputs/deblacklisted/maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.combined.passed.maf"
            hla_alleles: "results/mhc_hammer-1.0/99-outputs/hla_alleles/{seq_type}--{genome_build}/{patient_id}.hla_alleles.txt"
            # Optional -- omit or leave as "" to run Class-I-only
            hla2_alleles: "results/mhc_hammer-1.0/99-outputs/hla2_alleles/{seq_type}--{genome_build}/{patient_id}.hla2_alleles.csv"
            vep_cache: "ref/ensembl_vep_cache/"
        options:
            vep_path: "/path/to/vep/bin"
            # Only relevant to the neoantigen MAF feature -- see above. No extra inputs needed:
            # it reuses vcf2maf_maf above directly.
            neoantigen_target_build: "grch37"
```

The example snakefile:

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

configfile: "../modules/pvactools/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/pvactools/1.0/pvactools.smk"

rule all:
    input:
        rules._pvactools_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/pvactools/CHANGELOG.md)
