# neo

# Purpose

The `neo` module is a Level 3 module that predicts and scores tumour neoepitopes using HMF Tools' NEO ([hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools/tree/master/neo)), combining somatic mutation calls, LINX fusion neoepitopes, and LILAC HLA typing into a ranked set of candidate MHC class I neoepitope peptides for `capture, genome` (WES/WGS) data.

NEO is a genuinely cross-tool downstream analysis, not a preprocessing detail of any one upstream tool -- it structurally requires:
- **Somatic variants**: PURPLE's own somatic VCF (variant copy number and subclonal likelihood, real, standard PURPLE VCF annotations -- confirmed to *not* require the same `cnv.gene.tsv` gene-copy-number file whose schema drift broke `modules/lilac/1.0`'s own `gene_copy_number` input).
- **Gene fusion neoepitopes**: LINX's own `-write_neo_epitopes` output.
- **HLA typing**: `modules/lilac/1.0`'s own allele coverage/CN output.

This is exactly why `neo` is its own standalone module rather than being folded into `modules/lilac/1.0` -- its real inputs are a cross-module-input problem (the same convention `lilac` itself already uses for `gene_copy_number`/`somatic_vcf`), not something any one of those other modules' own internal steps needs. A standalone module can point at *any* compatible PURPLE/LINX/LILAC output, not just this repo's own.

RNA/expression annotation (Isofox, cohort TPM medians) is out of scope for v1 -- DNA-only mutation-derived predictions only, matching the same DNA-only precedent already set for `modules/mhc_hammer/1.0` and `modules/lilac/1.0`. NEO's own README confirms this degrades gracefully (un-adjusted presentation likelihood, no expression annotation) rather than failing.

## Prerequisites (read before use)

This module needs real output from **three other modules/pipeline steps**, none of which it can produce itself:

1. **`modules/hmftools/1.1`'s own PURPLE somatic VCF**, exposed via a new rule (`_hmftools_output_purple_somatic_vcf`) added to that module specifically for this -- PURPLE already writes this file as a side-effect of being given `-somatic_vcf`, it just wasn't tracked before. No config changes needed in `hmftools/1.1` itself beyond having this module updated.
2. **`modules/hmftools/1.1`'s own LINX neoepitope output**, via a new, separate rule (`_hmftools_linx_neo_epitopes`) added to that module -- deliberately *not* merged into the existing `_hmftools_linx` rule, so requesting it never forces an already-completed sample's LINX run to redo. Confirmed the pinned `hmftools-linx=1.15` already supports `-write_neo_epitopes` (no LINX version bump needed).
3. **`modules/lilac/1.0`'s own HLA typing output** (`inputs.lilac_tsv_file`, pointing at `lilac`'s *internal* `02-lilac/` directory, not `99-outputs/` -- see the long comment on that config value for why).

`_neo_download_reference` fetches NEO's own bundled scoring reference data (position-weight matrices, likelihood distributions) from `www.bcgsc.ca/downloads/morinlab/hmftools-references/neo/`, same pre-extracted mirror convention as `modules/lilac/1.0`/`modules/hmftools/1.1`. The exact 8 filenames (`neo_train_pos_weight.csv`, `neo_train_flank_pos_weight.csv`, `neo_train_likelihood.csv`, `neo_train_expression_dist.csv`, `neo_train_recognition.csv`, `neo_train_rand_dist.csv`, `neo_train_likelihood_rand_dist.csv`, `neo_train_exp_likelihood_rand_dist.csv`) are now confirmed against both `BindCommon.formFilename()` and the real, current `hartwigmedical/hmftools` `pipeline/README_RESOURCES.md` NEO section -- source these 8 files from there (via the [Oncoanalyser reference-data page](https://nf-co.re/oncoanalyser/docs/usage/#reference-data-urls) it links to) and place them, unrenamed, at that mirror path. `-score_file_id` is deliberately never passed to `NeoScorer` -- see `CHANGELOG.md` for why.

This module runs its own PAVE annotation step (`_neo_pave_annotate`) on the raw PURPLE somatic VCF before NEO ever sees it -- **not optional**: NEO reads the same `IMPACT` VCF INFO tag `modules/lilac/1.0` needed (confirmed by reading `PointMutationData.isRelevantMutation()`/`NeoSampleTask.java` directly), and a raw, unannotated VCF would silently exclude every variant with no error. `_neo_download_ensembl_cache` fetches PAVE's required Ensembl gene/transcript cache -- reusable via `inputs.ensembl_data_dir` if you already have a copy (e.g. from `modules/hmftools/1.1`), same pattern as `modules/lilac/1.0`.

## Real CLI verification notes

Confirmed directly against the real `neo-v1.3` release jar (downloaded from GitHub, run locally with `java -cp`/the installed `neo` wrapper -- same technique already used for `lilac`/`pave`/`purple`):

- The module's own README documents stale, nonexistent class names -- the real entry points are `com.hartwig.hmftools.neo.epitope.NeoEpitopeFinder` (Step 1) and `com.hartwig.hmftools.neo.score.NeoScorer` (Step 2), not `neo.cohort.NeoEpitopeAnnotator`/`neo.scorer.NeoScorer`.
- The installed `neo` wrapper script dispatches to either entry point from a single binary: pass the fully-qualified class name as the first argument to switch from `-jar neo.jar` (`NeoEpitopeFinder`, the default) to `-cp neo.jar <class>` (`NeoScorer`) -- confirmed `-Xmx` placement doesn't matter either way.
- `NeoScorer`'s own registered CLI config has zero CLI-flagged-`REQUIRED` items at all -- its true requirement `-score_file_dir` is enforced at runtime; confirmed via a real crash (`NullPointerException` in `BindScoreMatrix.loadFromCsv`) without it.
- `-score_file_id`'s real effect on filename construction was initially guessed wrong (`ScoreConfig.java`/`TrainConfig.java` alone don't show it) -- reading `BindCommon.formFilename()` directly showed the real pattern is `neo_<fileType>` id-less (`neo_train_pos_weight.csv`) when no id is given, or `neo_<fileType>_<id>.csv` (underscore, not the originally-assumed `.`) when one is. Since the real, standard HMF resource bundle (per `pipeline/README_RESOURCES.md`) ships the id-less filenames, this module never passes `-score_file_id` at all -- confirmed working with `_neo_download_reference`'s own (now-corrected) filenames.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    neo:
        inputs:
            linx_neo_epitope_file: "results/hmftools-1.1/06-linx_neo_epitopes/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.linx.neoepitope.tsv"
            somatic_vcf: "results/hmftools-1.1/99-outputs/purple_somatic_vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.purple.somatic.vcf.gz"
            lilac_tsv_file: "results/lilac-1.0/02-lilac/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.lilac.tsv"
            purple_purity_file: "results/hmftools-1.1/04-purple/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.purple.purity.tsv"
```

The example snakefile:

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

configfile: "../modules/hmftools/1.1/config/default.yaml"
configfile: "../modules/lilac/1.0/config/default.yaml"
configfile: "../modules/neo/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/hmftools/1.1/hmftools.smk"
include: "../modules/lilac/1.0/lilac.smk"
include: "../modules/neo/1.0/neo.smk"

rule all:
    input:
        rules._hmftools_all.input,
        rules._lilac_all.input,
        rules._neo_all.input,
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/neo/CHANGELOG.md)
