# mhc_hammer

# Purpose

The `mhc_hammer` module is a Level 2 module that operates on `BAM` files to use MHC Hammer ([McGranahanLab/mhc-hammer](https://github.com/McGranahanLab/mhc-hammer)) and detect HLA class I gene disruption (somatic mutations, loss of heterozygosity, and allelic imbalance) for `capture, genome` (WES/WGS) data. It generates `CSV` files as outputs.

This is a **DNA-only v1** port of MHC Hammer's DNA analysis arm (HLA typing, personalised HLA reference construction, Novoalign alignment, allele-specific BAM splitting, copy-number/allelic-imbalance/LOH detection, and Mutect2+VEP mutation calling), reimplemented as native Snakemake rules instead of wrapping upstream's Nextflow pipeline. Upstream's RNA allelic expression/imbalance/repression analysis and the alternative-splicing arm are **not** included in this version.

In addition to the HLA class I disruption-detection pipeline above (which is entirely upstream MHC Hammer's own scope), this module also includes a separate, parallel **HLA class II germline typing** path (`_mhc_hammer_hla2_*` rules) -- not part of upstream MHC Hammer at all. It types the classical peptide-presenting genes (`DRA`, `DRB1`, `DRB3`, `DRB4`, `DRB5`, `DQA1`, `DQB1`, `DPA1`, `DPB1`) plus the non-classical peptide-loading/editing genes `DMA`, `DMB`, `DOA`, `DOB` (`options.hla2_genes`) from each patient's germline sample via HLA-HD, using this module's own `src/parse_hlahd_output.R` script rather than upstream's `bin/hlahd_parse_output.R` (which needs a class-I-only GTF this path has no equivalent of -- see the comment at the top of that script). This is **germline genotyping only** -- no personalised reference, Novoalign alignment, copy-number/allelic-imbalance, or mutation calling for class II, unlike the class I path. Output: `99-outputs/hla2_alleles/{seq_type}--{genome_build}/{patient_id}.hla2_alleles.csv`.

## Prerequisites (read before use)

This module requires three things the user must obtain and configure separately -- it cannot install them automatically:

1. **Novoalign** (`options.novoalign_dir`): download from [novocraft.com](https://www.novocraft.com/support/download/) (only V3.x and earlier are free for non-commercial/academic use) and point this at the directory containing the `novoalign`/`novoindex` binaries.
2. **HLA-HD** (`options.hlahd_dir`): request a download from the [HLA-HD website](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/download-request/), install it, and point this at the installed directory (must contain `bin/`, `dictionary/`, `freq_data/`, `HLA_gene.split.txt`). Its database version must match `options.imgt_release` below.
3. **MHC Hammer scripts** (`options.mhc_hammer_scripts_dir`): MHC Hammer is distributed by Cancer Research Horizons under an **academic-use licence that prohibits redistribution and modification** (see the [LICENSE](https://github.com/McGranahanLab/mhc-hammer/blob/main/LICENSE) in the upstream repo). This module therefore never bundles or copies any of upstream's `bin/*.R`/`bin/*.sh` scripts -- it calls them from a path you provide. `git clone https://github.com/McGranahanLab/mhc-hammer.git` yourself, read and accept its licence, and point `mhc_hammer_scripts_dir` at that clone. If you publish results produced with this module, cite the MHC Hammer publication per the licence's attribution requirement.
4. **VEP** (`options.vep_path`, `inputs.vep_cache`): also user-supplied, following this repo's own `vcf2maf` module convention, to avoid the bioconda/Perl dependency conflicts already hit there (see `CHANGELOG.md`).
5. **A reference genome FASTA** (`inputs.reference_genome_fasta`): any complete human genome FASTA, doesn't need to match any sample's `genome_build` -- used only to identify non-specific/repetitive kmers genome-wide.

Every patient in your sample table must have **exactly one** germline WES sample (`tissue_status: normal`) -- unlike upstream, this module does not tolerate multiple germline samples per patient by silently picking one, and it will error instead. Tumour WES samples without a matched germline in the same patient are never processed (see `CFG["paired_runs"]` in the module code) since HLA typing and the personalised reference both require the patient's own germline sample.

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
