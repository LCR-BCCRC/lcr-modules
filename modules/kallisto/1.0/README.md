# LCR Module: kallisto

## Authors

- **Module:** Helena Winata
- **Original Snakefile:** Helena Winata

## Inputs
Paired-run FASTA files.

## Outputs
`abundances.tsv` - Estimated gene counts
`abundances.h5` - HDF5 file
`run_info.json` - run information

### Optional
`pseudoalignments.bam` - BAM file generated when using `--pseudobam` flag


## Requirements

### Software

<!-- TODO: If you can't use conda, update this section accordingly -->

#### Conda

This module automatically installs software requirements using conda. The core packages for each conda environment in `envs/kallisto-0.46` are listed below.

kallisto v0.46

```bash
# kallisto environment
TODO
```

#### Singularity

This module does not rely on Singularity.

### Reference data

<!-- TODO: Add required references below (including the reference key) -->
References can be generated using `kallisto index`
- Kallisto index 
- Kallisto gtf file and chromosomes list

## Release Notes

### Version 1.0 (Helena Winata)

<!-- TODO: Add items below explaining each decision -->

- TODO
