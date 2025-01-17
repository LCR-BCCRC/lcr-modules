# lcr-modules: Standardizing genomic analyses

This project aims to become a collection of standard analytical modules for genomic and transcriptomic data. Too often do we copy-paste from each other’s pipelines, which has several pitfalls. Fortunately, all of these problems can be solved with standardized analytical modules, and the benefits are many. 

**Documentation:** https://lcr-modules.readthedocs.io/

**License:** [LICENSE](LICENSE)

## Installing compatible Snakemake and Oncopipe

Run the following command using `conda` or `mamba` to create the `opv12` environment. 
```bash
conda env create -f demo/env.yaml
```

Always activate this environment before running any pipelines that use LCR-modules. 
```bash
conda activate opv12
```

## Modules

![Module levels](images/module_levels.png)
