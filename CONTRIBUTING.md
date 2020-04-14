# Contribution Guidelines for LCR Modules

This document provides a set of guidelines for contributing to the `lcr-modules` repository of standard analytical pipelines for genomic and transcriptomic data.

These guidelines rely on Git and GitHub for features like branching and pull requests. If you want to learn Git, I recommend working your way through [this tutorial](https://hamwaves.com/collaboration/doc/rypress.com/index.html).

## Adding a New Module

### Getting set up

**Note:** Don't forget to update any values in angle brackets (`<...>`).

1. First, check the [repository](https://github.com/LCR-BCCRC/lcr-modules/tree/master/modules/) to see if the module already exists. Then, check the [list of open issues](https://github.com/LCR-BCCRC/lcr-modules/issues?q=is%3Aopen+is%3Aissue+label%3Anew-module) to see if the module has already been proposed. If it has and hasn't been assigned yet, feel free to assign yourself. Otherwise, reach out to the assignee and find out how you can help.

   **Important:** Please make sure you're assigned to a GitHub issue before you start developing the module to avoid duplicated efforts.

2. Install `snakemake` (5.4 or later), `pandas`, `cookiecutter`, and the custom `modutils` Python packages into your conda environment. If you also need Git, you can add `git` to the `conda install` command.

   ```bash
   # `snakemake-minimal` lacks extraneous dependencies and is faster to install
   conda install --satisfied-skip-solve 'snakemake-minimal>=5.4' 'pandas' 'cookiecutter'
   # This is installing from the `lcr-modules` repository clone
   pip install -e lcr-modules/modutils
   ```

3. Clone the `lcr-modules` repository **recursively** if you don't already have a local copy.

   ```bash
   git clone --recursive https://github.com/LCR-BCCRC/lcr-modules.git
   ```

4. Create a new branch from the `master` branch with the format `module/<tool_name>/1.0`. For example, the branch name for contributing a STAR alignment module would be `module/star/1.0`.

   **Important:** Your `<tool_name>` should only contain alphanumerical characters or underscores (_i.e._ no spaces).

   ```bash
   cd lcr-modules
   git checkout master
   git pull --ff-only
   git checkout -b "module/<tool_name>/1.0"
   ```

   You can confirm which branch you are on by running `git branch`.

5. Initialize your module with the template with the following command:

   ```bash
   cookiecutter template --output-dir 'modules/'
   ```

   You will be asked for the following information:

   - `module_name`: Short name consisting of only alphanumerical characters or underscores (_i.e._ no spaces).
   - `module_author`: Full name of the person who is writing this module.
   - `original_snakefile_author`: Full name of the person who wrote the Snakefile used to create this module. If this Snakefile doesn't exist, then this is simply the module author.
   - `module_run_per`: Possible values are "tumour" and "sample". This determines whether the module is meant to be run once per tumour (_e.g._ variant calling modules) or once per sample regardless of tissue status (_e.g._ BAM alignment and processing). Additional options will be added later, such as "tumour_cohort" and "sample_cohort" for [level-3 modules](README.md#module-levels).
   - `seq_type.genome`, `seq_type.capture`, and `seq_type.mrna`: Possible values are "paired", "unpaired", and "omit". This determines which sequencing data types (`seq_type`) are meant as input for this module (whole genome, hybrid capture-based, and RNA sequencing, respectively). Select "omit" if a `seq_type` is not meant to be analyzed by this module. If you selected "sample" for `module_run_per`, then you should use "unpaired" here. Lastly, if you selected "tumour" for `module_run_per`, you can select "paired" or "unpaired" depending on whether the module is meant to be run on tumour-normal pairs or not. By default, if you select "paired" and the user has tumours without matched normal samples, they will be expected to provide an unmatched normal sample to be used instead. This behaviour is determined by `run_unpaired_tumours_with` in `config/default.yaml` of the generated module.

6. Stay tuned! This document is a work in progress.
