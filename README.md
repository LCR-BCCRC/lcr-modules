# LCR Snakemake Analysis Modules

**Note:** This repository includes a custom Python package called `modutils`. This README assumes that it is loaded into Python using `import modutils as md`. Hence, any functions therein are referred to with the `md.` prefix.

## Setup Instructions

These instructions assume your working directory is your project root. They also assume that you have conda set up and your current environment has Python 3.6 or later.

1. Clone the `lcr-modules` repository **recursively**.

```bash
git clone --recursive https://github.com/LCR-BCCRC/lcr-modules.git
```

2. Install `snakemake` (5.4 or later), `pandas`, and the custom `modutils` Python packages into your conda environment.

```bash
# `snakemake-minimal` lacks extraneous dependencies and is faster to install
conda install snakemake-minimal>=5.4 pandas
pip install -e lcr-modules/modutils
```

3. (Optional) Test your environment with the demo project in this repository. You shouldn't get any error after running the following `snakemake` command:

```bash
cd lcr-modules/demo
snakemake -n _manta_all
cd ../..
```

4. Create a samples table as a tabular file with the following columns: `seq_type`, `sample_id`, `patient_id`, `tissue_status`, and optionally, `genome_build`. **Important:** You must follow the file format described [below](#samples-table).

5. Create a project configuration YAML file (if you don't already have one), add the following section, and load it using `configfile:` in your project Snakefile. See [Project Configuration](#project-configuration) for more details.

```yaml
lcr-modules:
  _shared:
    repository: "lcr-modules/"
    root_output_dir: "results/"
```

6. Add the following lines **once** near the beginning of your project Snakefile, updating any values in angle brackets (`<...>`).

```python
import modutils as md
configfile: "lcr-modules/references/<hg38>.yaml"
SAMPLES = md.load_samples("<path/to/samples.tsv>")
config["lcr-modules"]["_shared"]["samples"] = SAMPLES
```

7. Add the following lines to your project Snakefile **for each module**, updating any values in angle brackets (`<...>`). This specific order is required. **Important:** Read each module README for any module-specific configuration.

```python
configfile: "lcr-modules/modules/<manta/1.0>/config/default.yaml"
config["lcr-modules"]["<manta>"]["inputs"]["<sample_bam>"] = <SAMPLE_BAM>
# Repeat previous line for any other required input files
include: "lcr-modules/modules/<manta/1.0>/manta.smk"
```

8. Launch Snakemake for the target rule of any module you added. See [Snakemake Commands](#snakemake-commands) for suggestions on how to run Snakemake.

9. ???

10. Profit! And reproducible research!

If you feel comfortable with the above steps, consider reading through the suggestions laid out in [Advanced Usage](#advanced-usage).

## Samples Table

This section is a human-friendly summary of what is described in the base schema (see `lcr-modules/schemas/base.yaml`). In case of any discrepancies, the schema file takes precedence. The samples table can be stored on disk using any column delimiter, but `modutils` expects tab-delimited files by default.

### Entity–relationship model

Before describing the required columns, it is useful to consider the entities related to each sample, namely `patient`, `biopsy`, `sample`, `library`, `dataset`, and `alignment`. The relationships between each entity are spelled out in the blockquote below. While the term "sample" can easily refer to any of these entities except for `patient`, we use it to indicate the entities that will be analyzed, which usually are datasets/alignments.

> Each patient has one or more biopsies (_e.g._ a tumour biopsy and a blood draw; tumour FF and FFPE biopsies). Each biopsy has one or more nucleic acid samples (_e.g._ DNA and RNA). Each sample has one or more sequencing libraries constructed from its nucleic acid samples (_e.g._ whole genome and RNA sequencing libraries for a tumour FF sample). Each sequenced library produces a a set of sequencing reads (_i.e._ a dataset) with one or more alignments (_e.g._ an hg19 and hg38 alignments), although there is generally a "canonical" alignment if more than one exists and thus a one-to-one relationship between datasets and alignments.

### Required columns

- **`seq_type`:** Sequencing data type. The most common values for this column are `genome` for whole genome sequencing, `mrna` for RNA sequencing, `capture` for any form of hybridization-capture sequencing (including exome sequencing), and `mirna` for miRNA sequencing. While `lcr-modules` was designed to handle any value for `seq_type`, the modules are currently configured for these common values. New values for `seq_type` will need to be added to the [pairing configuration](#pairing-configuration) of each module.
- **`sample_id`:** Sample identifiers. These must be unique within each value for `seq_type`. For example, your sample IDs can be the patient ID suffixed with "T" and "N" for tumour and normal samples, and the genome and RNA sequencing alignments for a given sample can use the same ID. However, if you have more than one sample for a given tumour (_e.g._ FF and FFPE), their IDs will need to differ.
- **`patient_id`:** Patient identifiers. In essence, this column groups all samples that share the same germline sequence. This information is used in conjunction with the `tissue_status` column to generate all possible tumour-normal pairs for paired analyses.
- **`tissue_status`:** Tumour or normal. This column is used in conjunction with the `patient_id` column to generate all possible tumour-normal pairs for paired analyses. It is also used to generate the list of tumours for tumour-focused analyses.
- **`genome_build`:** Reference genome build. This column is only required if you have alignment (_i.e._ samples) using different genome builds. Otherwise, `lcr-modules` will assume that the single set of reference data (_e.g._ `lcr-modules/references/hg38.yaml`) that you load is the one to use.

### Renaming columns

If you already have a file with this information but using different column names, you can use the `load_samples()` function in `modutils` to rename your columns. For example, if you use `sample` and `patient` as the column headers for your sample and patient IDs, you can rename them as follows:

```python
md.load_samples("samples.tsv", sample_id="sample", patient_id="patient")
```

Alternatively, you can rename the column using a function with the `renamer` argument. For instance, if you use two-letter prefixes to indicate which entity a column describes (_e.g._ `pt.` for patient-related columns, `lb.` for library-related columns, etc.), you can remove the prefix from all columns using a regular expression with the following code:

```python
import re
remove_prefix = lambda x: re.sub(r"[a-z]{2}\.", "", x)
SAMPLES = load_samples("samples.tsv", renamer=remove_prefix)
```

## Project Configuration

All configuration relating to `lcr-modules` is stored under the `'lcr-modules'` key in the Snakemake `config` variable. The only exception to this rule is the reference data, which is stored under the `'reference'` key. The configuration for each module will be loaded under `config['lcr-modules']`. For example, the `manta` configuration will be loaded to `config['lcr-modules']['manta']`.

While most configuration is done at the module level, there are some values that are required at the project level. To avoid clashing with future module names, the project-level configuration is stored under the `'_shared'` key. (The underscore prefix stems from a Python convention.) It is worth noting that everything under `'_shared'` is set as the default value for each module unless that module has a specific value, which will overwrite the shared value.

You will need to specify a value for `repository` and `root_output_dir`. If you have unpaired tumour samples, you will probably need to list the IDs for the samples to be used as unmatched normal samples in paired analyses under `pairing_config`. See the example project configuration below for the required format.

- **`repository`:** File path for the cloned `lcr-modules` repository relative to your project Snakefile. **This parameter is required.**
- **`root_output_dir`:** Directory where all of the module output subdirectories will be created (_e.g._ `results/manta-1.0/`). Technically, this shared parameter is optional and will default to `'results/'`. I include it because I expect most users will want to customize this parameter.
- **`pairing_config`:** Optional unless you have unpaired tumours, in which case you need to specify which samples to use as unmatched normal samples for each `seq_type` under `unmatched_normal_id`. See below for the required format.

### Example project configuration

```yaml
lcr-modules:
  _shared:
    repository: "lcr-modules/"
    root_output_dir: "results/"
    pairing_config:
      genome:
        unmatched_normal_id: "BLGSP-71-06-00286-99A-01D"
      capture:
        unmatched_normal_id: "BLGSP-71-08-00508-10A-01D"
```

## Pairing Configuration

Each module has a pairing configuration (_i.e._ `pairing_config`). This configuration dictates what the module can handle in terms of paired and/or unpaired analyses for each sequencing data type (_i.e._ `seq_type`). This information is used by the `md.generate_runs_for_patient()` function in `modutils`.

Specifically, the following parameters are required for each `seq_type`. The descriptions were taken from `help(md.generate_runs_for_patient)`. An example pairing configuration can be found below.

- **`run_paired_tumours`:** `True` or `False`, specifying whether to run paired tumours. Setting this to `False` is useful for naturally unpaired or tumour-only analyses (_e.g._ for RNA-seq).
- **`run_unpaired_tumours_with`:** `None`, `'unmatched_normal'`, or `'no_normal'`, specifying what to pair with unpaired tumours. This cannot be set to `None` if `run_paired_tumours_as_unpaired` is `True`. Provide value for `unmatched_normal_id` (see below) if this is set to `'unmatched_normal'`.
- **`unmatched_normal_id`:** Identifier for the normal sample to be used with unpaired tumours when `run_unpaired_tumours_with` is set to `'unmatched_normal'`. This is only required if you have unpaired samples, even if `run_unpaired_tumours_with` is set to `'unmatched_normal'`. See [Project Configuration](#project-configuration) for how to configure this parameter for your project.
- **`run_paired_tumours_as_unpaired`:** `True` or `False`, specifying whether paired tumours should be run as unpaired (_i.e._ separate from their matched normal sample). This is useful for benchmarking purposes or preventing unwanted paired analyses (_e.g._ in RNA-seq analyses intended to be tumour-only).

### Example pairing configuration

This `pairing_config` was taken from the `manta` module. As you can see, the module can handle `genome`, `capture`, and `mrna` data. It treats `genome` and `capture` data the same way, namely by allowing unpaired tumours to be analyzed using unmatched normals (as opposed to a truly unpaired analysis without a normal sample). Also, paired tumours are not unnecessarily run as unpaired. In contrast, `mrna` data is run specifically in an unpaired fashion without a normal sample because tumour RNA-seq alignments generally do not have matched normal RNA-seq data. This can be overriden on a project-by-project basis.

```yaml
# Taken from lcr-modules/modules/manta/1.0/config/default.yaml
pairing_config:
  genome:
    run_unpaired_tumours_with: "unmatched_normal"
    run_paired_tumours: True
    run_paired_tumours_as_unpaired: False
  capture:
    run_unpaired_tumours_with: "unmatched_normal"
    run_paired_tumours: True
    run_paired_tumours_as_unpaired: False
  mrna:
    run_unpaired_tumours_with: "no_normal"
    run_paired_tumours: False
    run_paired_tumours_as_unpaired: True
```

## Snakemake Commands

**Note:** Don't forget to update any values in angle brackets (`<...>`).

### Snakemake profiles

The most convenient way of running Snakemake is using [Snakemake profiles](https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles). Each profile contains a YAML file that dictates the default command-line options to use. This way, you don't have to remember all those Snakemake options.

#### GSC Snakemake profiles

Make sure you first install the custom [GSC Snakemake profiles](https://github.com/LCR-BCCRC/snakemake-profiles.git) using [these instructions](https://github.com/LCR-BCCRC/snakemake-profiles#installation). Then, you can use each profile using the commands below.

```bash
# Run this command on numbers (ideally, the n104 node)
nice snakemake --profile numbers --keep-going <targets>
# Run this command on the gphosts (see below for determining <cores>)
nice snakemake --profile gphosts --keep-going --cores <cores> <targets>
```

### Explicit commands

If you prefer to spell out all of the command-line options in your Snakemake commands, example commands are included below. These may eventually become out of sync with the above Snakemake profiles. Feel free to compare with the list of arguments for [local usage](https://github.com/LCR-BCCRC/snakemake-profiles/blob/master/gphosts/config.yaml) or [cluster usage](https://github.com/LCR-BCCRC/snakemake-profiles/blob/master/numbers/config.yaml).

#### Local usage

```bash
# See below for determining <cores>
nice snakemake --printshellcmds --use-conda --keep-going --cores <cores> <targets>
```

#### Cluster usage

```bash
nice snakemake --cluster-sync "srun --partition=all --ntasks=1 --nodes=1 --output=none --error=none --job-name={rule} --cpus-per-task={threads} --mem={resources.mem_mb}" --max-jobs-per-second=5 --max-status-checks-per-second=10 --local-cores=1 --latency-wait=120 --jobs=1000 --default-resources="mem_mb=2000" --printshellcmds --use-conda --keep-going <targets>
```

### Extra information

#### Determining value for `--cores`

To determine the number of cores to grant to Snakemake, compare the number of installed cores and the current load on the server. These values can either be obtained precisely using the commands below, or they can be estimated by looking at the output of the [`htop` command](https://hisham.hm/htop/index.php?page=screenshots). I generally select a value for `--cores` equal to the number of installed cores minus the server load minus 10-20 to leave some buffer.

```bash
# Get the number of installed logical cores
nproc
# Get the average server load over the past 5 minutes
cut -d " " -f 2 /proc/loadavg
```

#### Increasing `ulimit`

Snakemake tends to spawn A LOT of processes and open A LOT of files depending on the number of running and pending jobs. You may eventually start running into cryptic errors about processors not being able to start or files not being able to be opened. This happens when you run into user limits. You can get around this issue by increasing the user limits with the `ulimit` command. However, there are hard limits set by administrators that determine the maximum permitted for non-admin users. You can always ask your administrators to increase these hard limits for certain machines to run Snakemake.

##### GSC `ulimit` setup

GSC users can include the following code in their `.bashrc` file to increase their ulimits based on the server. Notice how the `n104` numbers head node has a much higher hard limit than the other head nodes. This is because it was manually increased when `n104` was the only head node. For this reason, it is recommended that GSC users specically log into `n104` instead of `numbers`, which will assign you to a random head node.

```bash
# Only change these values for interactive shells
if [[ $- == *i* ]]; then
  if [[ "$HOSTNAME" == "n104" ]]; then
    # Change the max number of processes
    ulimit -u 32768
    # Change the max number of file descriptors
    ulimit -n 100000
  fi
fi
```

#### Creating `nice` processes

You will notice that the `snakemake` commands below are all prepended with `nice`. Briefly, this has the effect of lowering the priority of your Snakemake process. Now, you're probably wondering why would you ever want to do that. Granted, compute resources should be utilized on a first come, first served basis, but in practice, not every user will pay close attention to who is already running jobs on a server.

Ultimately, it doesn't matter whether this act is intentional, an accident, or due to insufficient knowledge of how to manage shared compute resources. If someone launches a job that uses more cores than are available, your Snakemake process will be competing for CPU time, and this will make both processes take longer to complete.

In this situation, we should fall back on the motto from the wise Michelle Obama: "When they go low, we go high." In this case, we follow this rule quite literally, because the `nice` command will increase the "niceness" value of your Snakemake process, which will cede CPU time to competing processes with lower (usually default) "niceness" values until they're done.

#### Submitting cluster jobs remotely

It is possible to submit jobs to a cluster remotely via SSH. This could be useful in situations where you have quick jobs that you don't want to submit to the cluster, but you also don't want to run locally on the cluster head node. **Important:** This section assumes that you have SSH keys set up, allowing SSH login to the head node without entering a password.

The command below differs from the explicit command above simply by prepending the `srun` command in `--cluster-sync` with `ssh <head_node>`, where `<head_node>` is the cluster head node where you run `srun` normally. You can now increase the value for `--local-cores` (see above for how to determine this value).

```bash
nice snakemake --local-cores=<cores> --cluster-sync "ssh <head_node> srun --partition=all --ntasks=1 --nodes=1 --output=none --error=none --job-name={rule} --cpus-per-task={threads} --mem={resources.mem_mb}" --max-jobs-per-second=5 --max-status-checks-per-second=10 --latency-wait=120 --jobs=1000 --default-resources="mem_mb=2000" --printshellcmds --use-conda --keep-going <targets>
```

## Advanced Usage

### Directory shorthands

When specifying any value in the module configuration, you can use the following shorthands as placeholders in the string. They will be replaced with the actual values dynamically. See the [Parameterization](#parameterization) section below for example usage.

- **`REPODIR`:** `lcr-modules` repository directory. This is equivalent to the value stored in `config['_shared']['repository']` (see [Project Configuration](#project-configuration)).
- **`MODSDIR`:** `lcr-modules` module subdirectory. This corresponds to the module subdirectory within `REPODIR`, _e.g._ `lcr-modules/modules/manta/1.0/`.

### Convenience set functions

The [Setup Instructions](#setup-instructions) demonstrate that everything is configured using the same Snakemake `config` nested dictionary object, generally under the `'lcr-modules'` key. While transparent, it results in verbose code, such as:

```python
config["lcr-modules"]["manta"]["inputs"]["sample_bam"] = SAMPLE_BAM
```

Alternatively, you can use the so-called convenience "set functions" to simplify the code somewhat. In order to use them, you must first enable them. Behind the scenes, the Snakemake `config` object is stored internally for easy access.

```python
md.enable_set_functions(config)
```

The first set function you can use is `md.set_samples()`, which sets the samples you want to use at the shared or module level. This function automatically concatenates the data frames that are provided.

```python
md.set_samples("_shared", SAMPLES)
md.set_samples("_shared", GENOMES, CAPTURES)
```

The second function you can use is `md.set_input()`, which set the given input for a module.

```python
md.set_input("manta", "sample_bam", SAMPLE_BAM)
```

While also possible with the more verbose approach,

### Parameterization

Sometimes, a parameter or input file depends on some sample attribute. This sample attribute can be stored in the file as a wildcard or in the samples tables as a column. Two functions are available to parameterize virtually anything, namely `md.switch_on_wildcard()` and `md.switch_on_column()`. These functions are useful for both module users and module developers. Read their documentation for more details, _e.g._ `help(md.switch_on_wildcard)`.

In the example below, I want to override the default Manta configuration and provide the high-sensitivity version for `mrna` and `capture` tumour samples. This piece of code would be added after loading the module configuration but before including the module Snakefile.

```python
MANTA_CONFIG_OPTIONS = {
    "_default": "{MODSDIR}/etc/manta_config.default.ini",
    "mrna": "{MODSDIR}/etc/manta_config.high_sensitivity.ini",
    "capture": "{MODSDIR}/etc/manta_config.high_sensitivity.ini",
}
MANTA_CONFIG_SWITCH = md.switch_on_wildcard("seq_type", MANTA_CONFIG_OPTIONS)
md.set_input("manta", "manta_config", MANTA_CONFIG_SWITCH)
```
