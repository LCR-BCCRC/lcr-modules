.. _getting-started-for-users:

Getting Started for Users
=========================

.. Links (begin)

.. _demo snakefile: https://github.com/LCR-BCCRC/lcr-modules/blob/dev-bgrande/demo/snakefile

.. _miniconda: https://docs.conda.io/en/latest/miniconda.html

.. _anaconda: https://docs.anaconda.com/anaconda/install/

.. _lcr-modules: https://github.com/LCR-BCCRC/lcr-modules

.. _lcr-scripts: https://github.com/LCR-BCCRC/lcr-scripts

.. _test data: https://www.bcgsc.ca/downloads/lcr-modules/test_data/

.. _snakemake validation: https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html?highlight=schema#validation

.. _pandas data frame: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html

.. _pandas: https://pandas.pydata.org/docs/index.html

.. Links (end)

As a companion to these instructions, you can check out the `demo snakefile`_, where the placeholders contain actual values that work.

1. Install ``conda`` (Python version 3.6 or later) with miniconda_ or anaconda_. Ideally, create a project-specific conda environment. You can install Git using conda if need be.

   .. code:: bash

      conda install git

2. Clone the lcr-modules_ and lcr-scripts_ repositories.

   .. code:: bash

      git clone https://github.com/LCR-BCCRC/lcr-modules.git
      git clone https://github.com/LCR-BCCRC/lcr-scripts.git

3. Install the custom :py:mod:`oncopipe` python package included in the lcr-modules_ repository, which will also install dependencies such as ``snakemake`` and ``pandas``.

   .. code:: bash

      cd lcr-modules
      pip install -e oncopipe/

4. Test your environment with the `demo snakefile`_ using the following ``snakemake`` command. You can omit the ``--dry-run`` option if you first replace the placeholders in ``demo/data/`` with the actual `test data`_.

   .. code:: bash

      cd demo
      nice snakemake --dry-run --use-conda --cores <cores> _manta_all

5. Create a :ref:`samples-table` as a tab-delimited file with these :ref:`required-columns`. The following samples table is an example taken from the `demo snakefile`_ (``demo/data/samples.tsv``).

   +---------------------------+----------+-------------------+---------------+--------------+
   | sample_id                 | seq_type | patient_id        | tissue_status | genome_build |
   +===========================+==========+===================+===============+==============+
   | BLGSP-71-08-00508-10A-01D | capture  | BLGSP-71-08-00508 | normal        | hg38         |
   +---------------------------+----------+-------------------+---------------+--------------+
   | BLGSP-71-08-00508-01B-01D | capture  | BLGSP-71-08-00508 | tumour        | hg38         |
   +---------------------------+----------+-------------------+---------------+--------------+
   | BLGSP-71-08-00508-01B-01R | mrna     | BLGSP-71-08-00508 | tumour        | hg38         |
   +---------------------------+----------+-------------------+---------------+--------------+
   | BLGSP-71-06-00166-01B-01D | capture  | BLGSP-71-06-00166 | tumour        | hg38         |
   +---------------------------+----------+-------------------+---------------+--------------+
   | BLGSP-71-06-00166-01B-01R | mrna     | BLGSP-71-06-00166 | tumour        | hg38         |
   +---------------------------+----------+-------------------+---------------+--------------+

6. Add the following section to your :ref:`project-configuration` YAML file (*e.g.* ``config.yaml``) while updating values in angle brackets (``<...>``).

   .. code:: yaml

      lcr-modules:
         _shared:
            repository: "<path/to/lcr-modules>"
            lcr-scripts: "<path/to/lcr-scripts>"
            root_output_dir: "results/"
            scratch_directory: null

7. Load the :ref:`samples-table` and the :ref:`reference-files` workflow by adding the following lines near the beginning of your project snakefile.

   .. code:: python

      import oncopipe as op

      SAMPLES = op.load_samples("<path/to/samples.tsv>")
      config["lcr-modules"]["_shared"]["samples"] = SAMPLES

      subworkflow reference_files:
         workdir:
            "</path/to/reference_directory/>"
         snakefile:
            "<path/to/lcr-modules/workflows/reference_files/1.0/reference_files.smk>"
         configfile:
            "<path/to/lcr-modules/workflows/reference_files/1.0/config/default.yaml>"

8. Include and configure the modules you want to run by adding the following lines to your project snakefile. **Important:** Any values that need to be updated by the user will be indicated in the default module configuration by a ``# UPDATE`` comment. The following example updates the values specific to the ``manta`` and ``star`` modules. 

   .. code:: python

      # Load the default configuration file for each module
      configfile: "<path/to/lcr-modules/modules/manta/2.0/config/default.yaml>"
      configfile: "<path/to/lcr-modules/modules/star/1.0/config/default.yaml>"
      # ...

      # Load your project-specific configuration
      configfile: "<config.yaml>"

      # Update any configuration values (such as those indicated by `# UPDATE`)
      SAMPLE_BAM = "<data/{seq_type}_bams_{genome_build}/{sample_id}.bam>"
      config["lcr-modules"]["<manta>"]["inputs"]["<sample_bam>"] = SAMPLE_BAM
      config["lcr-modules"]["<manta>"]["inputs"]["<sample_bai>"] = SAMPLE_BAM + "bai"
      SAMPLE_FASTQ = "<data/{seq_type}_fastq_{genome_build}/{sample_id}>"
      config["lcr-modules"]["<star>"]["inputs"]["<sample_fastq_1>"] = SAMPLE_FASTQ + ".R1.fastq.gz"
      config["lcr-modules"]["<star>"]["inputs"]["<sample_fastq_2>"] = SAMPLE_FASTQ + ".R2.fastq.gz"
      config["lcr-modules"]["<star>"]["reference_params"]["<star_overhang>"] = "75"
      # ...

      # Include the snakefile for each module
      include: "<path/to/lcr-modules/modules/manta/2.0/manta.smk>"
      include: "<path/to/lcr-modules/modules/star/1.0/star.smk>"
      # ...

9. Launch snakemake by specifying the module target rule(s). See :ref:`snakemake-commands` for suggestions on how to run snakemake.

   .. code:: bash

      nice snakemake --dry-run --use-conda --cores <cores> _<module_name>_all

10. If you feel comfortable with the above steps, consider reading through the :ref:`advanced-usage`.

.. _samples-table:

Samples Table
=============

One of the requirements for using lcr-modules_ is a samples table. This format was selected for its flexibility. Each sample can be annotated with any amount of metadata, but for the purposes of lcr-modules_, there are only a few :ref:`required-columns`. These columns allow the modules to understand the relationship between the samples, especially for tumour-normal pairing. 

These requirements are encoded in schemas, which are stored and versioned in ``schemas/``. These schemas are used in conjunction with `snakemake validation`_. If the samples table doesn't confirm to a schemas that is required by a module, the user will given an informative error message. For example, the list of :ref:`required-columns` below is encoded in the ``base-1.0.yaml`` schema (located in ``schemas/base/``). The list of schemas will grow as modules are added with specific metadata requirements (*e.g.* strandedness of an RNA-seq library for expression quantification). 

The only format requirement for the samples table is that it is a `pandas data frame`_ (*i.e.* ``pandas.DataFrame``). Hence, the format of the file on disk doesn't matter. If you wish to use the :py:func:`oncopipe.load_samples` convenience function, note that it defaults to parsing tab-delimited files, but this can be overriden using the ``sep`` argument. The advantage of using :py:func:`oncopipe.load_samples` is that it offers a straightforward method for renaming your columns to comply with the schema(s). See :ref:`renaming-columns` for examples. 

Entity–relationship model
-------------------------

Before describing the required columns, it is useful to consider the entities related to each sample, namely ``patient``, ``biopsy``, ``sample``, ``library``, ``dataset``, and ``alignment``. These entities relate to one another in the following ways:

   **Relationships between entities:** Each patient has one or more biopsies (*e.g.* a tumour biopsy and a blood draw; tumour FF and FFPE biopsies). Each biopsy has one or more nucleic acid samples (*e.g.* DNA and RNA). Each sample has one or more sequencing libraries constructed from its nucleic acid samples (*e.g.* whole genome and RNA sequencing libraries for a tumour FF sample). Each sequenced library produces a a set of sequencing reads (*i.e.* a dataset) with one or more alignments (*e.g.* an hg19 and hg38 alignments), although there is generally a “canonical” alignment if more than one exists and thus a one-to-one relationship between datasets and alignments.

While the term "sample" generally refers to nucleic acid samples, lcr-modules_ uses the term to refer to the units of data that serve as input for the module, *i.e.* usually sequencing data in the form of FASTQ or BAM files. In most projects, there is a simple one-to-one relationship between these files and nucleic acid samples. In more complex projects where nucleic acid samples have more than one data file, the sample IDs will need to incorporate information to prevent duplicates.

.. _required-columns:

Required columns
----------------

Check out the :ref:`renaming-columns` section if your samples table has some of the required columns under different names. It also features a demonstration of the :py:func:`oncopipe.load_samples` convenience function you can use to load your TSV/CSV samples table. The :ref:`adding-columns` section is useful if you lack some of the required columns or can derive them from existing columns.

``seq_type`` – Sequencing data type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most common values for this column are ``genome`` (whole genome sequencing), ``mrna`` (RNA sequencing), ``capture`` (hybridization-capture or exome sequencing), and ``mirna`` (miRNA sequencing). While ``lcr-modules`` can handle any value for ``seq_type``, the modules are pre-configured for these common values. New values for ``seq_type`` will need to be added to the :ref:`pairing-configuration` of each module. If the pairing configuration would be same across multiple modules, it might be easier to set it under the ``_shared`` key in your :ref:`project-configuration`.

``sample_id`` – Sample identifiers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every ``seq_type`` and ``sample_id`` pair must be unique. In other words, if a tumour sample was sequenced using different technologies (*e.g.* whole genome and RNA sequencing), you can use the same sample ID since eachdata file will have a different ``seq_type`` (*e.g.* ``genome`` and ``mrna``, respectively). On the other hand, if you have been naming your samples based on patient ID and you have tumour-normal pairs, you will need to differentiate their sample IDs (*e.g.* with "T" and "N" suffixes). Similarly, if the same tumour sample has both FF and FFPE data files, you will also need to differentiate their sample IDs (*e.g.* with "FF" and "FFPE" suffixes). 

``tissue_status`` – Tumour or normal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This column classifies the samples as either ``tumour`` (or ``tumor``) and ``normal``. This information is required for tumour-normal paired analyses such as somatic variant calling. If you lack a matched normal samples, most modules support being run with an unmatched normal sample with the obvious caveats that the results will not be as clean. Check out :ref:`pairing-configuration` for more information on how to achieve this.

``patient_id`` – Patient identifiers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This column groups samples that originate from the same patient, *i.e.* that share the same underlying germline sequence. This information is primarily used in conjunction with the ``tissue_status`` column to generate all possible tumour-normal pairs from the list of samples.

``genome_build`` – Reference genome build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This column is only required if you have alignment (*i.e.* samples) using different genome builds. Otherwise, ``lcr-modules`` will assume that the single set of reference data (*e.g.* ``lcr-modules/references/hg38.yaml``) that you load is the one to use.

.. _renaming-columns:

Loading and renaming columns
----------------------------

For your convenience, the :py:func:`oncopipe.load_samples` function is provided to easily load your samples as a `pandas data frame`_. By default, the function assumes tab-delimited files, but you can change this using the ``sep`` argument. The function can also convert some of your columns to lowercase using the ``to_lowercase`` argument, which is useful to comply with some of the schemas. By default, it converts the ``tissue_status`` column to lowercase. It thus becomes trivial to load a samples table.

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv")

If your samples table uses different column names than those listed in :ref:`required-columns`, you can use the :py:func:`oncopipe.load_samples` function to rename your columns. For example, let's say you already have a samples table, but the sample ID and patient ID columns are named ``sample`` and ``patient`` rather than ``sample_id`` and ``patient_id``. You can easily achieve this as follows:

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv", sample_id = "sample", patient_id = "patient")

Alternatively, if the column names in your samples table differ systematically from the expected column names, you can rename them by passing a function to the ``renamer`` argument. You can also pass an anonymous ``lambda`` function. For instance, if you use two-letter prefixes with a period delimiter to indicate which entity a column describes (*e.g.* ``pt.`` for patient-related columns, ``lb.`` for library-related columns, etc.), you can remove the prefix from all columns using a regular expression with the following code:

.. code:: python

   import re
   import oncopipe as op
   remove_prefix = lambda x: re.sub(r"[a-z]{2}\.", "", x)
   SAMPLES = load_samples("samples.tsv", renamer=remove_prefix)

.. _adding-columns:

Adding and transforming columns
-------------------------------

If your samples table is missing a required column that has the same value for every sample (*e.g.* ``genome_build``), you can easily add the missing column in your snakefile using standard pandas_ syntax as follows:

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv")
   SAMPLES["genome_build"] = "hg38"

On the other hand, if your samples table is missing a required column that has different values for different samples, you can handle this one of two ways. If you can derive the missing column from existing columns, you can use standard pandas_ syntax to fill in the missing column. Otherwise, you can always resort to manually adding the missing column in the samples table on disk. The example below shows how the pandas_ syntax can be used to derive a ``tissue_status`` column by checking whether the ``sample_id`` column ends with the letter "T".

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv")
   SAMPLES["tissue_status"] = SAMPLES["sample_id"].str.endswith("T").map({True: "Tumour", False: "Normal"})

A similar approach can be taken if you have the columns, but they are formatted differently. For instance, if you encoded your sequencing data types as ``WGS`` and ``Exome`` instead of ``genome`` and ``capture``, respectively, you can use the ``map()`` method to switch to the expected values, as follows:

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv")
   SAMPLES["seq_type"] = SAMPLES["seq_type"].map({"WGS": "genome", "Exome": "capture"})

.. _reference-files:

Reference files
===============

The ``reference_files`` workflow is designed to simplify deployment of ``lcr-modules`` for any reference genome and on any computer system. This is achieved by (1) downloading the genome FASTA files and any additional reference files; (2) converting the additional files to match the same chromosome system as the genome builds (*e.g.* UCSC vs NCBI vs Ensembl); and (3) generate the required reference files from what was downloaded using snakemake rules. This approach also ensures that the steps taken to generate any reference file are tracked, ensuring their reproducibility.

More details will be added later.

.. _project-configuration:

Project Configuration
=====================



TODO: It is assumed that you have a project-specific configuration.

All configuration relating to ``lcr-modules`` is stored under the ``'lcr-modules'`` key in the snakemake ``config`` variable. The only exception to this rule is the reference data, which is stored under the ``'reference'`` key. The configuration for each module will be loaded under ``config['lcr-modules']``. For example, the ``manta`` configuration will be loaded to ``config['lcr-modules']['manta']``.

While most configuration is done at the module level, there are some values that are required at the project level. To avoid clashing with future module names, the project-level configuration is stored under the ``'_shared'`` key. (The underscore prefix stems from a Python convention.) It is worth noting that everything under ``'_shared'`` is set as the default value for each module unless that module has a specific value, which will overwrite the shared value.

You will need to specify a value for ``repository`` and ``root_output_dir``. If you have unpaired tumour samples, you will probably need to list the IDs for the samples to be used as unmatched normal samples in paired analyses under ``pairing_config``. See the example project configuration below for the required format.

-  ``repository``: File path for the cloned ``lcr-modules`` repository relative to your project snakefile. **This parameter is required.**
-  ``lcr-scripts``: File path for the cloned ``lcr-scripts`` repository relative to your project snakefile. **This parameter is required.**
-  ``root_output_dir``: Directory where all of the module output subdirectories will be created (*e.g.* ``results/manta-1.0/``). Technically, this shared parameter is optional and will default to ``'results/'``. I include it because I expect most users will want to customize this parameter.
- ``scratch_directory``: Directory where large temporary files can be written without worry of running out of space or clogging snapshots/backups. If set to ``null``, the files will be output locally.
-  ``pairing_config``: Optional unless you have unpaired tumours, in which case you need to specify which samples to use as unmatched normal samples for each ``seq_type`` under ``unmatched_normal_id``. See below for the required format.

Example project configuration
-----------------------------

.. code:: yaml

   lcr-modules:
       _shared:
           repository: "lcr-modules/"
           lcr-scripts: "lcr-scripts/"
           root_output_dir: "results/"
           pairing_config:
               genome:
                   unmatched_normal_id: "BLGSP-71-06-00286-99A-01D"
               capture:
                   unmatched_normal_id: "BLGSP-71-08-00508-10A-01D"

.. _pairing-configuration:

Pairing Configuration
=====================

SHARED

unmatched normal ID

Each module has a pairing configuration (*i.e.* ``pairing_config``). This configuration dictates what the module can handle in terms of paired and/or unpaired analyses for each sequencing data type (*i.e.* ``seq_type``). This information is used by the ``op.generate_runs_for_patient()`` function in ``oncopipe``.

Specifically, the following parameters are required for each ``seq_type``. The descriptions were taken from ``help(op.generate_runs_for_patient)``. An example pairing configuration can be found below.

-  ``run_paired_tumours``: ``True`` or ``False``, specifying whether to run paired tumours. Setting this to ``False`` is useful for naturally unpaired or tumour-only analyses (*e.g.* for RNA-seq).
-  ``run_unpaired_tumours_with``: ``None``, ``'unmatched_normal'``, or ``'no_normal'``, specifying what to pair with unpaired tumours. This cannot be set to ``None`` if ``run_paired_tumours_as_unpaired`` is ``True``. Provide value for ``unmatched_normal_id`` (see below) if this is set to ``'unmatched_normal'``.
-  ``unmatched_normal_id``: Identifier for the normal sample to be used with unpaired tumours when ``run_unpaired_tumours_with`` is set to ``'unmatched_normal'``. This is only required if you have unpaired samples, even if ``run_unpaired_tumours_with`` is set to ``'unmatched_normal'``. See `Project Configuration <#project-configuration>`__ for how to configure this parameter for your project.
-  ``run_paired_tumours_as_unpaired``: ``True`` or ``False``, specifying whether paired tumours should be run as unpaired (*i.e.* separate from their matched normal sample). This is useful for benchmarking purposes or preventing unwanted paired analyses (*e.g.* in RNA-seq analyses intended to be tumour-only).

Example pairing configuration
-----------------------------

This ``pairing_config`` was taken from the ``manta`` module. As you can see, the module can handle ``genome``, ``capture``, and ``mrna`` data. It treats ``genome`` and ``capture`` data the same way, namely by allowing unpaired tumours to be analyzed using unmatched normals (as opposed to a truly unpaired analysis without a normal sample). Also, paired tumours are not unnecessarily run as unpaired. In contrast, ``mrna`` data is run specifically in an unpaired fashion without a normal sample because tumour RNA-seq alignments generally do not have matched normal RNA-seq data. This can be overriden on a project-by-project basis.

.. code:: yaml

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

.. _snakemake-commands:

Snakemake Commands
==================

**Note:** Don’t forget to update any values in angle brackets (``<...>``).

snakemake profiles
------------------

The most convenient way of running snakemake is using `snakemake profiles <https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles>`__. Each profile contains a YAML file that dictates the default command-line options to use. This way, you don’t have to remember all those snakemake options.

GSC snakemake profiles
~~~~~~~~~~~~~~~~~~~~~~

Make sure you first install the custom GSC snakemake profiles using `these instructions <https://github.com/LCR-BCCRC/snakemake-profiles#installation>`__. Then, you can use each profile using `these commands <https://github.com/LCR-BCCRC/snakemake-profiles#usage>`__.

Explicit commands
-----------------

If you prefer to spell out all of the command-line options in your snakemake commands, example commands are included below. These may eventually become out of sync with the above snakemake profiles. Feel free to compare with the list of arguments for `local usage <https://github.com/LCR-BCCRC/snakemake-profiles/blob/master/gphosts/config.yaml>`__ or `cluster usage <https://github.com/LCR-BCCRC/snakemake-profiles/blob/master/numbers/config.yaml>`__.

Local usage
~~~~~~~~~~~

.. code:: bash

   # See below for determining <cores>
   nice snakemake --printshellcmds --use-conda --cores <cores> <targets>

Cluster usage
~~~~~~~~~~~~~

.. code:: bash

   nice snakemake --cluster-sync "srun --partition=all --ntasks=1 --nodes=1 --output=none --error=none --job-name={rule} --cpus-per-task={threads} --mem={resources.mem_mb}" --max-jobs-per-second=5 --max-status-checks-per-second=10 --local-cores=1 --latency-wait=120 --jobs=1000 --default-resources="mem_mb=2000" --printshellcmds --use-conda <targets>

Extra information
-----------------

Determining value for ``--cores``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To determine the number of cores to grant to snakemake, compare the number of installed cores and the current load on the server. These values can either be obtained precisely using the commands below, or they can be estimated by looking at the output of the ``htop`` `command <https://hisham.hm/htop/index.php?page=screenshots>`__. I generally select a value for ``--cores`` equal to the number of installed cores minus the server load minus 10-20 to leave some buffer.

.. code:: bash

   # Get the number of installed logical cores
   nproc
   # Get the average server load over the past 5 minutes
   cut -d " " -f 2 /proc/loadavg

Increasing ``ulimit``
~~~~~~~~~~~~~~~~~~~~~

snakemake tends to spawn A LOT of processes and open A LOT of files depending on the number of running and pending jobs. You may eventually start running into cryptic errors about processors not being able to start or files not being able to be opened. This happens when you run into user limits. You can get around this issue by increasing the user limits with the ``ulimit`` command. However, there are hard limits set by administrators that determine the maximum permitted for non-admin users. You can always ask your administrators to increase these hard limits for certain machines to run snakemake.

GSC ``ulimit`` setup
^^^^^^^^^^^^^^^^^^^^

GSC users can include the following code in their ``.bashrc`` file to increase their ulimits based on the server. Notice how the ``n104`` numbers head node has a much higher hard limit than the other head nodes. This is because it was manually increased when ``n104`` was the only head node. For this reason, it is recommended that GSC users specically log into ``n104`` instead of ``numbers``, which will assign you to a random head node.

.. code:: bash

   # Only change these values for interactive shells
   if [[ $- == *i* ]]; then
     if [[ "$HOSTNAME" == "n104" ]]; then
       # Change the max number of processes
       ulimit -u 32768
       # Change the max number of file descriptors
       ulimit -n 100000
     fi
   fi

Creating ``nice`` processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You will notice that the ``snakemake`` commands below are all prepended with ``nice``. Briefly, this has the effect of lowering the priority of your snakemake process. Now, you’re probably wondering why would you ever want to do that. Granted, compute resources should be utilized on a first come, first served basis, but in practice, not every user will pay close attention to who is already running jobs on a server.

Ultimately, it doesn’t matter whether this act is intentional, an accident, or due to insufficient knowledge of how to manage shared compute resources. If someone launches a job that uses more cores than are available, your snakemake process will be competing for CPU time, and this will make both processes take longer to complete.

In this situation, we should fall back on the motto from the wise Michelle Obama: “When they go low, we go high.” In this case, we follow this rule quite literally, because the ``nice`` command will increase the “niceness” value of your snakemake process, which will cede CPU time to competing processes with lower (usually default) “niceness” values until they’re done.

Submitting cluster jobs remotely
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to submit jobs to a cluster remotely via SSH. This could be useful in situations where you have quick jobs that you don’t want to submit to the cluster, but you also don’t want to run locally on the cluster head node. **Important:** This section assumes that you have SSH keys set up, allowing SSH login to the head node without entering a password.

The command below differs from the explicit command above simply by prepending the ``srun`` command in ``--cluster-sync`` with ``ssh <head_node>``, where ``<head_node>`` is the cluster head node where you run ``srun`` normally. You can now increase the value for ``--local-cores`` (see above for how to determine this value).

.. code:: bash

   nice snakemake --local-cores=<cores> --cluster-sync "ssh <head_node> srun --partition=all --ntasks=1 --nodes=1 --output=none --error=none --job-name={rule} --cpus-per-task={threads} --mem={resources.mem_mb}" --max-jobs-per-second=5 --max-status-checks-per-second=10 --latency-wait=120 --jobs=1000 --default-resources="mem_mb=2000" --printshellcmds --use-conda <targets>

.. _advanced-usage:

Advanced Usage
==============

Directory placeholders
----------------------

When specifying any value in the module configuration, you can use the following shorthands as placeholders in the string. They will be replaced with the actual values dynamically. See the `Parameterization <#parameterization>`__ section below for example usage.

-  ``{REPODIR}``: The ``lcr-modules`` repository directory. This corresponds to the ``repository`` value under ``_shared`` in the ``lcr-modules`` configuration.
-  ``{MODSDIR}``: The current module subdirectory. This corresponds to ``{REPODIR}/modules/<name>/<version>``.
-  ``{SCRIPTSDIR}``: The ``lcr-scripts`` repository directory. This corresponds to the ``lcr-scripts`` value under ``_shared`` in the ``lcr-modules`` configuration.

Convenience set functions
-------------------------

The `Setup Instructions <#setup-instructions>`__ demonstrate that everything is configured using the same snakemake ``config`` nested dictionary object, generally under the ``'lcr-modules'`` key. While transparent, it results in verbose code, such as:

.. code:: python

   config["lcr-modules"]["manta"]["inputs"]["sample_bam"] = SAMPLE_BAM

Alternatively, you can use the so-called convenience “set functions” to simplify the code somewhat. In order to use them, you must first enable them. Behind the scenes, the snakemake ``config`` object is stored internally for easy access.

.. code:: python

   op.enable_set_functions(config)

The first set function you can use is :py:func:`oncopipe.set_samples()`, which sets the samples you want to use at the shared or module level. The first argument corresponds to the module name (or ``"_shared"``), and all subsequent arguments should be sample tables each formatted as a `pandas data frame`_. This function automatically concatenates the data frames that are provided. Here, ``SAMPLES`` is the complete samples table, whereas ``GENOMES`` and ``CAPTURES`` are sample subsets generated from ``SAMPLES`` using :py:func:`oncopipe.filter_samples()`.

.. code:: python

   op.set_samples("_shared", SAMPLES)
   op.set_samples("_shared", GENOMES, CAPTURES)

The second function you can use is :py:func:`oncopipe.set_input()`, which sets the given input for a module. Just like ``op.set_samples()``, the first argument is the module name, but this function should not be used for ``"shared"``. The second argument is the name of the input file as listed in the module’s configuration file. Lastly, the third argument is the value you wish to provide for that input file, which generally is a string value containing the available wildcards (once again, as listed in the module’s configuration file). That said, you could provide a conditional value as described below in `Parameterization <#parameterization>`__.

.. code:: python

   op.set_input("manta", "sample_bam", SAMPLE_BAM)

Parameterization
----------------

Sometimes, a parameter or input file depends on some sample attribute. This sample attribute can be stored in the file as a wildcard or in the samples tables as a column. Two functions are available to parameterize virtually anything, namely ``op.switch_on_wildcard()`` and ``op.switch_on_column()``. These functions are useful for both module users and module developers. Read their documentation for more details, *e.g.* ``help(op.switch_on_wildcard)``.

In the example below, I want to override the default Manta configuration and provide the high-sensitivity version for ``mrna`` and ``capture`` tumour samples. This piece of code would be added after loading the module configuration but before including the module snakefile.

.. code:: python

   MANTA_CONFIG_OPTIONS = {
       "_default": "{MODSDIR}/etc/manta_config.default.ini",
       "mrna": "{MODSDIR}/etc/manta_config.high_sensitivity.ini",
       "capture": "{MODSDIR}/etc/manta_config.high_sensitivity.ini",
   }
   MANTA_CONFIG_SWITCH = op.switch_on_wildcard("seq_type", MANTA_CONFIG_OPTIONS)
   op.set_input("manta", "manta_config", MANTA_CONFIG_SWITCH)

.. _faq:

Frequently Asked Questions
==========================

How do I handle a conda environment that fails to build?
--------------------------------------------------------

While conda brings us much closer to computational reproducibility, it isn’t perfect. Issues arise when conda packages are removed from `Anaconda Cloud <https://anaconda.org/>`__ or when the dependency resolution algorithm changes. We suggest you try the following steps in order:

1. Remove the build IDs from the conda environment YAML file, although this should already be the case for all environments in ``lcr-modules``.
2. Remove the versions for the offending package(s) (*i.e.* the one(s) mentioned in the error message).
3. Remove the offending packages altogether.
4. Remove the dependency packages, leaving only the “target packages”. This generally means subsetting to the core conda packages listed in a module’s README for the environment in question. While extreme, the hope is that the versions of the dependency packages are not crucial for maintaining scientific reproducibility.
5. Remove the versions for the target packages.
6. If you reach this point, it usually means that a target package is problematic. If possible, replace that package with the same (or similar) version from another Anaconda channel. Ideally, restore the YAML file first and cycle through the previous steps.
7. Install the software tools manually (ideally the versions specified in the YAML file) and ensure they are available in your ``PATH`` environment variable.

What is up with the underscore prefix (*e.g.* in rule names)?
-------------------------------------------------------------

This underscore prefix stems from a Python convention. In ``lcr-modules``, it is generally meant to avoid name conflits. For example, in the ``manta`` module, the final target rule is called ``_manta_all`` just in case the user already has a rule called ``manta_all``. While this is unlikely, as modules are loaded, the risk for a conflict increases. Hence, the underscore prefix is a precautionary measure.
