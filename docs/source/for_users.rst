.. include:: links.rst

.. _getting-started-user:

Getting Started
===============

**Important:** Be sure to update any values in angle brackets (``<...>``). File paths can either be absolute or relative to the working directory (usually where you run your ``snakemake`` command). 

1. Install ``conda`` (Python version 3.6 or later) with Miniconda_ or Anaconda_. Ideally, create a project-specific conda environment. You can install Git using conda if need be.

   .. code:: bash

      # Optional: Create a project-specific conda environment
      # conda create -n <project_name> "python>=3.6"
      # conda activate <project_name>
      
      conda install git

2. Clone the `lcr-modules repository`_ and the `lcr-scripts repository`_.

   .. code:: bash

      git clone https://github.com/LCR-BCCRC/lcr-modules.git
      git clone https://github.com/LCR-BCCRC/lcr-scripts.git

3. Install the custom :py:mod:`oncopipe` python package included in the `lcr-modules repository`_, which will also install dependencies such as ``snakemake`` and ``pandas``.

   .. code:: bash

      cd lcr-modules/
      pip install -e oncopipe/

4. Test your environment with the `Demo Snakefile`_ using the following ``snakemake`` command. If you want to actually run the test, check out :ref:`running-the-demo-project`.

   .. code:: bash

      cd demo
      nice snakemake --dry-run --use-conda all

5. Create a :ref:`sample-table` as a tab-delimited file with these :ref:`required-columns`. The following sample table is an example taken from the `Demo Project`_ (``demo/data/samples.tsv``).

   +---------------+----------+------------+---------------+--------------+
   | sample_id     | seq_type | patient_id | tissue_status | genome_build |
   +===============+==========+============+===============+==============+
   | TCRBOA7-N-WEX | capture  | TCRBOA7    | normal        | grch37       |
   +---------------+----------+------------+---------------+--------------+
   | TCRBOA7-T-WEX | capture  | TCRBOA7    | tumour        | grch37       |
   +---------------+----------+------------+---------------+--------------+
   | TCRBOA7-T-RNA | mrna     | TCRBOA7    | tumour        | grch37       |
   +---------------+----------+------------+---------------+--------------+

6. Add the :ref:`lcr-modules-configuration` to your project configuration YAML file (*e.g.* ``config.yaml``). If you have unpaired tumour samples, you will need to add a ``pairing_config``. Check out :ref:`handling-unpaired-tumours` section for more information and :ref:`within-the-configuration-file` for an example of how this is done.

   .. code:: yaml

      # Update `scratch_directory` if you have a space to store large intermediate files
      lcr-modules:
         _shared:
            repository: "<path/to/lcr-modules>"
            lcr-scripts: "<path/to/lcr-scripts>"
            root_output_dir: "results/"
            scratch_directory: null

7. Load the :ref:`sample-table`, :ref:`specify-the-input-samples`, and add the :ref:`reference-files-workflow` by including the following lines in your project snakefile anywhere before the module snakefiles (see next step).

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
   
   .. gap

   .. note::

      | **BCGSC Users**
      | You can set the ``reference_files`` working directory (``workdir``) to the following file path in order to benefit from pregenerated reference files:
      
      .. code::

         /projects/bgrande/reference_files

8. Include and configure the modules you want to run by adding the following lines to your project snakefile. 

   **Important** 
   
   - Any values that need to be updated by the user will be indicated in the default :ref:`module-configuration` by ``UPDATE`` comments. 

   - We recommend following the order shown below: (1) load the default module configuration files; (2) load the project-specific configuration file; and (3) include the module snakefiles. 
   
   - The following example assumes that these values are updated in the project-specific configuration file. For more information, check out :ref:`updating-configuration-values`.

   .. code:: python

      # Load the default configuration file for each module
      configfile: "<path/to/lcr-modules/modules/manta/2.0/config/default.yaml>"
      configfile: "<path/to/lcr-modules/modules/star/1.0/config/default.yaml>"
      # ...

      # Load your project-specific configuration
      configfile: "<config.yaml>"

      # Include the snakefile for each module
      include: "<path/to/lcr-modules/modules/manta/2.0/manta.smk>"
      include: "<path/to/lcr-modules/modules/star/1.0/star.smk>"
      # ...

9. Launch snakemake by specifying the module target rule(s). See :ref:`snakemake-commands` for suggestions on how to run snakemake.

   .. code:: bash

      nice snakemake --dry-run --use-conda --cores <cores> _<star>_all _<manta>_all

10. If you feel comfortable with the above steps, consider reading through the :ref:`advanced-usage`. For example, you can use :ref:`conditional-module-behaviour-user` to set different input file paths for different sequencing data types (*e.g.* ``genome`` and ``mrna``).

.. _running-the-demo-project:

Running the Demo Project
========================

Before running the `Demo Project`_, you will need to download the `Test Data`_, which is composed of exome and RNA sequencing data for a follicular lymphoma case. Specifically, the tumour sample has both exome and RNA data, whereas the normal sample only has exome data. We acknowledge the Texas Cancer Research Biobank, Baylor College of Medicine Human Genome Sequencing Center, and the study participants for graciously providing these open-access data. These data are described in more detail in the following paper. **Important:** Before downloading the data, you must first read and agree to the terms laid out in the `Test Data README`_, which is a requirement for redistributing the data. 

   Becnel, L. B. et al. An open access pilot freely sharing cancer genomic data from participants in Texas. Sci. Data 3:160010 doi: 10.1038/sdata.2016.10 (2016).

.. note::

   | **BCGSC Users**
   | You can replace the empty placeholder files with symbolic links to the local copies of the test data using the following command:
   
   .. code:: bash

      ln -sf /projects/bgrande/lcr-modules/test_data/*.{bam,bai,fastq.gz} ./demo/data/

Once the `Test Data`_ is downloaded, you will need to update the placeholders in ``demo/data/`` with the downloaded files (or symbolic links to the files). At that point, you can technically run the `Demo Snakefile`_ by omitting the ``--dry-run`` option from the command in the :ref:`getting-started-user` instructions, but you might want to update the value set in ``demo/config.yaml`` under ``scratch_directory`` to an space where you can readily store large intermediate files (*e.g.* a directory without snapshots or backups).

If you are interested in learning how you can conditionally use the STAR BAM files for RNA-seq samples while using the BAM files in ``data/`` for other samples, check out the :ref:`conditional-module-behaviour-user` section.

.. _reference-files-workflow:

Reference Files Workflow
========================

The ``reference_files`` workflow serves many purposes. In general, it simplifies the deployment of lcr-modules for any reference genome on any computer system. This approach ensures that the steps required to generate any reference file are tracked in a snakefile, ensuring their reproducibility. This is achieved by:

1. Downloading the genome FASTA files and any additional reference files (*e.g.* Gencode transcript annotations). 

2. Converting the additional files to match the same chromosome system as the genome builds (*e.g.* UCSC vs NCBI vs Ensembl)

3. Generate the required reference files from what was downloaded using snakemake rules. 

The snippet below needs to be added to your project snakefile before you include any of the individual module snakefiles. It adds the ``reference_files`` workflow as a sub-workflow in your project snakefile. Essentially, you gain access to the snakemake rules in that workflow, which you can trigger by passing files to the ``reference_files()`` function. In this case, the module will ask for specific reference files (*e.g.* genome FASTA file, STAR index), and if the file doesn't exist, the ``reference_files`` sub-workflow will create them. For more details, check out the `Snakemake Sub-Workflows`_ documentation. 

.. code:: python

   subworkflow reference_files:
         snakefile:
            "<path/to/lcr-modules/workflows/reference_files/1.0/reference_files.smk>"
         configfile:
            "<path/to/lcr-modules/workflows/reference_files/1.0/config/default.yaml>"
         workdir:
            "</path/to/reference_directory/>"

The ``configfile`` field specifies the path to the default ``reference_files`` configuration YAML file, which contains the information required for the workflow to run. The ``genome_builds`` section is the most important part of this configuration file. It defines the details for each available genome build, including the download URL. This enables the portability and thus reproducibility of lcr-modules. Each genome build also has a version (*i.e.* ``grch37`` or ``grch38``) and a provider (*e.g.* ``ensembl``, ``ucsc``). This metadata allows the ``reference_files`` workflow to automatically convert between the chromosome names of different providers (*e.g.* with and without the ``chr`` prefix). 

The ``workdir`` field specifies where the reference files will be created. Optionally, you can set this to a location shared between multiple lcr-modules users to avoid duplicating reference data. If you are considering building a shared reference directory, you might want to consider pre-populating it using the ``prepare_reference_files.smk`` snakefile. This will generate all of the output files that the ``reference_files`` workflow can produce for every genome build listed in the configuration file mentioned above. If you only need one genome build, you can remove the unnecessary genome builds from the configuration file. A bash script is included in the repository to perform this task:

.. code:: bash

   workflows/reference_files/1.0/prepare_reference_files.sh </path/to/reference_directory> <num_cores>

One caveat with the ``reference_files`` workflow is that the rules therein don't have any names. This is due to a Snakemake limitation. Rule names had to be omitted because this workflow could be included in more than one module loaded by the user and Snakemake doesn't allow duplicate rule names. As a result, you will be numbered rules (*e.g.* ``1``, ``2``, etc.) in your snakemake logs, such as the example shown below:

.. code::

   Job counts:
	count	jobs
	1	1
	1	2
	1	3
	7	_manta_augment_vcf
	2	_manta_configure
	2	_manta_dispatch
	1	_manta_index_bed
	3	_manta_input_bam
	1	_manta_input_bam_none
	3	_manta_output_bedpe
	7	_manta_output_vcf
	2	_manta_run
	3	_manta_vcf_to_bedpe
	1	_star_input_fastq
	1	_star_output_bam
	1	_star_run
	1	_star_symlink_sorted_bam
	1	_star_symlink_star_bam
	1	all
	40

.. _lcr-modules-configuration:

lcr-modules Configuration
=========================

.. code:: python

   # lcr-modules configuration
   config["lcr-modules"]

One of snakemake's most useful features is the ability to separate the workflow logic in the snakefiles from the tuneable parameters in the configuration files. lcr-modules is configured using :ref:`module-configuration` and :ref:`shared-configuration`. All configuration relating to lcr-modules is stored under the ``"lcr-modules"`` key in the snakemake ``config`` variable. This way, there is no risk of messing up any existing configuration created by the user.

**Important:** For brevity, the configuration under the ``"lcr-modules"`` key will be referred to as "lcr-modules configuration".

.. _module-configuration:

Module Configuration
--------------------

.. code:: python

   # Module configuration
   config["lcr-modules"]["<module_name>"]

Each module in the `lcr-modules repository`_ comes bundled with a default configuration to get users started. This module configuration is stored in the ``config/default.yaml`` YAML file in module subdirectory. These YAML files load the module configuration under module name in the lcr-modules configuration. You can see this in the excerpt below, which is taken from the ``manta`` module:

.. code:: yaml

   lcr-modules:
      manta:
         inputs:
            # Available wildcards: {seq_type} {genome_build} {sample_id}
            sample_bam: null  # UPDATE
            sample_bai: null  # UPDATE
            augment_manta_vcf: "{SCRIPTSDIR}/augment_manta_vcf/1.0/augment_manta_vcf.py"
         # ...

The intent behind these module configuration files is that any field can be (and often should be) updated by the user. In fact, some fields **must** be updated before the module can be run. These are indicated by ``UPDATE`` comments in the default configuration file. In the above excerpt, the two input files ``sample_bam`` and ``sample_bai`` are set to ``null`` and labelled with ``UPDATE`` comments, indicating that these must be updated by the user.

**Important:** Before running any module, you must search for any ``UPDATE`` comments in the default configuration file. See :ref:`updating-configuration-values` for different approaches on how to override the default configuration for each module.

.. _shared-configuration:

Shared Configuration
---------------------

.. code:: python
 
   # Shared configuration
   config["lcr-modules"]["_shared"]

One of the components of the :ref:`lcr-modules-configuration` is the shared configuration. As the name implies, the purpose of this shared configuration is to provide some common and relevant information to all modules. To avoid clashing with module names, this configuration is stored under the ``"_shared"`` key. (:ref:`faq-underscore`) 

**Important:** The configuration of each module is "merged" with the shared configuration, and when there are conflicts, the module configuration takes precedence. In other words, everything under ``"_shared"`` is used as default values when configuring each module unless the :ref:`module-configuration` already has a value, which will override the shared value. To demonstrate this point, consider the following shared configuration and the module configuration before merging. Once they have been merged, you can see that ``key2`` now appears in the module configuration using the value from the shared configuration, whereas the value for ``key1`` didn't change.

.. code:: yaml

   # Shared configuration
   _shared:
      key1: "a"
      key2: "b"
   
   # Module configuration (before merging)
   module_x:
      key1: "x"
      key3: "z"
   
   # Module configuration (after merging)
   module_x:
      key1: "x"
      key2: "b"
      key3: "z"

This behaviour can be leveraged in a number of ways. For example, by setting ``unmatched_normal_id`` under the ``"_shared"`` key, you avoid having to specify that value for every module that performs paired analyses (assuming you have unpaired tumours). That said, if you want to use a different ``unmatched_normal_id`` for a subset of modules, you can override the shared value. Another useful instance is the sharing of the :ref:`sample-table` between modules. This way, the user doesn't have to repeatedly provide the same sample table to each module. Note that this is only possible because each module automatically filters the samples based on the sequencing data types (``seq_type``) listed in their :ref:`pairing-configuration`. 

Common Shared Configuration Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``repository``: This field specifies the file path for the cloned ``lcr-modules`` repository. This path can be relative to your project snakefile or absolute. **This parameter is required.**

- ``lcr-scripts``: This field specifies the file path for the cloned ``lcr-scripts`` repository. This path can be relative to your project snakefile or absolute. **This parameter is required.**

- ``root_output_dir``: This field specifies the directory (*e.g.* ``"results/"``) where the modules will produce their output files in their respective subdirectories (*e.g.* ``results/star-1.0/``, ``results/manta-2.0/``). Technically, this shared parameter is optional, as it will default to ``"results/"``. 

- ``scratch_directory``: This field specifies the directory where large intermediate files can be written without worry of running out of space or clogging snapshots/backups. If set to ``null``, the files will be output locally.

- ``pairing_config``: This field is optional, but it's useful for specifying the normal samples to use in paired analyses with unpaired tumours. See :ref:`handling-unpaired-tumours` for more details and :ref:`updating-configuration-values` for an example configuration file where this is provided.

.. _handling-unpaired-tumours:

Handling Unpaired Tumours
^^^^^^^^^^^^^^^^^^^^^^^^^

Each module has a pairing configuration (``pairing_config``) in their default configuration file. This ``pairing_config`` dictates which sequencing data types (``seq_type``) are supported by the module, whether the module runs in paired or unpaired mode for each ``seq_type``, and if so, how it performs these analyses for each ``seq_type``. This information is ultimately used by the :py:func:`oncopipe.generate_runs_for_patient` function when producing (or not) tumour-normal pairs. If you want to learn more, check out the :ref:`pairing-configuration` section in the developer documentation. 

That said, the user doesn't need to worry about the ``pairing_config`` unless they have unpaired tumour samples or they wish to configure modules for new sequencing data types (``seq_type``). If they have unpaired tumours, for each ``seq_type``, they need to specify which normal sample to use for paired analyses where an unmatched normal sample will be used instead of a matched normal sample. This is done by providing values for ``unmatched_normal_id``, as demonstrated in the :ref:`within-the-configuration-file` section. 

If the user wants to configure new sequencing data types, they should check out the :ref:`configuring-new-seqtypes` section. 

.. _updating-configuration-values:

Updating Configuration Values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _within-the-configuration-file:

Within a Configuration File
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you followed the :ref:`getting-started-user` instructions, you should have a section in your project configuration file for ``lcr-modules`` (with at least the ``_shared`` sub-section). One approach to updating configuration values is to add to this section. **Important:** One requirement for this to work is that you need to load your project configuration file **after** the default module configuration files. Again, if you followed the :ref:`getting-started-user` instructions, this should already be the case. 

By the way, there is nothing forcing you to store your project-specific configuration in the same file as the lcr-modules configuration. You can easily have a ``project.yaml`` file loaded near the beginning of your snakefile and a ``lcr-modules.yaml`` file loaded as described in the :ref:`getting-started-user` instructions.

One of the main limitations of this approach is that you are restricted to value types that can be encoded in YAML format. For the most part, this means numbers, strings and booleans organized into lists or dictionaries. In other words, this precludes the use of functions as values, such as `Input File Functions`_. If you need to specify functions, you will have to update configuration values :ref:`within-the-snakefile`, or use a hybrid approach.

The example YAML file below is taken from the `Demo Configuration`_. You can see that it includes a ``pairing_config`` under ``_shared`` to indicate which normal samples to use for unpaired tumours for paired analyses (see :ref:`handling-unpaired-tumours`). It also updates a number of configuration values for the ``star`` and ``manta`` modules. All of these fields were labelled with an ``UPDATE`` comment in the modules' respective default configuration file. The only exception is the ``scratch_subdirectories`` field for the ``star`` module, which was updated here to include the ``"mark_dups"`` subdirectory such that the final BAM files from the module are also stored in the scratch directory.

.. code:: yaml

   # Taken from lcr-modules/demo/config.yaml
   lcr-modules:

      _shared:
         lcr-modules: "../"
         lcr-scripts: "../../lcr-scripts"
         root_output_dir: "results/"
         scratch_directory: "scratch/"
         pairing_config:
            capture:
                  unmatched_normal_id: "TCRBOA7-N-WEX"

      star:
         inputs:
            sample_fastq_1: "data/{sample_id}.read1.fastq.gz"
            sample_fastq_2: "data/{sample_id}.read2.fastq.gz"
         reference_params:
            star_overhang: "99"
         scratch_subdirectories: ["star", "sort_bam", "mark_dups"]

      manta:
         inputs:
            sample_bam: "data/{sample_id}.bam"
            sample_bai: "data/{sample_id}.bam.bai"

.. _within-the-snakefile:

Within the Snakefile
^^^^^^^^^^^^^^^^^^^^

You can always update configurations values within the snakefile **after** the default configuration files have been loaded. The advantage of this approach is that you can update the value to anything that Python allows, including functions. This is incredibly powerful in snakemake thanks to `Input File Functions`_ and `Parameter Functions`_. Also, :py:mod:`oncopipe` includes some useful functions that make use of these snakemake features (*e.g.* :ref:`conditional-module-behaviour-user`).

The main drawback of this approach is that it can be rather verbose, not very elegant, and as a result, not as readable. For instance, the equivalent of the above YAML file using this approach would look like this:

.. code:: python

   config["lcr-modules"]["_shared"]["pairing_config"] = {
      "capture": {
         "unmatched_normal_id": "TCRBOA7-N-WEX"
      }
   }

   config["lcr-modules"]["star"]["inputs"]["sample_fastq_1"] = "data/{sample_id}.read1.fastq.gz"
   config["lcr-modules"]["star"]["inputs"]["sample_fastq_2"] = "data/{sample_id}.read2.fastq.gz"
   config["lcr-modules"]["star"]["reference_params"]["star_overhang"] = "99"
   config["lcr-modules"]["star"]["scratch_subdirectories"] = ["star", "sort_bam", "mark_dups"]

   config["lcr-modules"]["manta"]["inputs"]["sample_bam"] = "data/{sample_id}.bam"
   config["lcr-modules"]["manta"]["inputs"]["sample_bai"] = "data/{sample_id}.bam.bai"

Alternatively, some of the redundancy can be avoided by using the `Snakemake update_config() Function`_, as follows. However, this alternative approach isn't much better. It takes up as many (if not more) lines, especially if you format the code to be readable.

.. code:: python

   import snakemake as smk

   smk.utils.update_config(config["lcr-modules"]["_shared"], {
      "pairing_config": {
         "capture": {
            "unmatched_normal_id": "TCRBOA7-N-WEX"
         }
      }
   })

   smk.utils.update_config(config["lcr-modules"]["star"], {
      "inputs": {
         "sample_fastq_1": "data/{sample_id}.read1.fastq.gz",
         "sample_fastq_2": "data/{sample_id}.read2.fastq.gz"
      },
      "reference_params": {"star_overhang": "99"},
      "scratch_subdirectories": ["star", "sort_bam", "mark_dups"]
   })

   smk.utils.update_config(config["lcr-modules"]["star"], {
      "inputs": {
         "sample_bam": "data/{sample_id}.bam",
         "sample_bai": "data/{sample_id}.bam.bai"
      }
   })

If you want a simpler syntax, you can consider using the :ref:`convenience-set-functions`. That said, a good compromise might be to store as much of these configuration updates :ref:`within-the-configuration-file` (*i.e.* anything that isn't a function), and you can update values with functions :ref:`within-the-snakefile`.

.. _sample-table:

Sample Table
============

+---------------+----------+------------+---------------+--------------+
| sample_id     | seq_type | patient_id | tissue_status | genome_build |
+===============+==========+============+===============+==============+
| TCRBOA7-N-WEX | capture  | TCRBOA7    | normal        | grch37       |
+---------------+----------+------------+---------------+--------------+
| TCRBOA7-T-WEX | capture  | TCRBOA7    | tumour        | grch37       |
+---------------+----------+------------+---------------+--------------+
| TCRBOA7-T-RNA | mrna     | TCRBOA7    | tumour        | grch37       |
+---------------+----------+------------+---------------+--------------+

One of the requirements for using lcr-modules is a sample table. This format was selected for its flexibility. Each sample can be annotated with any amount of metadata, but for the purposes of lcr-modules, there are only a few :ref:`required-columns`. These columns allow the modules to understand the relationship between the samples, especially for tumour-normal pairing. 

These requirements are encoded in schemas, which are stored and versioned in ``schemas/``. These schemas are used in conjunction with `Snakemake Validation`_. If the sample table doesn't confirm to a schemas that is required by a module, the user will given an informative error message. For example, the list of :ref:`required-columns` below is encoded in the ``base-1.0.yaml`` schema (located in ``schemas/base/``). The list of schemas will grow as modules are added with specific metadata requirements (*e.g.* strandedness of an RNA-seq library for expression quantification). 

The only format requirement for the sample table is that it is a `pandas DataFrame`_ (*i.e.* ``pandas.DataFrame``). Hence, the format of the file on disk doesn't matter. If you wish to use the :py:func:`oncopipe.load_samples` convenience function, note that it defaults to parsing tab-delimited files, but this can be overriden using the ``sep`` argument. The advantage of using :py:func:`oncopipe.load_samples` is that it offers a straightforward method for renaming your columns to comply with the schema(s). See :ref:`renaming-columns` for examples. 

Entity–Relationship Model
-------------------------

Before describing the required columns, it is useful to consider the entities related to each sample, namely ``patient``, ``biopsy``, ``sample``, ``library``, ``dataset``, and ``alignment``. These entities relate to one another in the following ways:

   **Relationships between entities:** Each patient has one or more biopsies (*e.g.* a tumour biopsy and a blood draw; tumour FF and FFPE biopsies). Each biopsy has one or more nucleic acid samples (*e.g.* DNA and RNA). Each sample has one or more sequencing libraries constructed from its nucleic acid samples (*e.g.* whole genome and RNA sequencing libraries for a tumour FF sample). Each sequenced library produces a a set of sequencing reads (*i.e.* a dataset) with one or more alignments (*e.g.* an hg19 and hg38 alignments), although there is generally a “canonical” alignment if more than one exists and thus a one-to-one relationship between datasets and alignments.

While the term "sample" generally refers to nucleic acid samples, lcr-modules uses the term to refer to the units of data that serve as input for the module, *i.e.* usually sequencing data in the form of FASTQ or BAM files. In most projects, there is a simple one-to-one relationship between these files and nucleic acid samples. In more complex projects where nucleic acid samples have more than one data file, the sample IDs will need to incorporate information to prevent duplicates.

.. _required-columns:

Required Columns
----------------

Check out the :ref:`renaming-columns` section if your sample table has some of the required columns under different names. It also features a demonstration of the :py:func:`oncopipe.load_samples` convenience function you can use to load your TSV/CSV sample table. The :ref:`adding-columns` section is useful if you lack some of the required columns or can derive them from existing columns.

``seq_type`` – Sequencing data type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most common values for this column are ``genome`` (whole genome sequencing), ``mrna`` (RNA sequencing), ``capture`` (hybridization-capture or exome sequencing), and ``mirna`` (miRNA sequencing). While ``lcr-modules`` can handle any value for ``seq_type``, the modules are pre-configured for these common values. If you have new ``seq_type`` values, you can configure them for modules of interest; this is explained in the :ref:`configuring-new-seqtypes` section under :ref:`advanced-usage`.

``sample_id`` – Sample identifiers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every ``seq_type`` and ``sample_id`` pair must be unique. In other words, if a tumour sample was sequenced using different technologies (*e.g.* whole genome and RNA sequencing), you can use the same sample ID since eachdata file will have a different ``seq_type`` (*e.g.* ``genome`` and ``mrna``, respectively). On the other hand, if you have been naming your samples based on patient ID and you have tumour-normal pairs, you will need to differentiate their sample IDs (*e.g.* with "T" and "N" suffixes). Similarly, if the same tumour sample has both FF and FFPE data files, you will also need to differentiate their sample IDs (*e.g.* with "FF" and "FFPE" suffixes). 

``tissue_status`` – Tumour or normal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This column classifies the samples as either ``tumour`` (or ``tumor``) and ``normal``. This information is required for tumour-normal paired analyses such as somatic variant calling. If you lack a matched normal samples, most modules support being run with an unmatched normal sample with the obvious caveats that the results will not be as clean. Check out :ref:`handling-unpaired-tumours` for more information on how to achieve this.

``patient_id`` – Patient identifiers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This column groups samples that originate from the same patient, *i.e.* that share the same underlying germline sequence. This information is primarily used in conjunction with the ``tissue_status`` column to generate all possible tumour-normal pairs from the list of samples.

``genome_build`` – Reference genome build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This column is only required if you have alignment (*i.e.* samples) using different genome builds. Otherwise, ``lcr-modules`` will assume that the single set of reference data (*e.g.* ``lcr-modules/references/hg38.yaml``) that you load is the one to use.

.. _renaming-columns:

Loading and Renaming Columns
----------------------------

For your convenience, the :py:func:`oncopipe.load_samples` function is provided to easily load your samples as a `pandas DataFrame`_. By default, the function assumes tab-delimited files, but you can change this using the ``sep`` argument. The function can also convert some of your columns to lowercase using the ``to_lowercase`` argument, which is useful to comply with some of the schemas. By default, it converts the ``tissue_status`` column to lowercase. It thus becomes trivial to load a sample table.

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv")

If your sample table uses different column names than those listed in :ref:`required-columns`, you can use the :py:func:`oncopipe.load_samples` function to rename your columns. For example, let's say you already have a sample table, but the sample ID and patient ID columns are named ``sample`` and ``patient`` rather than ``sample_id`` and ``patient_id``. You can easily achieve this as follows:

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv", sample_id = "sample", patient_id = "patient")

Alternatively, if the column names in your sample table differ systematically from the expected column names, you can rename them by passing a function to the ``renamer`` argument. You can also pass an anonymous ``lambda`` function. For instance, if you use two-letter prefixes with a period delimiter to indicate which entity a column describes (*e.g.* ``pt.`` for patient-related columns, ``lb.`` for library-related columns, etc.), you can remove the prefix from all columns using a regular expression with the following code:

.. code:: python

   import re
   import oncopipe as op
   remove_prefix = lambda x: re.sub(r"[a-z]{2}\.", "", x)
   SAMPLES = load_samples("samples.tsv", renamer=remove_prefix)

.. _adding-columns:

Adding and Transforming Columns
-------------------------------

If your sample table is missing a required column that has the same value for every sample (*e.g.* ``genome_build``), you can easily add the missing column in your snakefile using standard pandas_ syntax as follows:

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv")
   SAMPLES["genome_build"] = "hg38"

On the other hand, if your sample table is missing a required column that has different values for different samples, you can handle this one of two ways. If you can derive the missing column from existing columns, you can use standard pandas_ syntax to fill in the missing column. Otherwise, you can always resort to manually adding the missing column in the sample table on disk. The example below shows how the pandas_ syntax can be used to derive a ``tissue_status`` column by checking whether the ``sample_id`` column ends with the letter "T".

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv")
   SAMPLES["tissue_status"] = SAMPLES["sample_id"].str.endswith("T").map({True: "Tumour", False: "Normal"})

A similar approach can be taken if you have the columns, but they are formatted differently. For instance, if you encoded your sequencing data types as ``WGS`` and ``Exome`` instead of ``genome`` and ``capture``, respectively, you can use the ``map()`` method to switch to the expected values, as follows:

.. code:: python

   import oncopipe as op
   SAMPLES = op.load_samples("samples.tsv")
   SAMPLES["seq_type"] = SAMPLES["seq_type"].map({"WGS": "genome", "Exome": "capture"})

.. _specify-the-input-samples:

Specify the Input Samples
-------------------------

Once you have a sample table, you need to inform the modules which samples they need to run on. Normally, this is accomplished by storing the customized sample table under the ``"samples"`` key in the :ref:`module-configuration`. Because each module automatically filters for samples whose ``seq_type`` appear in their respective :ref:`pairing-configuration`, the user doesn't need to worry about pre-filtering the samples. For example, the user doesn't need to filter for RNA-seq samples for the ``star`` RNA-seq alignment module. That said, if the user had RNA-seq samples they didn't want processed by the ``star`` module, they can remove the samples in question and set this pre-filtered sample table as the value for the ``"samples"`` key.

However, since most users will probably want to run all samples through all applicable modules, it is possible to avoid the step of setting the sample table for each module. To skip this step, you can simply set the full sample table under the ``"samples"`` key in the :ref:`shared-configuration` (``"_shared"``). This is the method used in the :ref:`getting-started-user` instructions. The :ref:`shared-configuration` section explains why this works. 

.. _snakemake-commands:

Snakemake Commands
==================

**Note:** Don’t forget to update any values in angle brackets (``<...>``).

Snakemake Profiles
------------------

The most convenient way of running snakemake is using `snakemake profiles <https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles>`__. Each profile contains a YAML file that dictates the default command-line options to use. This way, you don’t have to remember all those snakemake options.

GSC Snakemake Profiles
~~~~~~~~~~~~~~~~~~~~~~

Make sure you first install the custom GSC snakemake profiles using `these instructions <https://github.com/LCR-BCCRC/snakemake-profiles#installation>`__. Then, you can use each profile using `these commands <https://github.com/LCR-BCCRC/snakemake-profiles#usage>`__.

Explicit Commands
-----------------

If you prefer to spell out all of the command-line options in your snakemake commands, example commands are included below. These may eventually become out of sync with the above snakemake profiles. Feel free to compare with the list of arguments for `local usage <https://github.com/LCR-BCCRC/snakemake-profiles/blob/master/gphosts/config.yaml>`__ or `cluster usage <https://github.com/LCR-BCCRC/snakemake-profiles/blob/master/numbers/config.yaml>`__.

Local Usage
~~~~~~~~~~~

.. code:: bash

   # See below for determining <cores>
   nice snakemake --printshellcmds --use-conda --cores <cores> <targets>

Cluster Usage
~~~~~~~~~~~~~

.. code:: bash

   nice snakemake --cluster-sync "srun --partition=all --ntasks=1 --nodes=1 --output=none --error=none --job-name={rule} --cpus-per-task={threads} --mem={resources.mem_mb}" --max-jobs-per-second=5 --max-status-checks-per-second=10 --local-cores=1 --latency-wait=120 --jobs=1000 --default-resources="mem_mb=2000" --printshellcmds --use-conda <targets>

Extra information
-----------------

Determining Value for ``--cores``
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

GSC ``ulimit`` Setup
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

Creating ``nice`` Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You will notice that the ``snakemake`` commands below are all prepended with ``nice``. Briefly, this has the effect of lowering the priority of your snakemake process. Now, you’re probably wondering why would you ever want to do that. Granted, compute resources should be utilized on a first come, first served basis, but in practice, not every user will pay close attention to who is already running jobs on a server.

Ultimately, it doesn’t matter whether this act is intentional, an accident, or due to insufficient knowledge of how to manage shared compute resources. If someone launches a job that uses more cores than are available, your snakemake process will be competing for CPU time, and this will make both processes take longer to complete.

In this situation, we should fall back on the motto from the wise Michelle Obama: “When they go low, we go high.” In this case, we follow this rule quite literally, because the ``nice`` command will increase the “niceness” value of your snakemake process, which will cede CPU time to competing processes with lower (usually default) “niceness” values until they’re done.

Submitting Cluster Jobs Remotely
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to submit jobs to a cluster remotely via SSH. This could be useful in situations where you have quick jobs that you don’t want to submit to the cluster, but you also don’t want to run locally on the cluster head node. **Important:** This section assumes that you have SSH keys set up, allowing SSH login to the head node without entering a password.

The command below differs from the explicit command above simply by prepending the ``srun`` command in ``--cluster-sync`` with ``ssh <head_node>``, where ``<head_node>`` is the cluster head node where you run ``srun`` normally. You can now increase the value for ``--local-cores`` (see above for how to determine this value).

.. code:: bash

   nice snakemake --local-cores=<cores> --cluster-sync "ssh <head_node> srun --partition=all --ntasks=1 --nodes=1 --output=none --error=none --job-name={rule} --cpus-per-task={threads} --mem={resources.mem_mb}" --max-jobs-per-second=5 --max-status-checks-per-second=10 --latency-wait=120 --jobs=1000 --default-resources="mem_mb=2000" --printshellcmds --use-conda <targets>

.. _advanced-usage:

Advanced Usage
==============

.. _directory-placeholders-users:

Directory Placeholders
----------------------

When specifying any value in the module configuration, you can use the following shorthands as placeholders in the string. They will be replaced with the actual values dynamically. See the :ref:`conditional-module-behaviour-user`. section below for example usage.

-  ``{REPODIR}``: The ``lcr-modules`` repository directory. This corresponds to the ``repository`` value under ``_shared`` in the ``lcr-modules`` configuration.
-  ``{MODSDIR}``: The current module subdirectory. This corresponds to ``{REPODIR}/modules/<name>/<version>``.
-  ``{SCRIPTSDIR}``: The ``lcr-scripts`` repository directory. This corresponds to the ``lcr-scripts`` value under ``_shared`` in the ``lcr-modules`` configuration.

.. _convenience-set-functions:

Convenience Set Functions
-------------------------

The :ref:`getting-started-user` instructions demonstrate that everything is configured using the same snakemake ``config`` nested dictionary object, generally under the ``"lcr-modules"`` key. While transparent, it results in verbose code, such as:

.. code:: python

   config["lcr-modules"]["manta"]["inputs"]["sample_bam"] = "data/{sample_id}.bam"

Alternatively, you can use the so-called convenience “set functions” to simplify the code somewhat. In order to use them, you must first enable them. Behind the scenes, the snakemake ``config`` object is stored internally for easy access. Note that you only need to use the :py:func:`oncopipe.enable_set_functions` function once.

.. code:: python

   import oncopipe as op

   op.enable_set_functions(config)

The first set function you can use is :py:func:`oncopipe.set_samples()`, which sets the samples you want to use at the shared or module level. The first argument corresponds to the module name (or ``"_shared"``), and all subsequent arguments should be sample tables each formatted as a `pandas DataFrame`_. This function automatically concatenates the data frames that are provided. Here, ``SAMPLES`` is the complete sample table, whereas ``GENOMES`` and ``CAPTURES`` are sample subsets generated from ``SAMPLES`` using :py:func:`oncopipe.filter_samples()`.

.. code:: python

   import oncopipe as op
   
   op.set_samples("_shared", SAMPLES)
   op.set_samples("_shared", GENOMES, CAPTURES)

The second function you can use is :py:func:`oncopipe.set_input()`, which sets the given input for a module. Just like :py:func:`oncopipe.set_samples`, the first argument is the module name, but this function should not be used for ``"shared"``. The second argument is the name of the input file as listed in the module’s configuration file. Lastly, the third argument is the value you wish to provide for that input file, which generally is a string value containing the available wildcards (once again, as listed in the module’s configuration file). That said, you could provide a conditional value as described below in :ref:`conditional-module-behaviour-user`.

.. code:: python

   import oncopipe as op

   op.set_input("manta", "sample_bam", "data/{sample_id}.bam")

Currently, only ``samples`` and ``inputs`` can be updated using these convenience set functions, but a more general set function will be implemented soon.

.. _conditional-module-behaviour-user:

Conditional Module Behaviour
----------------------------

Sometimes, a parameter or input file depends on some sample attribute. In snakemake, this kind of conditional behaviour is usually achieved with `Input File Functions`_ and `Parameter Functions`_. :py:mod:`oncopipe` offers two convenience functions that utilize these features to conditionally provide values depending on the value of a wildcard or in a column of the sample table, namely :py:func:`oncopipe.switch_on_wildcard` and :py:func:`oncopipe.switch_on_column`. These functions are useful for both module users and module developers. Read their documentation for more details, which you can access by clicking the function names above. 

Essentially, both functions create switches, which return a value based on a sample attribute, which can be stored in a wildcard or a column in the sample table, respectively. For example, :py:func:`oncopipe.switch_on_wildcard` takes in two arguments. The first argument is the name of the wildcard that contains the sample attribute on which the behaviour should be based. The second argument is a dictionary mapping possible values for that wildcard to the values that should be returned in the event of each wildcard. This is best illustrated by the example below. For more information, you can also check out the :ref:`conditional-module-behaviour-dev` section in the developer documentation.

Here, I am updating the `Demo Snakefile`_ such that the ``manta`` module uses BAM files from different sources depending on the sequencing data type (``seq_type``). Specifically, I want to configure the module to use the BAM files in the ``data/`` directory unless the ``seq_type`` is ``mrna`` (*i.e.* RNA-seq data), in which case the module should use the BAM file produced by the ``star`` module. Since ``seq_type`` is a wildcard in our file names, we can use the :py:func:`oncopipe.switch_on_wildcard` function. 

For simplicity, I am only updating the ``sample_bam`` input file, but the same would have to be done for the ``sample_bai`` input file. Also, this code assumes you have ``import oncopipe as op`` somewhere beforehand in your snakefile.

.. code:: python

   MANTA_BAM_SWITCH = op.switch_on_wildcard("seq_type", {
      "genome": "data/{sample_id}.bam",
      "capture": "data/{sample_id}.bam",
      "mrna": "results/star-1.0/99-outputs/bam/mrna--{genome_build}/{sample_id}.bam",
   })

   config["lcr-modules"]["manta"]["inputs"]["sample_bam"] = MANTA_BAM_SWITCH

While this is a good start, the redundancy between ``genome`` and ``capture`` is less than ideal. To address this, we can use the ``_default`` key to, well, set a default value. In other words, we can achieve the same thing using the code below. 

.. code:: python

   MANTA_BAM_SWITCH = op.switch_on_wildcard("seq_type", {
      "_default": "data/{sample_id}.bam",
      "mrna": "results/star-1.0/99-outputs/bam/mrna--{genome_build}/{sample_id}.bam"
   })

   config["lcr-modules"]["manta"]["inputs"]["sample_bam"] = MANTA_BAM_SWITCH

In the code above, I knew where the ``star`` module was going to output the BAM file. Alternatively, you can use the useful ``rules`` variable in snakemake. This highlights another reason why it's sometimes useful to update configuration values :ref:`within-the-snakefile`. **Important:** This will only work if you include the ``star`` module before the code below. In other words, the code snippet will only work if you add it between the ``star`` include and the ``manta`` include. 

.. code:: python

   MANTA_BAM_SWITCH = op.switch_on_wildcard("seq_type", {
      "_default": "data/{sample_id}.bam",
      "mrna": rules._star_output_bam.output.bam
   })

   config["lcr-modules"]["manta"]["inputs"]["sample_bam"] = MANTA_BAM_SWITCH

Lastly, we can make use of the :ref:`convenience-set-functions` and simplify the code just a bit more. If we also update ``sample_bai``, the final code would look like the following, added after the ``star`` include and before the ``manta`` include. 

.. code:: python

   op.enable_set_functions(config)
   
   MANTA_BAM_SWITCH = op.switch_on_wildcard("seq_type", {
      "_default": "data/{sample_id}.bam",
      "mrna": rules._star_output_bam.output.bam
   })   
   MANTA_BAI_SWITCH = op.switch_on_wildcard("seq_type", {
      "_default": "data/{sample_id}.bam.bai",
      "mrna": rules._star_output_bam.output.bam + ".bai
   })

   op.set_input("manta", "sample_bam", MANTA_BAM_SWITCH)
   op.set_input("manta", "sample_bai", MANTA_BAI_SWITCH)

For more information, check out :ref:`conditional-module-behaviour-dev`.

.. _configuring-new-seqtypes:

Configuring New Sequencing Data Types
-------------------------------------

This section is a work in progress, but you should be able to get started by reading the :ref:`pairing-configuration` section in the developer documentation. It's not recommended to add new sequencing data types (``seq_type``) under the ``_shared`` key because that will trigger all modules to try to run on the new ``seq_type``. It's best to configure the ``seq_type`` on a module-by-module basis.
