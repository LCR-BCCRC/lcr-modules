.. _getting-started-for-developers:

Getting Started for Developers
==============================

These guidelines rely on Git and GitHub for features like branching and pull requests. If you want to learn Git, it’s suggested that you work through `this tutorial <https://hamwaves.com/collaboration/doc/rypress.com/index.html>`__.

**Note:** Don’t forget to update any values in angle brackets (``<...>``).

1. First, check the `repository <https://github.com/LCR-BCCRC/lcr-modules/tree/master/modules/>`__ to see if the module already exists. Then, check the `list of open issues <https://github.com/LCR-BCCRC/lcr-modules/issues?q=is%3Aopen+is%3Aissue+label%3Anew-module>`__ to see if the module has already been proposed. If it has and hasn’t been assigned yet, feel free to assign yourself. Otherwise, reach out to the assignee and find out how you can help.

      **Important:** Please make sure you’re assigned to a GitHub issue before you start developing the module to avoid duplicating efforts.

2. Research how to best design the module. There are many ways to design any given module. While there is no expectation that the first version is perfect, it is preferable that an honest attempt is made to collect feedback from people who’ve run the tool before and from the literature (*e.g.* benchmark studies).

3. Clone the ``lcr-modules`` repository if you don’t already have a local copy.

   .. code:: bash

      git clone https://github.com/LCR-BCCRC/lcr-modules.git

   ..

      | **Note:** If you have SSH keys set up with GitHub, you can clone with SSH instead. This will ensure you don’t need to provide your username and password every time you push to GitHub. For this, use this link instead:
      | ``git@github.com:LCR-BCCRC/lcr-modules.git``.

4. Install ``snakemake`` (5.4 or later), ``pandas``, ``cookiecutter``, and the custom ``oncopipe`` Python packages into your conda environment. If you also need Git, you can add ``git`` to the ``conda install`` command.

   .. code:: bash

      # Move into the `lcr-modules` repository
      cd lcr-modules
      # `snakemake-minimal` lacks extraneous dependencies and is faster to install
      conda install --satisfied-skip-solve 'snakemake-minimal>=5.4' 'pandas' 'cookiecutter'
      # This is installing from the `lcr-modules` repository clone
      pip install -e lcr-modules/oncopipe

5. Create a new branch from the ``master`` branch with the format ``module/<module_name>/1.0``. Generally, ``<module_name>`` should refer to the core or defining software tool being used in the module. For example, the branch name for contributing a STAR alignment module would be ``module/star/1.0``, even though samtools would be used to sort the BAM file and sambamba for marking duplicates.

   **Important:** Your ``<module_name>`` should only contain lowercase alphanumerical characters or underscores (*i.e.* no spaces).

   .. code:: bash

      git checkout master
      git pull --ff-only
      git checkout -b "module/<module_name>/1.0"
      git branch  # To confirm which branch you are on

6. Initialize your module with the template with the following command:

   .. code:: bash

      cookiecutter "template/" --output-dir 'modules/'

   You will be asked for the following information:

   -  ``module_name``: Short name consisting of only lowercase alphanumerical characters or underscores (*i.e.* no spaces). This should match ``<module_name>`` in the branch name.

   -  ``module_author``: Full name of the person who is writing this module.

   -  ``original_author``: Full name of the person who wrote the Snakefile or script used to create this module. If this Snakefile/script doesn’t exist, then this can simply be set to ``N/A``.

   -  ``input_file_type`` and ``output_file_type``: List the file type (preferably the file extension) of the input (*e.g.* ``bam``) and output (*e.g.* ``vcf``) files, respectively. If there is more than one input file type and/or more than one output file type, just list one of them for each option. This is simply used to format the template as a starting point. You’ll be able to add more output file types in the Snakefile. Each of these should only consist of lowercase alphanumerical characters or underscores (*i.e.* no spaces).

   -  ``module_run_per``: Possible values are “tumour” and “sample”. This determines whether the module is meant to be run once per tumour (*e.g.* variant calling modules) or once per sample regardless of tissue status (*e.g.* BAM alignment and processing). Additional options will be added later, such as “tumour_cohort” and “sample_cohort” for `level-3 modules <README.md#module-levels>`__.

   -  ``seq_type.genome``, ``seq_type.capture``, and ``seq_type.mrna``: Possible values are “paired”, “unpaired”, and “omit”. This determines which sequencing data types (``seq_type``) are meant as input for this module (whole genome, hybrid capture-based, and RNA sequencing, respectively). Select “omit” if a ``seq_type`` is not applicable for your module. If you selected “sample” for ``module_run_per``, then you should use “unpaired” here. If this is a “paired” analysis, you should select “tumour” for ``module_run_per``. Lastly, if you selected “tumour” for ``module_run_per``, you can select “paired” or “unpaired” depending on whether the module is meant to be run on tumour-normal pairs or not. By default, if you select “paired” and the user has tumours without matched normal samples, they will be expected to provide an unmatched normal sample to be used instead. This behaviour is determined by ``run_unpaired_tumours_with`` in ``config/default.yaml`` of the generated module.

   **Important:** While technically possible to create a new module without using the cookiecutter template, we recommend against it. The template is maintained to follow the latest best practices for ``lcr-modules``.

7. Once you’ve generated your module from the cookiecutter template, you should be able to find it under ``modules/<module_name>/1.0/``. The parts you need to update are annotated with ``TODO``. These can be found in the ``<module_name>.smk`` file and the ``CHANGELOG.md`` file. A more detailed checklist can be found `here <.github/PULL_REQUEST_TEMPLATE.md>`__. You will need to work through this checklist when you submit your module to ``lcr-modules`` through a pull request (described below).

Module Description
==================

Module structure
----------------

When you create a new module `using the template <#getting-set-up>`__, you obtain the following files:

.. code:: bash

   ❯ tree modules/<module_name>
   modules/<module_name>
   ├── 1.0
   │   ├── <module_name>.smk
   │   ├── config
   │   │   └── default.yaml
   │   ├── envs
   │   │   └── samtools-1.9.yaml -> ../../../../envs/samtools/samtools-1.9.yaml
   │   ├── etc
   │   └── schemas
   │       └── base-1.0.yaml -> ../../../../schemas/base/base-1.0.yaml
   └── CHANGELOG.md

-  **``<module_name>.smk``:** This Snakefile contains the rules defining the module. See `Module snakefile <#module-snakefile>`__ below for more details.
-  **``config/default.yaml``:** This configuration YAML file contains all of the user-configurable options, such as input files, conda environments, command-line options, cluster parameters, and the pairing configuration (*i.e.* whether/how to run samples as tumour-normal pairs).
-  **``envs/``:** This folder contains symlinks to individual conda environment YAML files from the ``envs/`` directory, which is found in the root of the repository. These conda environment are generally tool-specific (*e.g.* ``samtools``, ``star``). Symlinks are used to keep the repository lightweight and promote reuse of conda environments between modules.
-  **``etc/``:** This folder can contain any accessory files required to run the module, such as configuration files (see ``manta`` module for an example).
-  **``schemas/``:** This folder contains symlinks to individual schema YAML files from the ``schemas/`` directory in the root of the repository. These schemas determine the required columns in the samples table. Every module should have the ``base-1.0.yaml`` schema as a minimum requirement. For more information, check out the `Required sample metadata <#required-sample-metadata>`__ section below. Symlinks are used to keep the repository lightweight and promote reuse of schemas between modules.
-  **``CHANGELOG.md``:** This file contains the release notes for the module. These release notes should list the changes and the rationale for each change.

Module snakefile
----------------

| This section will describe the key components of a module snakefile. It uses the ``star`` module as an example. Note that ``CFG`` refers to the module-specific configuration. In the case of the ``star`` module, this would correspond to:
| ``config["lcr-modules"]["star"]``.

Module attribution
~~~~~~~~~~~~~~~~~~

This section simply lists the individuals who have contributed to the module in one way or another. The ``Original Author`` refers to the person who wrote the Snakefile or script that was adapted for the module. The ``Module Author`` refers to the person who either adapted a previously written Snakefile/script or created the module from scratch. Finally, the ``Contributors`` refers to the list of individuals who have contributed to the module over time, mainly through incremental version updates.

.. code:: python

   ##### ATTRIBUTION #####


   # Original Author:   Nicole Thomas
   # Module Author:     Bruno Grande
   # Contributors:      N/A

Module setup
~~~~~~~~~~~~

There are a few standard components for the module setup and some optional components. Importing standard modules such as ``os`` (for the ``os.remove()`` function) is optional. On the other hand, importing the ``oncopipe`` module is required because it offers a suite of functions that greatly simplify the process of developing modules and facilitate configuration by the user. For brevity, the module is commonly imported with ``import oncopipe as op``, which allows the functions to be accessible using the ``op`` prefix/namespace (*e.g.* ``op.as_one_line()``).

The ``op.setup_module()`` function call is also required. This function does most of the heavy-lifting behind the scenes to streamline the process of developing modules. To find out more about what it does, you can check out the function docstring with ``help(op.setup_module)`` after importing ``oncopipe``. The arguments are self-explanatory: ``name`` is the module name, ``version`` is the module version, and ``subdirectories`` is the output subdirectories, which will be numbered automatically by ``op.setup_module()``.

The first and last subdirectories must be ``inputs`` and ``outputs``, and they will be numbered as ``00-inputs`` and ``99-outputs``, respectively. You should name the subdirectories after the tool name or the process, whatever is more evocative and specific (*e.g.* ``star`` over ``align``, or ``mark_dups`` over ``picard``).

Also, it’s worth noting that ``lcr-modules`` use a variant of semantic versioning where major versions represent changes in the number of rules in the module (or changes in the relationship between rules), whereas minor versions reprsent changes in the configuration of the module (*e.g.* command-line parameters).

The ``include`` statement for the ``utils`` module is optional. For more information on the ``include`` statement, you can refer to the `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes>`__. The ``utils`` module contains rules that are generally useful (*e.g.* BAM file sorting, BAM file indexing). It is meant to be included into another module after it has been configured with ``op.setup_module()``. The reason for this is that ``utils.smk`` makes use of the ``CFG`` variable to make sure it doesn’t interfere with other modules.

Finally, the ``localrules`` statement is technically optional, but it is recommended to include it in every module. For more information, you can refer to the `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#local-rules>`__. Essentially, when snakemake submits jobs to a cluster, these rules are run locally instead. It is meant for quick rules (*e.g.* symlinking rules) that aren’t computationally intensive and could potentially get stuck in the cluster queue for much longer than they take to run.

.. code:: python

   ##### SETUP #####


   # Import standard modules
   import os

   # Import package with useful functions for developing analysis modules
   import oncopipe as op

   # Setup module and store module-specific configuration in `CFG`
   # `CFG` is a shortcut to `config["lcr-modules"]["star"]`
   CFG = op.setup_module(
       name = "star",
       version = "1.0",
       subdirectories = ["inputs", "star", "sort_bam", "mark_dups", "outputs"],
   )

   # Include `utils` module
   include: "../../utils/1.0/utils.smk"

   # Define rules to be run locally when using a compute cluster
   localrules:
       _star_input_fastq,
       _star_symlink_in_sort_bam,
       _star_symlink_in_mark_dups,
       _star_output_bam,
       _star_all,

Module rules
~~~~~~~~~~~~

Input and output rules
^^^^^^^^^^^^^^^^^^^^^^

The input and output rules serve a few purposes. First, they clearly define the entry and exit points of the module, making the module more modular and easier to tie different modules together. Second, they make it clear to anyone exploring the module output directory what the input files were and what the most useful output files (or deliverables) are. Third, by symlinking the most important files in subdirectories with the same name (*i.e.* ``99-outputs``), it makes it easier to archive those files (*e.g.* from scratch space to backed-up storage).

You will notice that the ``op.relative_symlink()`` function (from the ``oncopipe`` module) is used in the rules below rather than ``os.symlink()`` (from the ``os`` module). The different between the two function is explained `below <#what-is-the-difference-between-oprelative_symlink-and-ossymlink>`__.

Below is the input and output rules for the ``star`` module. Because STAR operates on paired FASTQ files, we actually need to symlink two files per sample. While this could have been achieved in two rules, it was simpler to implement as one shared rule. The output file symlinks both the BAM and BAM index (BAI) files at the same time since they need to travel together. Otherwise, I find it useful to output different file types in different subdirectories in ``99-outputs``; see the ``manta`` module for an example, where VCF and BEDPE files are stored separately.

.. code:: python

   rule _star_input_fastq:
       input:
           fastq_1 = CFG["inputs"]["sample_fastq_1"],
           fastq_2 = CFG["inputs"]["sample_fastq_2"],
       output:
           fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R1.fastq.gz",
           fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R2.fastq.gz",
       run:
           op.relative_symlink(input.fastq_1, output.fastq_1)
           op.relative_symlink(input.fastq_2, output.fastq_2)

   # The other rules, which are normally in between, were omitted

   rule _star_output_bam:
       input:
           bam = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.mdups.bam",
           bai = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.mdups.bam.bai"
       output:
           bam = CFG["dirs"]["outputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
           bai = CFG["dirs"]["outputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
       run:
           op.relative_symlink(input.bam, output.bam)
           op.relative_symlink(input.bai, output.bai)

In some situations, it is useful to have more than one input or output symlinking rule. The example below is taken from the ``manta`` module, which can run on paired and unpaired tumour samples. To minimize duplicated code, the same rule is used for running paired and unpaired tumour samples, but the ``--normalBam`` argument is omitted for unpaired tumours. Unfortunately, the value for the ``normal_id`` in the filename is ``None``, which causes snakemake to look for a ``None.bam`` input file. While it is technically possible to omit the normal input file for unpaired tumours using duplicated rules, it’s harder to maintain. Hence, I added a second input rule that simply created an empty ``None.bam`` file.

.. code:: python

   # Symlinks the input BAM files into the module output directory (under '00-inputs/').
   rule _manta_input_bam:
       input:
           sample_bam = CFG["inputs"]["sample_bam"],
           sample_bai = CFG["inputs"]["sample_bai"]
       output:
           sample_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
           sample_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
       run:
           op.relative_symlink(input.sample_bam, output.sample_bam)
           op.relative_symlink(input.sample_bai, output.sample_bai)


   # Create empty file for "no normal" runs (but this is ultimately omitted from downstream rules)
   rule _manta_input_bam_none:
       output:
           touch(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/None.bam")

Target rules
^^^^^^^^^^^^

Generally, the last rule of the module snakefile is the “master target rule”. This rule is usually named ``_<module_name>_all`` (*e.g.* ``_star_all``), and expands all of the output files (the files symlinked into ``99-outputs``) using either the samples table (``CFG["samples"]``) or the runs table (``CFG["runs"]``) depending on whether the module is run once per sample or once per tumour. The two examples below show a preview of each table and how each can be used in the target rule.

Using the samples table
'''''''''''''''''''''''

+---------------------------+----------+-------------------+---------------+--------------+
| sample_id                 | seq_type | patient_id        | tissue_status | genome_build |
+===========================+==========+===================+===============+==============+
| BLGSP-71-08-00508-01B-01R | mrna     | BLGSP-71-08-00508 | tumour        | hg38         |
+---------------------------+----------+-------------------+---------------+--------------+
| BLGSP-71-06-00166-01B-01R | mrna     | BLGSP-71-06-00166 | tumour        | hg38         |
+---------------------------+----------+-------------------+---------------+--------------+

In the example below, since STAR is run on all RNA-seq BAM file, we are using the samples table, which has been automatically filtered for samples whose ``seq_type`` appears in the module’s ``pairing_config``. For more information on the ``pairing_config``, check out the `README <README.md#pairing-configuration>`__. Note the use of the ``rules`` variable that snakemake automatically generates for retrieving the output files from previous rules in the module.

.. code:: python

   rule _star_all:
       input:
           expand(
               [
                   rules._star_output_bam.output.bam,
                   rules._star_output_bam.output.bai,
               ]
               zip,  # Run expand() with zip(), not product()
               seq_type=CFG["samples"]["seq_type"],
               genome_build=CFG["samples"]["genome_build"],
               sample_id=CFG["samples"]["sample_id"])

Using the runs table
''''''''''''''''''''

+-------------+---------------------------+---------------------------+-----------------+-----------------+-------------------+-------------------+----------------------+----------------------+---------------------+---------------------+
| pair_status | tumour_sample_id          | normal_sample_id          | tumour_seq_type | normal_seq_type | tumour_patient_id | normal_patient_id | tumour_tissue_status | normal_tissue_status | tumour_genome_build | normal_genome_build |
+=============+===========================+===========================+=================+=================+===================+===================+======================+======================+=====================+=====================+
| matched     | BLGSP-71-08-00508-01B-01D | BLGSP-71-08-00508-10A-01D | capture         | capture         | BLGSP-71-08-00508 | BLGSP-71-08-00508 | tumour               | normal               | hg38                | hg38                |
+-------------+---------------------------+---------------------------+-----------------+-----------------+-------------------+-------------------+----------------------+----------------------+---------------------+---------------------+
| unmatched   | BLGSP-71-06-00166-01B-01D | BLGSP-71-08-00508-10A-01D | capture         | capture         | BLGSP-71-06-00166 | BLGSP-71-08-00508 | tumour               | normal               | hg38                | hg38                |
+-------------+---------------------------+---------------------------+-----------------+-----------------+-------------------+-------------------+----------------------+----------------------+---------------------+---------------------+
| no_normal   | BLGSP-71-08-00508-01B-01R |                           | mrna            |                 | BLGSP-71-08-00508 |                   | tumour               |                      | hg38                |                     |
+-------------+---------------------------+---------------------------+-----------------+-----------------+-------------------+-------------------+----------------------+----------------------+---------------------+---------------------+
| no_normal   | BLGSP-71-06-00166-01B-01R |                           | mrna            |                 | BLGSP-71-06-00166 |                   | tumour               |                      | hg38                |                     |
+-------------+---------------------------+---------------------------+-----------------+-----------------+-------------------+-------------------+----------------------+----------------------+---------------------+---------------------+


In this second example, taken from the ``manta`` module, we can see how the runs table (``CFG["runs"]``) is used to define the targets. Because the runs table lists tumour-normal pairs, each column from the samples table is present, but they are prefixed with ``tumour_`` and ``normal_``. The only column that isn’t taken from the samples table is ``pair_status``, which described the relationship between the tumour-normal pair. Generally, this can be ``matched`` if the tumour and normal samples come from the same patient; ``unmatched`` if the two samples come from different patients; and ``no_normal`` if there is no normal paired with the tumours.

It’s worth noting that the output rule being expanded is ``_manta_dispatch`` rather than ``_manta_output_vcf`` and ``_manta_output_bedpe``. The reason for this is technical, but briefly, it is because an input file function in the ``_manta_dispatch`` rule determines which files are converted into BEDPE format.

.. code:: python

   rule _manta_all:
       input:
           expand(
               [
                   rules._manta_dispatch.output.dispatched,
               ],
               zip,  # Run expand() with zip(), not product()
               seq_type=CFG["runs"]["tumour_seq_type"],
               genome_build=CFG["runs"]["tumour_genome_build"],
               tumour_id=CFG["runs"]["tumour_sample_id"],
               normal_id=CFG["runs"]["normal_sample_id"],
               pair_status=CFG["runs"]["pair_status"])

Other rules
^^^^^^^^^^^

Every other rule serve to complete the module. These other rules can vary considerably in scope. Therefore, below is a list of guiding principles to follow when designing these rules. These principles simply make it easier for users to achieve what they want. If one of these guidelines gets in the way of designing your module, feel free to employ a different approach, ideally not at the cost of flexibility for the user.

An example rule that follows most of these principles is included below (taken from the ``star`` module).

1.  Each rule should only consist of one command, unless the rule uses standard tools like ``gzip`` for additional commands. Otherwise, split into multiple rules, optionally connected using ``pipe()`` or ``temp()`` to avoid intermediate files.

       This guideline ensures that rules are modular and can easily be rearranged by the user. It also enables tool-specific conda environments (*e.g.* ``samtools``, ``star``) to be used, which is not possible is more than one tool is used in a rule.

2.  For ``input`` files, use ``rules`` references to previous output (or input) files wherever possible.

       These ``rules`` references minimizes the risk that two files get out of sync, *e.g.* if you update an upstream output file and forget to update every downstream occurrence of that file.

3.  Reference data should be provided as input files and ideally have rules in the ``reference_files`` workflow so they can be generated automatically. If a reference file has parameters, these can be exposed to the user under the ``reference_params`` section in the module configuration.

       Having reference data as input files ensures that rules are re-run if the reference data is updated. For more information on the ``reference_files`` workflow, check out the `Reference Files <#reference-files>`__ section below.

4.  The ``output`` (and ``input``) files should use values in the ``CFG["dirs"]``, which correspond to the subdirectory names provided to ``setup_module()``.

       This allows the user to easily adjust the output directory for the entire module.

5.  Avoid using non-standard wildcards. The standard wildcards for sample-based modules are: ``seq_type``, ``genome_build``, and ``sample_id``. The standard wildcards for tumour-based modules are: ``seq_type``, ``genome_build``, ``sample_id``, ``tumour_id``, and ``normal_id``.

       Adding new wildcards makes it hard to connect different modules together. For example, if module A adds an ``ffpe_status`` wildcard and module B depends on module A, module B will have to include ``ffpe_status`` as a wildcard, even though it’s not relevant to module B. You can thus see how this would result in the steady accumulation of wildcards. To change the behaviour of a module/rule based on sample metadata, see the `Condition rule behaviour <#condition-rule-behaviour>`__ section below.

6.  For ``log`` files, use the corresponding subdirectory names in ``CFG["logs"]``.

       The directories in ``CFG["logs"]`` are automatically timestamped, which allows the log files from each run to be stored separately for posterity.

7.  Store ``stdout`` and ``stderr`` in separate ``log`` files, unless the tool outputs to ``stdout``, in which case only ``stderr`` needs to be stored.

       Storing ``stdout`` and ``stderr`` in separate files makes it easier to know what output came from where, and it prevent potential issues with truncated log files.

8.  Create an ``opts`` entry under ``param`` for all command-line options that are not linked to a ``{...}`` value, which are configured in the ``default.yaml`` file.

       As you can see in the example below, every option under ``shell`` is associated with a value taken from the rule (*e.g.* ``--genomeDir {input.index}``), whereas it completely lacks “standalone options” (*e.g.* ``--runMode alignReads``). This guideline is to allow the user to have absolute control over the parameterization of the command-line tool.

9.  Re-use (or provide) tool-specific conda environments for each rule needing one, which are configured in the ``default.yaml`` file. This can be skipped if the rule only uses standard UNIX tools (*e.g.* ``gzip``, ``awk``) or if it uses the ``run`` directive (instead of the ``shell`` directive).

       Conda environments simplify software installation for a module and ensure reproducibility by specifying tool versions. Even if a rule only uses standard UNIX tools, it might still be worth using the ``coreutils`` conda environment to avoid OS variations (*e.g.* GNU vs BSD for ``sed``).

10. Add the ``threads`` and ``resources`` (``mem_mb``) directives for all non-local rules, which are configured in the ``default.yaml`` file.

       These directives are essential for running the module on a compute cluster. The values should be as low as possible while ensuring that most jobs are run within a reasonable amount of time (to minimize time spent in the queue).

11. Use the ``shell`` directive for rules with the ``conda`` directive. Use the ``run`` directive instead if more complicated logic is required.

       The ``op.as_one_line()`` function is meant to be used with the triple-quoted (``"""``) strings for long commands. The benefits of using this function are: (1) spaces are automatically added at the end of each line; (2) double-quotes do not need to be escaped; and (3) cleaner commands that are easier to organize using indentation. For example, any pipes (``|``) or double-ampersands (``&&``) can be indented to indicate the separation between two commands.

.. code:: python

   rule _star_run:
       input:
           fastq_1 = rules._star_input_fastq.output.fastq_1,
           fastq_2 = rules._star_input_fastq.output.fastq_2,
           index = reference_files(
               "genomes/{genome_build}/star_index/star-2.7.3a" +
                   "/gencode-" + CFG["reference_params"]["gencode_release"] +
                   "/overhang-" + CFG["reference_params"]["star_overhang"]
           ),
           gtf = reference_files(
               "genomes/{genome_build}/annotations" +
                   "/gencode_annotation-" + CFG["reference_params"]["gencode_release"] + ".gtf"
           )
       output:
           bam = CFG["dirs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/Aligned.out.bam"
       log:
           stdout = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/star.stdout.log",
           stderr = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/star.stderr.log"
       params:
           opts = CFG["options"]["star"],
           prefix = CFG["dirs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/"
       conda:
           CFG["conda_envs"]["star"]
       threads:
           CFG["threads"]["star"]
       resources:
           mem_mb = CFG["mem_mb"]["star"]
       shell:
           op.as_one_line("""
           STAR {params.opts} --readFilesIn {input.fastq_1} {input.fastq_2}
           --genomeDir {input.index} --outFileNamePrefix {params.prefix}
           --runThreadN {threads} --sjdbGTFfile {input.gtf} 
           > {log.stdout} 2> {log.stderr}
           """)

Module cleanup
~~~~~~~~~~~~~~

Every module ends with a clean-up step. At the moment, this mainly consists of outputting the module configuration, including the samples and runs, to disk for future reference. These files are output in a timestampted directory in the ``logs/`` subdirectory. Additionally, this function will delete the ``CFG`` variable from the environment to ensure it does not interfere with other modules.

.. code:: python

   # Perform some clean-up tasks, including storing the module-specific
   # configuration on disk and deleting the `CFG` variable
   op.cleanup_module(CFG)

Module configuration
--------------------

One of the core principles of ``lcr-modules`` is configurability, and this is primarily achieved by storing anything that can be adjusted in a configuration file separate from the Snakefile. For most modules, there will be a single configuration file called ``default.yaml``. On the other hand, some modules might have multiple configuration files to account for different scenarios. For this reason, there is a ``config/`` subdirectory for each module where all of these configuration files live.

In theory, configuration YAML files can take on any structure. However, it helps both module users and developers to start with a standard structure. This also facilitates feature development. Below is a description of each section of a typical ``default.yaml`` file using the ``star`` module as an example.

Configuration features
~~~~~~~~~~~~~~~~~~~~~~

Configuration comments
^^^^^^^^^^^^^^^^^^^^^^

It’s important to note the comment system used in the configuration files, which is explained at the top of every configuration file generated by the cookiecutter template. This comment system is intended to promote self-documentation as opposed to having the developer maintain a separate ``README.md`` file describing the ``default.yaml`` file. This latter approach is prone to files becoming out of sync.

Instead, every user is expected to read through the module configuration file and pay special attention to any lines commented out with ``#!``. They generally mean that some form of intervention is required from the user before the user can run the module. An example can be seen below in `Configuring options <#configuring-options>`__. On the other hand, ``#?`` comments generally do not require user intervention, but they might provide a means to adjust the behaviour of the module. Lastly, the ``##`` comments are regular comments, generally explaining the line(s) below them, including ``#!`` or ``#?`` comments.

.. code:: yaml

   ## Lines commented out with `#!` are required for the module to run
   ## Lines commented out with `#?` can optionally be user-configured
   ## Lines commented out with `##` act as regular comments

Directory placeholders
^^^^^^^^^^^^^^^^^^^^^^

Since the module developer won’t know where the ``lcr-modules`` (and ``lcr-scripts``, if applicable) repository will be located, one of the features of the ``setup_module()`` function in ``oncopipe`` is to replace the following directory placeholders with their actual values. This way, you can specify file paths relative to these directories. See the README for the list of `available placeholders <README.md#directory-placeholders>`__.

Configuring header
~~~~~~~~~~~~~~~~~~

Each module configuration should fall under the ``lcr-modules`` and ``<module_name>`` (*e.g.* ``star``) keys. The ``lcr-modules`` top-level configuration key is considered reserved for use by modules in this project and the ``oncopipe`` package. This ensures that the module configuration is properly siloed and avoids clashes with other configuration set by the user.

.. code:: yaml

   lcr-modules:
       star:

Configuring input and reference files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Virtually all modules will have input files, and many will also require reference files. These are defined using the ``inputs`` and ``reference_params`` keys, respectively.

The input files will generally be set to ``null`` and labelled with ``# UPDATE`` since they need to be specified by the user. This can be done in the configuration file or in the Snakefile (see the `demo Snakefile <demo/Snakefile>`__ for an example). Either way, the available wildcards are usually listed in a comment. If not, you can always look at the wildcards in the output files of the rule using the ``inputs`` configuration section. In general, these are ``{seq_type}``, ``{genome_build}``, and ``{sample_id}``.

   One advantage of specifying the input files in the Snakefile (as opposed to in the configuration file) is that the user can provide an `input file function <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#functions-as-input-files>`__ rather than a string.

While conceptually similar to input files, reference files are handled differently in ``lcr-modules``. They are generally genome build–specific rather than sample-specific. Accordingly, they need to be generated separately. In the past, this was often done in a time-consuming ad hoc way where the commands used to generate the reference files were often not tracked. A ``reference_files`` workflow was developed as part of ``lcr-modules`` to streamline this process and promote reproducibility. Most reference files depend only on the genome build and thus required no intervention from the user since the ``genome_build`` is a standard wildcard. However, some reference files require additional parameterization (*e.g.* the amount of splice-junction overhang when building a STAR index). These parameters are exposed to the user under the ``reference_params`` section. Some parameters are so important that they will be commented out with ``#!`` to require user intervention, such as the ``star_overhang`` parameter in the example below.

For more information on the approach taken in ``reference_files`` and its benefits and limitations, check out the section on `Reference files <#reference-files>`__.

.. code:: yaml

           inputs:
               # The inputs can be configured here or in the Snakefile
               # Available wildcards: {seq_type} {genome_build} {sample_id}
               sample_fastq_1: "<path/to/sample.R1.fastq.gz>"  # UPDATE
               sample_fastq_2: "<path/to/sample.R2.fastq.gz>"  # UPDATE

           reference_params:
               # Ideally, `star_overhang` = max(read_length) - 1
               # STAR indices were precomputed for "74" and "99"
               star_overhang: "99"  # UPDATE
               # The Gencode release to use for the transcript annotation
               gencode_release: "33"

Configuring scratch subdirectories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``scratch_subdirectories`` section provides the user with the ability of storing intermediate files in a scratch directory. Essentially, the listed subdirectories, which must match the names provided to the ``subdirectories`` argument in ``op.setup_module()``, will be made into symlinks to corresponding directories in a scratch space. This scratch space is also specified by the user, generally with the ``scratch_directory`` key under ``_shared``.

Note that if you’ve already run your Snakefile, the subdirectories will already exist as actual directories and not symlinks. Accordingly, you will have to delete them before adding another entry to ``scratch_subdirectories``. Otherwise, you will run into an error.

.. code:: yaml

           scratch_subdirectories: ["star", "sort_bam"]

Configuring options
~~~~~~~~~~~~~~~~~~~

The ``options`` section specifies the command-line options for each tool used in the module (where such options exist). Generally, any command-line option not linked to a placeholder (*e.g.* ``{input}``, ``{output}``, ``{params}``) should be listed under the tool’s corresponding entry in ``options``. This provides the user with ultimate control over how the tool is run without having to deal with the Snakefile.

Even if a tool has no command-line options beyond those already used in the Snakefile, it is useful to include an entry under ``options`` with an empty string in case options appear in future versions of the tool. For example, if the user wants to use a command-line option available in a later version of a tool, they can update the conda environment (see `below <#configuring-conda-environments>`__) and replace the empty string under ``options`` with the new option, thus avoiding any editing of the underlying Snakefile.

In the example below, the command-line options for STAR are commented out using ``#!`` because they require user intervention. Specifically, the value provided to the ``--sjdbOverhang`` argument should match the value provided to the ``star_overhang`` key under ``reference_params`` earlier in the configuration file (see `above <#configuring-input-and-reference-files>`__). A comment explains the user intervention that is required.

.. code:: yaml

           options:
               ## The value for `--sjdbOverhang` must match `star_overhang` above
               #! star: >
               #!     --runMode alignReads
               #!     --twopassMode Basic
               #!     --genomeLoad NoSharedMemory
               #!     --readFilesCommand zcat
               #!     --outSAMtype BAM Unsorted
               #!     --outSAMattrIHstart 0
               #!     --chimOutType WithinBAM SoftClip
               #!     --chimSegmentMin 20
               #!     --sjdbOverhang <star_overhang>
               utils_bam_sort: ""
               utils_bam_markdups: ""
               utils_bam_index: "-b"

Configuring conda environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conda environments that power each module are listed under ``conda_envs``. These allow for specific versions of tools to be automatically installed, which facilitates reproducibility. Each module will specify a set of default versions of each tool. The user can update this conda environments (*e.g.* to use a more recent version), but this might break the module if there are backwards-incompatible changes to the tool’s command-line interface.

Each conda environment should ideally be tool-specific because that promotes re-use of environments between modules. Otherwise, commonly used tools such as ``samtools`` would be included in multiple module-specific environments. This also allows for easier tracking of the tool versions in the file names. This can only be achieved if each module rule is indeed only using one tool, which should be the case.

Note that Snakemake expects the paths to be relative to the Snakefile. This is automatically handled by the ``op.setup_module()`` function, so these paths are expected to be relative to the working directory. In the example below, you can see the ``{MODSDIR}`` `placeholder <#directory-placeholders>`__ being used such that the paths are portably regardless of where the user stores the ``lcr-modules`` repository (as long as ``repository`` is specified under ``_shared``).

.. code:: yaml

           conda_envs:
               star: "{MODSDIR}/envs/star-2.7.3a.yaml"
               samtools: "{MODSDIR}/envs/samtools-1.9.yaml"
               sambamba: "{MODSDIR}/envs/sambamba-0.7.1.yaml"

Configuring compute resources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many users will be launching the modules on a high-performance computing cluster. Hence, all non-local rules should have sensible default values for resources such as CPU (``threads``) and memory (``mem_mb``). These settings should strike a balance between the time spent waiting in the queue (with higher resource values) and the time spent running (with lower resource values).

-  **``threads``:** The number of logical cores to allocate. This number is typically passed to a command-line argument such as ``--threads`` or ``--cores``. Make sure to check the tool’s actual CPU usage. If it’s consistently lower or higher than the specified amount, consider adjusting the value.
-  **``mem_mb``:** The amount of memory to allocate in megabytes (MB). This number is usually best determined empirically based on actual tool runs. This can be done in a number of ways, including monitoring ``top``/``htop`` or inspecting “Maximum resident set size” when the command is prepended with ``/usr/bin/time -v``.

.. code:: yaml

           threads:
               star: 12
               utils_bam_sort: 12
               utils_bam_markdups: 12
               utils_bam_index: 6

           mem_mb:
               star: 40000
               utils_bam_sort: 12000
               utils_bam_markdups: 8000
               utils_bam_index: 4000

Configuring sequencing data types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``pairing_config`` section is where the module is configured to run for each sequencing data type (``seq_type``). For example, in the STAR module, the pairing configuration obviously lists ``mrna`` for RNA-seq samples. The `user documentation <README.md#pairing-configuration>`__ on pairing configuration provides a description of each parameter (*e.g.* ``run_paired_tumours``).

.. code:: yaml

           pairing_config:
               mrna:
                   run_paired_tumours: False
                   run_unpaired_tumours_with: "no_normal"
                   run_paired_tumours_as_unpaired: True

Advanced module features
========================

Required sample metadata
------------------------

Every module requires the samples table, which contains metadata on the samples being analyzed. The minimum set of columns expected by ``lcr-modules`` are the ``sample_id``, ``patient_id``, ``seq_type``, and ``tissue_status`` columns (see `README <README.md#required-columns>`__ for more info). These requirements are spelled out using schemas in YAML format. The base requirements can be found in ``schemas/base/base-1.0.yaml``.

Some modules will need additional metadata (*e.g.* the strandedness of RNA-seq libraries). These extra requirements should also be described in schema files. To promote modularity, each required column should have its own file to promote modularity. An exception can be made for a set of columns should always be present together. The new schemas should be stored in the shared ``schemas/`` directory and then symlinked into individual modules. Symlinks are used to keep the repository lightweight and promote reuse of schemas between modules.

An example single-column schema file can be found in ``schemas/ffpe_status/ffpe_status-1.0.yaml``, where as a multi-column schema file should look like the base schema, *i.e.* ``schemas/base/base-1.0.yaml``.

**Important:** Read the section below on `Conditional module behaviour <#conditional-module-behaviour>`__ for an explanation on why you should avoid adding new wildcards beyond the standard ones described `above <#other-rules>`__.

Conditional module behaviour
----------------------------

One size doesn’t always fit all, so modules sometimes have to tailor their behaviour based on sample attributes. Snakemake offers more than one avenue to implement these conditional behaviours. The simplest approach is to create parallel rules, which will handle samples differently based on the file names, potentially using wildcard constraints. However, this approach has two major issues.

First, the resulting parallel rules are mostly identical except for a few, often minor differences (*e.g.* a single command-line argument). This redundancy violates the `DRY principle <https://en.wikipedia.org/wiki/Don%27t_repeat_yourself>`__, making the module harder to maintain and more vulnerable to bugs. This pitfall can be avoided by merging the two rules and using the `Switch on wildcard value <#switch-on-wildcard-value>`__ function from ``oncopipe`` described below.

Second, it requires the module developer to encode the sample attributes in the file names. While this is not a severe limitation on its own, it complicates the task of connecting modules together because the file names in downstream modules will need to include every wildcard from upstream modules. This would not only lead to unsustainably long file names, but the file names of a module shouldn’t depend on which modules are upstream to ensure modularity. The accumulation of module-specific wildcards can be avoided using the `Switch on sample metadata <#switch-on-sample-metadata>`__ function from ``oncopipe`` described below.

   To give a specific example, let’s say the ``salmon`` module requires the strandedness of the RNA-seq samples, so this information is encoded in the file name, *e.g.* ``{sample_id}.{strandedness}.quant``. Once we have quantified gene expression in all RNA-seq samples, we wish to perform cohort-wide correction for library size. Unfortunately, we need to pull the information about strandedness from the sample metadata in order to find the ``salmon`` output files because it’s part of the file names, even though that information isn’t relevant to our library size correction module.

**Important:** The ``op.switch_on_wildcard()`` and ``op.switch_on_column()`` functions do not currently support `Directory placeholders <#directory-placeholders>`__. This `issue <https://github.com/LCR-BCCRC/lcr-modules/issues/27>`__ will track the implementation.

Switch on wildcard value
~~~~~~~~~~~~~~~~~~~~~~~~

You can use the ``op.switch_on_wildcard()`` function to dynamically set the value of an input file or parameter for a snakemake rule based on the value of a wildcard. The first argument (``wildcard``) is the name of the wildcard, and the second argument (``options``) is a dictionary mapping possible values for the wildcard to the corresponding values that should be returned.

This dictionary can make use of special keys. The most important one to note is the ``"_default"`` special key, whose associated value is selected if the wildcard value isn’t among the other keys. You should check out the function docstring with ``help(op.switch_on_wildcard)`` to find out about the other special keys. See `below <#what-does-the-underscore-prefix-mean>`__ for an explanation for the underscore prefix.

By default, the ``op.switch_on_wildcard()`` will replace any placeholders (using the same format as the ``shell`` directive; *e.g.* ``{wildcards.seq_type}``) with the actual values. This beheviour can be tweaked with the ``format`` (default = ``True``) and ``strict`` (default = ``False``) optional arguments. See the function docstring for more information on these optional arguments.

An example taken from the ``manta`` module is included below (only relevant parts are shown). Here, the ``_manta_configure`` rule needs to use a different configuration file based on the sequencing data type (``seq_type``). Specifically, we wish to provide the high-sensitivity configuration if the ``seq_type`` is RNA-seq (``mrna``) or capture-based sequencing (``capture``), or the default configuration otherwise. Accordingly, the first argument is ``"seq_type"``.

.. code:: python

   rule _manta_configure:
       input:
           config = op.switch_on_wildcard("seq_type", CFG["switches"]["manta_config"])

The second argument is a reference to the module configuration (``CFG``), specifically the ``switches`` section. Since YAML files are parsed as nested dictionaries, it is straightforward to store the mapping between wildcard values and desired return values in the ``default.yaml`` configuration file. The relevant part from the YAML file is included below.

.. code:: yaml

   lcr-modules:
     manta:
       switches:
         manta_config:
           _default: "{MODSDIR}/etc/manta_config.default.ini"
           mrna: "{MODSDIR}/etc/manta_config.high_sensitivity.ini"
           capture: "{MODSDIR}/etc/manta_config.high_sensitivity.ini"

``CFG["switches"]["manta_config"]`` contains the dictionary representation of the ``manta_config`` section from the YAML file shown above. You can see how the ``"_default"`` special key is being used here (see `above <#switch-on-wildcard-value>`__ for more info) as well as the ``{MODSDIR}`` placeholder for the module subdirectory (see `above <#directory-placeholders>`__ for more info).

.. code:: python

   # This is the dictionary stored in `CFG["switches"]["manta_config"]`
   {
       '_default': '{MODSDIR}/etc/manta_config.default.ini',
       'mrna': '{MODSDIR}/etc/manta_config.high_sensitivity.ini',
       'capture': '{MODSDIR}/etc/manta_config.high_sensitivity.ini'
   }

Switch on sample metadata
~~~~~~~~~~~~~~~~~~~~~~~~~

As I mentioned `above <#conditional-module-behaviour>`__, adding wildcards for conditional behaviour in a Snakefile is unsustainable and goes against the core principle of modularity. One workaround is to query the metadata for each sample (or each tumour-normal pair) and to update the tool command accordingly. The approach is similar to a `Switch on wildcard value <#switch-on-wildcard-value>`__, but with a few notable differences.

The function to use is ``op.switch_on_column()``, where the first argument (``column``) is the column name, the second argument (``samples``) is the samples data frame (typically ``CFG["samples"]``), and the third argument (``options``) is a dictionary mapping possible values in the column to the corresponding values that should be returned. This dictionary follows the same structure as the `Switch on wildcard value <#switch-on-wildcard-value>`__. An additional albeit optional argument is called ``match_on``, which needs to be set to either ``"tumour"`` (default) or ``"normal"`` to determine whether the function uses the ``wildcards.tumour_id`` or ``wildcards.normal_id`` to look up a sample ID. The function will automatically use ``wildcards.seq_type`` to also filter on sequencing data type.

   At the moment, this function only works for tumour-based modules (*e.g.* paired variant calling). It should soon be generalized to also work with sample-based modules (*e.g.* STAR alignment). This issue is tracked `here <https://github.com/LCR-BCCRC/lcr-modules/issues/35>`__.

The code block below shows how we could achieve the same outcome using ``op.switch_on_column()`` for the example given in `Switch on wildcard value <#switch-on-wildcard-value>`__. The only difference other than the function name is the addition of the ``samples`` argument before providing the same ``options`` dictionary. By default, the function will use ``wildcards.tumour_id`` (and ``wildcards.seq_type``) to look up the sample in ``CFG["samples"]``. In practice, you would simply use ``op.switch_on_wildcard()`` since ``seq_type`` is available as a wildcard.

.. code:: python

   rule _manta_configure:
       input:
           config = op.switch_on_column("seq_type", CFG["samples"], CFG["switches"]["manta_config"])

Switch on file contents
-----------------------

The behaviour of some module depends on the contents (or existence) of input or intermediate files. The best way to address this is using `Snakemake checkpoints <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution>`__. They are a bit complicated to implement, but you can look at the ``manta`` module (version 1.0) for an example. Do note that checkpoints can be slow because the function using the checkpoint is run sequentially for each sample.

Frequently Asked Questions
==========================

What does the underscore prefix mean?
-------------------------------------

The underscore prefix is mainly used to avoid name conflicts. This convention is borrowed from Python. For instance, ``collections.namedtuple`` has an ``_asdict()`` method, where the underscore helps prevent clashes with user-defined attributes for the ``namedtuple``. For more examples in Python, check out this `blog post <https://medium.com/python-features/naming-conventions-with-underscores-in-python-791251ac7097>`__.

In ``lcr-modules``, the underscore prefix is used in a few areas. First, the name of every rule or function defined in a module starts with an underscore followed by the module name (*e.g.* ``_manta``). This minimizes the risk for clashing with other rule/function names defined elsewhere by the user, which isn’t allowed by Snakemake. Second, the underscore prefix is used for dictionary keys with special behaviour, such as the ``"_default"`` key in the ``op.switch_on_wildcard()`` `function <#switch-on-wildcard-value>`__. Third, the shared ``lcr-modules`` configuration is stored under the ``_shared`` key, which is done to avoid clashing with a potential module called ``shared``.

What is the difference between ``op.relative_symlink()`` and ``os.symlink()``?
------------------------------------------------------------------------------

Behind the scenes, ``op.relative_symlink()`` uses ``os.symlink()`` while ensuring that the symlinks are relative and correct regardless of the current working directory. This is equivalent to the ``-r`` option on modern version of the ``ln`` command-line tool.

Why am I running into a ``NameError: name 'CFG' is not defined`` exception?
---------------------------------------------------------------------------

Each module creates a ``CFG`` variable as a convenient but temporary pointer to the module configuration (*i.e.* ``config["lcr-modules"]["<module_name>"]``). Because each module uses this variable name, the ``op.cleanup_module()`` function deletes the variable to be safe. Hence, you will run into this ``NameError`` exception if some code tries to use ``CFG`` after it’s been deleted. If you use ``CFG`` in the rule directives that are evaluated when the module snakefile is parsed (*e.g.* ``input``, ``output``, ``log``, ``params``, etc.), it’s not an issue. However, if you use this variable in a function or ``run`` directive, *i.e.* code that is run after the ``op.cleanup_module()`` function is run, you will get the error above. You can fix this error by adding this line of code before using the ``CFG`` variable, which recreates the variable in a local scope:

.. code:: python

   # Replace <module_name> with the actual module name (e.g., `star`)
   CFG = config["lcr-modules"]["<module_name>"]

How do I specify the available memory per thread for a command-line tool?
-------------------------------------------------------------------------

The ``mem_mb`` resource is meant to represent the total amount of memory used by all threads of a given process. Some tools have command-line arguments allowing the user to specify the amount of memory they can use, such as any Java-based application (*i.e.* using ``-Xmx``). In some cases, the tool expects the amount of memory per thread (*e.g.* ``samtools sort``), whereas ``resources.mem_mb`` represents the total amount of memory. `Arithmetic expansion <https://www.shell-tips.com/bash/performing-math-calculation-in-bash#using-arithmetic-expansion-with-or>`__ in Bash allows you to circumvent this issue as long as you are dealing with integers, which should be the case with ``threads`` and ``mem_mb``. For example, here’s how you would divide two integers and print the result: ``echo $((12000 / 12))``. We can leverage the same syntax within the ``shell`` directive of a Snakemake rule. The example below is taken from the ``samtools sort`` rule in the ``utils`` module.

.. code:: bash

   rule:
       ...
       shell:
           op.as_one_line("""
           samtools sort {params.opts} -@ {threads} -m $(({resources.mem_mb} / {threads}))M
           -T {params.prefix} -o {output.bam} {input.bam} > {log.stdout} 2> {log.stderr}
           """)
