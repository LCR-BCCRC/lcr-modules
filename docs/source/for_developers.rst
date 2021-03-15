.. include:: links.rst

.. _getting-started-dev:

Getting Started
===============

**Important:** Be sure to update any values in angle brackets (``<...>``). If you are new to Git, we recommend reading through this excellent `Git Tutorial`_.

1. First, check the `lcr-modules repository`_ to see if the module already exists under the ``modules/`` directory. Then, check the `lcr-modules open issues`_ to see if the module has already been proposed. If so, reach out to the assignee if there is one, or assign yourself if it's unassigned. Otherwise, create a new issue for the module and assign yourself.

      **Important:** Please make sure you’re assigned to a GitHub issue before you start developing the module to avoid duplicating efforts.

2. Research how to best design the module. While there is no expectation that the first version is perfect, it is preferable that an honest attempt is made to collect feedback both from people with experience with the tools involved and from the literature (*e.g.* benchmarking studies).

3. Install ``conda`` (Python version 3.6 or later) with Miniconda_ or Anaconda_. Ideally, create a conda environment specific for lcr-modules development. Then, install the ``cookiecutter`` package and, if need be, Git using conda.

   .. code:: bash

      # Optional: Create a conda environment for lcr-modules development
      # conda create -n lcr-modules "python>=3.6"
      # conda activate lcr-modules
      
      conda install -c conda-forge cookiecutter

4. Clone the `lcr-modules repository`_ and the `lcr-scripts repository`_.

   .. code:: bash

      git clone https://github.com/LCR-BCCRC/lcr-modules.git
      git clone https://github.com/LCR-BCCRC/lcr-scripts.git

4. Install the custom :py:mod:`oncopipe` python package included in the `lcr-modules repository`_, which will also install dependencies such as ``snakemake`` and ``pandas``.

   .. code:: bash

      cd lcr-modules/
      pip install -e oncopipe/

5. Create a new branch from the ``master`` branch with the format ``module/<module_name>/1.0``, where ``<module_name>`` refers to the core or defining software tool being used in the module. **Important:** Your ``<module_name>`` should only contain lowercase alphanumerical characters or underscores (*i.e.* no spaces).

   .. code:: bash

      git checkout master  # Make sure you're on the master branch
      git pull --ff-only   # Pull the latest changes from GitHub
      git checkout -b "module/<module_name>/1.0"  # Create new branch
      git branch  # Confirm you're on the new branch (with the asterisk)
      git push -u origin "module/<module_name>/1.0"

6. Create a new module based on the :ref:`module-template`. Check out the :ref:`module-template` section for details on the fields requested during module creation. 

   .. code:: bash

      cookiecutter "template/" --output-dir 'modules/'
      git add modules/<module_name>/1.0/
      git commit -m "Add initial draft of <module_name> version 1.0"
      git push origin "module/<module_name>/1.0"

7. Update the basic module created from the template, which can be found under ``modules/<module_name>/1.0/``. Parts that need to be updated are indicated by ``TODO`` comments. You can use the `New Module Checklist`_ as a guide. These changes should be regularly committed to Git and pushed to GitHub.

   .. code:: bash

      git add <files>
      git commit -m "<commit message>"
      git push origin "module/<module_name>/1.0"

   .. note::

      **Testing Your Module**

      There are different methods testing your module. One approach would be to leverage the `Demo Project`_ and the associated test data. Adding to the `Demo Snakefile`_ should be self-explanatory. This method works if your module operates on the kind of samples included in the `Test Data`_. 

      Another approach to consider is testing on your own data. This might be a good way to follow up on successful tests on the `Demo Project`_, which confirm that the syntax of the module works. Running it on a larger dataset will confirm that the commands work on a variety of samples and that the output is sensible.

8. When you are done with your module, commit any remaining changes and merge the `master` branch into your module branch. You shouldn't have any merge conflicts since any new files should be under new versions.

   .. code:: bash
   
      git merge master
      git push origin "module/<module_name>/1.0"

9. Submit a pull request (PR) with your module branch. After pushing to GitHub, you should be able to see a green button to create a new pull request on the `lcr-modules repository`_ page pushing the merge commit. If not, you should be able to create one for your branch on the `lcr-modules active branches`_ page.

10. Work through the checklist that will appear when you open the PR. Once this checklist is done, you can request someone to review your PR. They can test the module if they have time and/or provide feedback on its design. Finally, once the reviewer(s) are happy, the PR can be merged. Congratulations! 

.. _module-template:

Module Template
===============

While it is technically possible to create a new module without using the module template, it's not recommended because using the template will ensure that you are following the latest best practices for lcr-modules. 

When you run the command listed in the :ref:`getting-started-dev` instructions, you will be asked for the following information:

- ``module_name``: This field specifies the short name for the module. The value provided here should match the value used for ``<module_name>`` in the branch name when following the :ref:`getting-started-dev` instructions. 

   **Important:** This field should only consist of lowercase alphanumerical characters or underscores (*i.e.* no spaces).

- ``module_author``: This field specifies the full name of the person who will write the module (presumably the person entering this information).

- ``original_author``: This field specifies the full name of the person who originally wrote the Snakefile or script that is being used to inform the module. If the module is being written from scratch, this field can be set to ``N/A``.

- ``input_file_type`` and ``output_file_type``: These fields specify the file type of the input and output files, respectively. Generally, these values will be the file extensions (*e.g.* ``bam``, ``vcf``). 

   **Important**

   - If there is more than one input file type, just list one of them for now. The same applies for the output file type. You’ll be able to add more file types in the Snakefile based on the existing structure.

   - Each of these should only consist of lowercase alphanumerical characters or underscores (*i.e.* no spaces).

- ``module_run_per``: Possible values are ``tumour`` and ``sample``. This field determines whether the module is intended to be run once per tumour (*e.g.* variant calling modules) or once per sample regardless of tissue status (*e.g.* BAM alignment and processing).

   Additional options will be added later, such as ``tumour_cohort`` and ``sample_cohort`` for level-3 modules (see :ref:`what-are-modules` for more details).

-  ``seq_type.genome``, ``seq_type.capture``, and ``seq_type.mrna``: Possible values are ``omit``,``unpaired``, ``matched_only``, ``allow_unmatched``, and ``no_normal``, . These fields determine which sequencing data types (``seq_type``) are intended as input for the module and whether each ``seq_type`` is intended to be run in paired or unpaired mode, and if in paired mode, whether to allow unmatched pairs. Select ``omit`` if a ``seq_type`` is not applicable for the module or ``unpaired`` if you are running the module once per sample. For more information on the last three modes, check out the documentation for the :py:func:`oncopipe.generate_pairs` function.  The fields correspond to whole genome, hybrid capture-based, and RNA sequencing, respectively.

   **Important**

   - If you selected ``sample`` for ``module_run_per``, then you should use ``unpaired`` (or ``omit``) here. If this is a paired analysis, you should start over (cancel with Ctrl-C) and select ``tumour`` for ``module_run_per``. 
   
   - If you selected ``tumour`` for ``module_run_per``, you can select ``matched_only``, ``allow_unmatched``, or ``no_normal`` depending on whether the module is meant to be run on only matched tumour-normal pairs, on potentially unmatched tumour-normal pairs, or on tumours only. 

Module Description
==================

Module Structure
----------------

When you create a new module using the :ref:`getting-started-user` instructions, you obtain the following files:

.. code::

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

-  ``<module_name>.smk``: This Snakefile contains the rules defining the module. See :ref:`module-snakefile` below for more details.

-  ``config/default.yaml``: This configuration YAML file contains all of the user-configurable options, such as input files, conda environments, command-line options, cluster parameters, and the pairing configuration (*i.e.* whether/how to run samples as tumour-normal pairs).

-  ``envs/``: This folder contains symlinks to individual conda environment YAML files from the ``envs/`` directory, which is found in the root of the repository. These conda environment are generally tool-specific (*e.g.* ``samtools``, ``star``). Symlinks are used to keep the repository lightweight and promote reuse of conda environments between modules.

-  ``etc/``: This folder can contain any accessory files required to run the module, such as configuration files (see ``manta`` module version 2.0 for an example). For more details, check out the :ref:`module-accessory-files-and-scripts` section. 

-  ``schemas/``: This folder contains symlinks to individual schema YAML files from the ``schemas/`` directory in the root of the repository. These schemas determine the required columns in the samples table. Every module should have the ``base-1.0.yaml`` schema as a minimum requirement. For more information, check out the :ref:`required-sample-metadata` section below. Symlinks are used to keep the repository lightweight and promote reuse of schemas between modules.

-  ``CHANGELOG.md``: This file contains the release notes for the module. These release notes should list the changes and the rationale for each change.

.. _module-snakefile:

Module Snakefile
----------------

This section will describe the key components of a module snakefile. It uses the ``star`` module as an example. Note that ``CFG`` refers to the module-specific configuration. In the case of the ``star`` module, this would correspond to:

.. code:: python

   config["lcr-modules"]["star"]

Module Attribution
~~~~~~~~~~~~~~~~~~

This section simply lists the individuals who have contributed to the module in one way or another. The ``Original Author`` refers to the person who wrote the Snakefile or script that was adapted for the module. The ``Module Author`` refers to the person who either adapted a previously written Snakefile/script or created the module from scratch. Finally, the ``Contributors`` refers to the list of individuals who have contributed to the module over time, mainly through incremental version updates.

.. code:: python

   ##### ATTRIBUTION #####


   # Original Author:   Nicole Thomas
   # Module Author:     Bruno Grande
   # Contributors:      N/A

Module Setup
~~~~~~~~~~~~

There are a few standard components for the module setup and some optional components. Importing standard modules such as ``os`` (for the ``os.remove()`` function) is optional. On the other hand, importing the :py:mod:`oncopipe` module is required because it offers a suite of functions that greatly simplify the process of developing modules and facilitate configuration by the user. For brevity, the module is commonly imported with ``import oncopipe as op``, which allows the functions to be accessible using the ``op`` prefix/namespace (*e.g.* ``op.as_one_line()``).

The :py:func:`oncopipe.setup_module` function call is also required. This function does most of the heavy-lifting behind the scenes to streamline the process of developing modules. The arguments are self-explanatory: ``name`` is the module name, ``version`` is the module version, and ``subdirectories`` is the output subdirectories, which will be numbered automatically by :py:func:`oncopipe.setup_module`.

The first and last subdirectories must be ``inputs`` and ``outputs``, and they will be numbered as ``00-inputs`` and ``99-outputs``, respectively. You should name the subdirectories after the tool name or the process, whatever is more evocative and specific (*e.g.* ``star`` over ``align``, or ``mark_dups`` over ``picard``).

Also, it’s worth noting that ``lcr-modules`` use a variant of semantic versioning where major versions represent changes in the number of rules in the module (or changes in the relationship between rules), whereas minor versions reprsent changes in the configuration of the module (*e.g.* command-line parameters).

The ``include`` statement for the ``utils`` module is optional. For more information on the ``include`` statement, you can refer to the `Snakemake Includes`_ documentation. The ``utils`` module contains rules that are generally useful (*e.g.* BAM file sorting, BAM file indexing). It is meant to be included into another module after it has been configured with :py:func:`oncopipe.setup_module`. The reason for this is that ``utils.smk`` makes use of the ``CFG`` variable to make sure it doesn’t interfere with other modules.

Finally, the ``localrules`` statement is technically optional, but it is recommended to include it in every module. For more information, you can refer to the `Snakemake Local Rules`_ documentation. Essentially, when snakemake submits jobs to a cluster, these rules are run locally instead. It is meant for quick rules (*e.g.* symlinking rules) that aren’t computationally intensive and could potentially get stuck in the cluster queue for much longer than they take to run.

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

Module Rules
~~~~~~~~~~~~

Input and Output Rules
^^^^^^^^^^^^^^^^^^^^^^

The input and output rules serve a few purposes. First, they clearly define the entry and exit points of the module, making the module more modular and easier to tie different modules together. Second, they make it clear to anyone exploring the module output directory what the input files were and what the most useful output files (or deliverables) are. Third, by symlinking the most important files in subdirectories with the same name (*i.e.* ``99-outputs``), it makes it easier to archive those files (*e.g.* from scratch space to backed-up storage).

You will notice that the :py:func:`oncopipe.relative_symlink` function is used in the rules below rather than the standard ``os.symlink()`` function. The different between the two function is explained here: :ref:`faq-symlink`.

Below is the input and output rules for the ``star`` module. Because STAR operates on paired FASTQ files, we actually need to symlink two files per sample. While this could have been achieved in two rules, it was simpler to implement as one shared rule. The output file symlinks both the BAM and BAM index (BAI) files at the same time since they need to travel together. Otherwise, I find it useful to output different file types in different subdirectories in ``99-outputs``; see the ``manta`` module for an example, where VCF and BEDPE files are stored separately. In this specific example, the output file rule also deletes an intermediate file. This is being done here to ensure that the downstream file exists before deleting the upstream file.

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
         bai = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.mdups.bam.bai",
         sorted_bam = str(rules._star_symlink_sorted_bam.input.bam)
      output:
         bam = CFG["dirs"]["outputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
      run:
         op.relative_symlink(input.bam, output.bam)
         op.relative_symlink(input.bai, output.bam + ".bai")
         os.remove(input.sorted_bam)
         shell("touch {input.sorted_bam}.deleted")

Target Rules
^^^^^^^^^^^^

Generally, the last rule of the module snakefile is the “master target rule”. This rule is usually named ``_<module_name>_all`` (*e.g.* ``_star_all``), and expands all of the output files (the files symlinked into ``99-outputs``) using either the samples table (``CFG["samples"]``) or the runs table (``CFG["runs"]``) depending on whether the module is run once per sample or once per tumour. The two examples below show a preview of each table and how each can be used in the target rule.

Using the Samples Table
'''''''''''''''''''''''

+---------------+----------+------------+---------------+--------------+
| sample_id     | seq_type | patient_id | tissue_status | genome_build |
+===============+==========+============+===============+==============+
| TCRBOA7-T-RNA | mrna     | TCRBOA7    | tumour        | grch37       |
+---------------+----------+------------+---------------+--------------+

In the example below, since STAR is run on all RNA-seq BAM file, we are using the samples table, which has been automatically filtered for samples whose ``seq_type`` appears in the module’s ``pairing_config``. For more information on the ``pairing_config``, check out :ref:`pairing-configuration`. Note the use of the ``rules`` variable that snakemake automatically generates for retrieving the output files from previous rules in the module.

.. code:: python

   rule _star_all:
      input:
         expand(
            str(rules._star_output_bam.output.bam),
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])

Using the Runs Table
''''''''''''''''''''

+-------------+------------------+------------------+-----------------+-----------------+-------------------+-------------------+----------------------+----------------------+---------------------+---------------------+
| pair_status | tumour_sample_id | normal_sample_id | tumour_seq_type | normal_seq_type | tumour_patient_id | normal_patient_id | tumour_tissue_status | normal_tissue_status | tumour_genome_build | normal_genome_build |
+=============+==================+==================+=================+=================+===================+===================+======================+======================+=====================+=====================+
| matched     | TCRBOA7-T-WEX    | TCRBOA7-N-WEX    | capture         | capture         | TCRBOA7           | TCRBOA7           | tumour               | normal               | grch37              | grch37              |
+-------------+------------------+------------------+-----------------+-----------------+-------------------+-------------------+----------------------+----------------------+---------------------+---------------------+
| no_normal   | TCRBOA7-T-RNA    |                  | mrna            |                 | TCRBOA7           |                   | tumour               |                      | grch37              |                     |
+-------------+------------------+------------------+-----------------+-----------------+-------------------+-------------------+----------------------+----------------------+---------------------+---------------------+

In this second example, taken from the ``manta`` module, we can see how the runs table (``CFG["runs"]``) is used to define the targets. Because the runs table lists tumour-normal pairs, each column from the samples table is present, but they are prefixed with ``tumour_`` and ``normal_``. The only column that isn’t taken from the samples table is ``pair_status``, which described the relationship between the tumour-normal pair. Generally, this can be ``matched`` if the tumour and normal samples come from the same patient; ``unmatched`` if the two samples come from different patients; and ``no_normal`` if there is no normal paired with the tumours.

It’s worth noting that the output rule being expanded is ``_manta_dispatch`` rather than ``_manta_output_vcf`` and ``_manta_output_bedpe``. The reason for this is technical, but briefly, it is because an input file function in the ``_manta_dispatch`` rule determines which files are converted into BEDPE format.

.. code:: python

   rule _manta_all:
       input:
           expand(
               [
                   str(rules._manta_dispatch.output.dispatched),
               ],
               zip,  # Run expand() with zip(), not product()
               seq_type=CFG["runs"]["tumour_seq_type"],
               genome_build=CFG["runs"]["tumour_genome_build"],
               tumour_id=CFG["runs"]["tumour_sample_id"],
               normal_id=CFG["runs"]["normal_sample_id"],
               pair_status=CFG["runs"]["pair_status"])

.. _other-rules:

Other Rules
^^^^^^^^^^^

Every other rule serve to complete the module. These other rules can vary considerably in scope. Therefore, below is a list of guiding principles to follow when designing these rules. These principles simply make it easier for users to achieve what they want. If one of these guidelines gets in the way of designing your module, feel free to employ a different approach, ideally not at the cost of flexibility for the user.

An example rule that follows most of these principles is included below (taken from the ``star`` module).

1.  Each rule should only consist of one command, unless the rule uses standard tools like ``gzip`` for additional commands. Otherwise, split into multiple rules, optionally connected using ``pipe()`` or ``temp()`` to avoid intermediate files.

       This guideline ensures that rules are modular and can easily be rearranged by the user. It also enables tool-specific conda environments (*e.g.* ``samtools``, ``star``) to be used, which is not possible is more than one tool is used in a rule.

2.  For ``input`` files, use ``rules`` references to previous output (or input) files wherever possible. You should wrap any references to ``rules`` with ``str()``. 

       These ``rules`` references minimizes the risk that two files get out of sync, *e.g.* if you update an upstream output file and forget to update every downstream occurrence of that file. 
       
       The ``str()`` function ensures that the ``rules`` reference isn't considered as an explicit dependency on whatever rule is specified. Otherwise, users won't be able to provide an alternative rule to generate the input in question. 

3.  Reference data should be provided as input files and ideally have rules in the ``reference_files`` workflow so they can be generated automatically. If a reference file has parameters, these can be exposed to the user under the ``reference_params`` section in the module configuration.

       Having reference data as input files ensures that rules are re-run if the reference data is updated. For more information on the ``reference_files`` workflow, check out the :ref:`reference-files-workflow` section.

4.  The ``output`` (and ``input``) files should use values in the ``CFG["dirs"]``, which correspond to the subdirectory names provided to ``setup_module()``.

       This allows the user to easily adjust the output directory for the entire module.

5.  Avoid using non-standard wildcards. The standard wildcards for sample-based modules are: ``seq_type``, ``genome_build``, and ``sample_id``. The standard wildcards for tumour-based modules are: ``seq_type``, ``genome_build``, ``sample_id``, ``tumour_id``, and ``normal_id``.

       Adding new wildcards makes it hard to connect different modules together. For example, if module A adds an ``ffpe_status`` wildcard and module B depends on module A, module B will have to include ``ffpe_status`` as a wildcard, even though it’s not relevant to module B. You can thus see how this would result in the steady accumulation of wildcards. To change the behaviour of a module/rule based on sample metadata, see the :ref:`conditional-module-behaviour-dev` section below.

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

       The :py:func:`.as_one_line()` function is meant to be used with the triple-quoted (``"""``) strings for long commands. The benefits of using this function are: (1) spaces are automatically added at the end of each line; (2) double-quotes do not need to be escaped; and (3) cleaner commands that are easier to organize using indentation. For example, any pipes (``|``) or double-ampersands (``&&``) can be indented to indicate the separation between two commands.

.. code:: python

   rule _star_run:
      input:
         fastq_1 = str(rules._star_input_fastq.output.fastq_1),
         fastq_2 = str(rules._star_input_fastq.output.fastq_2),
         index = reference_files("genomes/{{genome_build}}/star_index/star-2.7.3a/gencode-{}/overhang-{}".format(
            CFG["reference_params"]["gencode_release"], CFG["reference_params"]["star_overhang"]
         )),
         gtf = reference_files("genomes/{{genome_build}}/annotations/gencode_annotation-{}.gtf".format(
            CFG["reference_params"]["gencode_release"]
         ))
      output:
         bam = CFG["dirs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/Aligned.out.bam"
      log:
         stdout = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/star.stdout.log",
         stderr = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/star.stderr.log"
      params:
         opts = CFG["options"]["star"],
         prefix = CFG["dirs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/",
         star_overhang = CFG["reference_params"]["star_overhang"]
      conda:
         CFG["conda_envs"]["star"]
      threads:
         CFG["threads"]["star"]
      resources:
         mem_mb = CFG["mem_mb"]["star"]
      shell:
         op.as_one_line("""
         STAR {params.opts} --readFilesIn {input.fastq_1} {input.fastq_2} --genomeDir {input.index} 
         --outFileNamePrefix {params.prefix} --runThreadN {threads} --sjdbGTFfile {input.gtf}
         --sjdbOverhang {params.star_overhang} > {log.stdout} 2> {log.stderr}
               &&
         rmdir {params.prefix}/_STARtmp
         """)

Module Cleanup
~~~~~~~~~~~~~~

Every module ends with a clean-up step. At the moment, this mainly consists of outputting the module configuration, including the samples and runs, to disk for future reference. These files are output in a timestampted directory in the ``logs/`` subdirectory. Additionally, this function will delete the ``CFG`` variable from the environment to ensure it does not interfere with other modules.

.. code:: python

   # Perform some clean-up tasks, including storing the module-specific
   # configuration on disk and deleting the `CFG` variable
   op.cleanup_module(CFG)

Module Configuration
--------------------

One of the core principles of ``lcr-modules`` is configurability, and this is primarily achieved by storing anything that can be adjusted in a configuration file separate from the Snakefile. For most modules, there will be a single configuration file called ``default.yaml``. On the other hand, some modules might have multiple configuration files to account for different scenarios. For this reason, there is a ``config/`` subdirectory for each module where all of these configuration files live.

In theory, configuration YAML files can take on any structure. However, it helps both module users and developers to start with a standard structure. This also facilitates feature development. Below is a description of each section of a typical ``default.yaml`` file using the ``star`` module as an example.

Configuration Features
~~~~~~~~~~~~~~~~~~~~~~

Requiring User Intervention
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure that anything that needs to be updated by the user contains ``__UPDATE__`` in the configuration file. You can see examples in the excerpts below taken from the ``star`` default configuration. If the string ``__UPDATE__`` is detected anywhere in the module configuration, an error will inform the user that they need to update a configuration field.

.. _directory-placeholders-dev:

Directory Placeholders
^^^^^^^^^^^^^^^^^^^^^^

Since the module developer won’t know where the ``lcr-modules`` (and ``lcr-scripts``, if applicable) repository will be located, one of the features of the ``setup_module()`` function in :py:mod:`oncopipe` is to replace the following directory placeholders with their actual values. This way, you can specify file paths relative to these directories. See the README for the list of :ref:`directory-placeholders-users`.

Configuring Header
~~~~~~~~~~~~~~~~~~

Each module configuration should fall under the ``lcr-modules`` and ``<module_name>`` (*e.g.* ``star``) keys. The ``lcr-modules`` top-level configuration key is considered reserved for use by modules in this project and the :py:mod:`oncopipe` package. This ensures that the module configuration is properly siloed and avoids clashes with other configuration set by the user.

.. code:: yaml

   lcr-modules:
      star:

Configuring Input and Reference Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Virtually all modules will have input files, and many will also require reference files. These are defined using the ``inputs`` and ``reference_params`` keys, respectively.

The input files will generally be set to ``__UPDATE__`` since they need to be specified by the user. This can be done in the configuration file or in the Snakefile (see the `Demo Snakefile`_ for an example). Either way, the available wildcards are usually listed in a comment. If not, you can always look at the wildcards in the output files of the rule using the ``inputs`` configuration section. In general, these are ``{seq_type}``, ``{genome_build}``, and ``{sample_id}``.

   One advantage of specifying the input files in the Snakefile (as opposed to in the configuration file) is that the user can provide `Input File Functions`_ rather than a string.

While conceptually similar to input files, reference files are handled differently in ``lcr-modules``. They are generally genome build–specific rather than sample-specific. Accordingly, they need to be generated separately. In the past, this was often done in a time-consuming ad hoc way where the commands used to generate the reference files were often not tracked. A ``reference_files`` workflow was developed as part of ``lcr-modules`` to streamline this process and promote reproducibility. Most reference files depend only on the genome build and thus required no intervention from the user since the ``genome_build`` is a standard wildcard. However, some reference files require additional parameterization (*e.g.* the amount of splice-junction overhang when building a STAR index). These parameters are exposed to the user under the ``reference_params`` section. Some parameters are so important that they will be commented out with ``#!`` to require user intervention, such as the ``star_overhang`` parameter in the example below.

For more information on the approach taken in ``reference_files`` and its benefits and limitations, check out the section on the :ref:`reference-files-workflow`.

.. code:: yaml

         inputs:
            # The inputs can be configured here or in the Snakefile
            # Available wildcards: {seq_type} {genome_build} {sample_id}
            sample_fastq_1: "__UPDATE__"
            sample_fastq_2: "__UPDATE__"

         reference_params:
            # Ideally, `star_overhang` = max(read_length) - 1
            # STAR indices were precomputed for "74" and "99"
            star_overhang: "__UPDATE__"
            # The Gencode release to use for the transcript annotation
            gencode_release: "33"

Configuring Scratch Subdirectories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``scratch_subdirectories`` section provides the user with the ability of storing intermediate files in a scratch directory. Essentially, the listed subdirectories, which must match the names provided to the ``subdirectories`` argument in :py:func:`oncopipe.setup_module`, will be made into symlinks to corresponding directories in a scratch space. This scratch space is also specified by the user, generally with the ``scratch_directory`` key under ``_shared``.

Note that if you’ve already run your Snakefile, the subdirectories will already exist as actual directories and not symlinks. Accordingly, you will have to delete them before adding another entry to ``scratch_subdirectories``. Otherwise, you will run into an error.

.. code:: yaml

         scratch_subdirectories: ["star", "sort_bam"]

Configuring Options
~~~~~~~~~~~~~~~~~~~

The ``options`` section specifies the command-line options for each tool used in the module (where such options exist). Generally, any command-line option not linked to a placeholder (*e.g.* ``{input}``, ``{output}``, ``{params}``) should be listed under the tool’s corresponding entry in ``options``. This provides the user with ultimate control over how the tool is run without having to deal with the Snakefile.

Even if a tool has no command-line options beyond those already used in the Snakefile, it is useful to include an entry under ``options`` with an empty string in case options appear in future versions of the tool. For example, if the user wants to use a command-line option available in a later version of a tool, they can update the conda environment (see :ref:`configuring-conda-environments`) and replace the empty string under ``options`` with the new option, thus avoiding any editing of the underlying Snakefile.

In the example below, you can see that any command-line options associated with a snakemake parameter (*e.g.* ``--sjdbOverhang``, ``--runThreadN``) or an input/output file (*e.g.* ``--readFilesIn``, ``--outFileNamePrefix``) are not included here. Instead, they reside in the ``star`` snakefile. 

.. code:: yaml

         options:
            star:
               --runMode alignReads
               --twopassMode Basic 
               --genomeLoad NoSharedMemory 
               --readFilesCommand zcat 
               --outSAMtype BAM Unsorted
               --outSAMattrIHstart 0
               --chimOutType WithinBAM SoftClip
               --chimSegmentMin 20 
            utils_bam_sort: ""
            utils_bam_markdups: ""
            utils_bam_index: "-b"

You will also notice the various ``utils_bam_*`` fields. These correspond to rules in the ``utils`` module. For example, the following ``utils`` rule can index a BAM file, and you can see how it has an ``opts`` parameter that looks up the ``utils_bam_index`` field under ``options``. If it doesn't find a value, it defaults to ``"-b"``. In this case, the module developer exposed these fields to the user by including them in the ``star`` default configuration. In this case, the same default values as in the ``utils`` module were used, but that might not always be the case depending on the module.

.. code:: yaml

   # _utils_bam_index: Index a BAM file
   rule:
      input:
         bam = CFG["dirs"]["_parent"] + "{prefix}/{suffix}.bam"
      output:
         bam = CFG["dirs"]["_parent"] + "{prefix}/{suffix}.bam.bai"
      # ...
      params:
         opts = CFG["options"].get("utils_bam_index", "-b"),
         prefix = CFG["dirs"]["_parent"] + "{prefix}/{suffix}"
      # ...
      shell:
         op.as_one_line("""
         samtools index {params.opts} -@ {threads} 
         {input.bam} > {log.stdout} 2> {log.stderr}
         """)

.. _configuring-conda-environments:

Configuring Conda Environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conda environments that power each module are listed under ``conda_envs``. These allow for specific versions of tools to be automatically installed, which facilitates reproducibility. Each module will specify a set of default versions of each tool. The user can update this conda environments (*e.g.* to use a more recent version), but this might break the module if there are backwards-incompatible changes to the tool’s command-line interface.

Each conda environment should ideally be tool-specific because that promotes re-use of environments between modules. Otherwise, commonly used tools such as ``samtools`` would be included in multiple module-specific environments. This also allows for easier tracking of the tool versions in the file names. This can only be achieved if each module rule is indeed only using one tool, which should be the case.

Note that Snakemake expects the paths to be relative to the Snakefile. This is automatically handled by the :py:func:`oncopipe.setup_module` function, so the paths provided under ``conda_envs`` in the module configuration are expected to be relative to the working directory (usually where you run the ``snakemake`` command). In the example below, you can see the ``{MODSDIR}`` directory placeholder being used such that the paths are portably regardless of where the user stores the ``lcr-modules`` repository (as long as ``repository`` is specified under ``_shared``). For more information, check out :ref:`directory-placeholders-dev`.

.. code:: yaml

         conda_envs:
            star: "{MODSDIR}/envs/star-2.7.3a.yaml"
            samtools: "{MODSDIR}/envs/samtools-1.9.yaml"
            sambamba: "{MODSDIR}/envs/sambamba-0.7.1.yaml"

Creating Conda Environments
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a suggested workflow for creating conda environments for your modules. The commands are based on the example of creating a conda environment for the STAR aligner. 

1. Create a conda environment for the tool in question. For example, in the case of STAR, the following command would create a conda environment called ``test-star`` and install the latest version of the ``star`` package from the ``bioconda`` Anaconda channel along with its dependencies

   .. code:: bash

      conda create -c bioconda -n test-<star> <star>

2. Optionally, activate this conda environment before testing any STAR commands. This testing can occur in bash scripts or in a snakefile as long as the rule(s) don't activate their own conda environment, thus relying on the global shell environment.

   .. code:: bash

      conda activate test-<star>
      bash run_<star>_test_commands.sh ...

3. Create a subdirectory in ``envs/`` named after the tool once you are ready to add the new conda environments to the lcr-modules repository. Again, for the STAR environment, this would look like:

   .. code:: bash

      mkdir envs/<star>/

4. Determine the version of the tool that was installed by conda. You can usually achieve this by grepping the tool from the output of ``conda env export``. 

   .. code:: bash

      conda env export -n test-<star> --no-builds | grep <star>

5. Output the conda environment specification into a file named after the tool and its version. Be sure to use the ``--no-builds`` command-line option to omit build ID, which tend to cause problems recreating environments later. 

   .. code:: bash

      conda env export -n test-<star> --no-builds > envs/<star>/<star>-<2.7.3a>.yaml

6. Symlink this new conda environment YAML file into your module's ``envs/`` directory. The symlinks prevent the repository from becoming bloated with duplicate files. This approach also promotes the re-use of conda environments across modules.

   .. code:: bash

      cd modules/<star>/<1.0>/envs/
      ln -s ../../../../envs/<star>/<star>-<2.7.3a>.yaml ./

Configuring Compute Resources
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

.. _pairing-configuration:

Pairing Configuration
~~~~~~~~~~~~~~~~~~~~~

The ``pairing_config`` section is where the module is configured to run for each sequencing data type (``seq_type``). Two examples are included below to illustrate how the ``pairing_config`` is used. Check out the :ref:`pairing-configuration-options` section for more details on each field (*e.g.* ``run_paired_tumours``).

In this first example, we continue with the ``star`` module. Here, the pairing configuration only lists ``mrna`` (*i.e.* RNA-seq data) as a supported ``seq_type``. In the future, additional sequencing data types could be added, such as ``mirna`` for miRNA sequencing data. For ``mrna``, the ``star`` module is configured to run on all samples in unpaired mode. This is achieved by first disabling paired mode (``run_paired_tumours`` as ``False``) and then ensuring that any paired tumours are forced to run in unpaired mode (``run_paired_tumours_as_unpaired`` as ``True``). Setting ``run_unpaired_tumours_with`` to ``"no_normal"`` is meant to clarify that the unpaired tumours should be included; otherwise, they would be omitted since the default for ``run_unpaired_tumours_with`` is ``None``.

.. code:: yaml

         pairing_config:
            mrna:
               run_paired_tumours: False
               run_unpaired_tumours_with: "no_normal"
               run_paired_tumours_as_unpaired: True

This second example was taken from the ``manta`` module. As you can see, the module can handle ``genome``, ``capture``, and ``mrna`` data. It treats ``genome`` and ``capture`` data the same way, namely by allowing unpaired tumours to be analyzed using unmatched normals (as opposed to a truly unpaired analysis without a normal sample). Also, paired tumours are not unnecessarily run as unpaired. In contrast, ``mrna`` data is run specifically in an unpaired fashion without a normal sample because tumour RNA-seq alignments generally do not have matched normal RNA-seq data. 

.. code:: yaml

         pairing_config:
            genome:
               run_paired_tumours: True
               run_unpaired_tumours_with: "unmatched_normal"
               run_paired_tumours_as_unpaired: False
            capture:
               run_paired_tumours: True
               run_unpaired_tumours_with: "unmatched_normal"
               run_paired_tumours_as_unpaired: False
            mrna:
               run_paired_tumours: False
               run_unpaired_tumours_with: "no_normal"
               run_paired_tumours_as_unpaired: True

.. _pairing-configuration-options:

Pairing Configuration Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here's a brief description of each of the options that go into a ``pairing_config``. Here, the term "unpaired tumour" refers to tumours that lack a matched normal sample with the same ``seq_type``.

- ``run_paired_tumours``: Possible values are ``True`` or ``False``. This option determines whether to run paired tumours. Setting this to ``False`` is useful for naturally unpaired or tumour-only analyses (*e.g.* for RNA-seq), which is normally done while setting ``run_paired_tumours_as_unpaired`` to True in case there are any paired tumours.

- ``run_unpaired_tumours_with``: Possible values are ``None``, ``"unmatched_normal"``, or ``"no_normal"``. This option determines what to pair with unpaired tumours. Specifying ``None`` means that unpaired tumours will be skipped for the given module. This option cannot be set to ``None`` if ``run_paired_tumours_as_unpaired`` is ``True``. Specifying ``"unmatched_normal"`` means that unpaired tumours will be run by being paired with the unmatched normal sample given by ``unmatched_normal_id`` (see below). Specifying ``"no_normal"`` means that unpaired tumours will be run without a normal sample. Note that modules need to be specifically configured to be run in paired and/or unpaired mode, since the commands of the underlying tools probably need to be tailored accordingly.

- ``run_paired_tumours_as_unpaired``: Possible values are ``True`` or ``False``. This option determines whether paired tumours should be run as unpaired (*i.e.* separate from their matched normal sample). This is useful for benchmarking purposes or preventing unwanted paired analyses (*e.g.* in RNA-seq analyses intended to be tumour-only).

.. _module-accessory-files-and-scripts:

Module Accessory Files and Scripts
----------------------------------

When you create a new module from the cookiecutter template, you will notice an empty ``etc/`` subdirectory in your module directory. This folder is meant to contain any additional file required to run your module. For example, the ``manta`` module requires configuration files, which are stored in ``etc/``. Scripts can also be stored in this directory. That said, if a script is generally useful, you might want to submit a pull request to the `lcr-scripts repository`_. The purpose of this separate repository is to avoid storing useful scripts in a nested directory within the `lcr-modules repository`_. Just as with lcr-modules, lcr-scripts are versioned and come with conda environment YAML files. For your convenience, there is a ``SCRIPTSDIR`` directory placeholder you can use in your default configuration file. 

You can look at how the ``augment_manta_vcf.py`` script from lcr-scripts is used in the ``manta`` module version 2.0. 

Advanced Module Features
========================

.. _required-sample-metadata:

Required Sample Metadata
------------------------

Every module requires the samples table, which contains metadata on the samples being analyzed. The minimum set of columns expected by ``lcr-modules`` are the ``sample_id``, ``patient_id``, ``seq_type``, and ``tissue_status`` columns (see :ref:`required-columns` for more info). These requirements are spelled out using schemas in YAML format. The base requirements can be found in ``schemas/base/base-1.0.yaml``.

Some modules will need additional metadata (*e.g.* the strandedness of RNA-seq libraries). These extra requirements should also be described in schema files. To promote modularity, each required column should have its own file to promote modularity. An exception can be made for a set of columns should always be present together. The new schemas should be stored in the shared ``schemas/`` directory and then symlinked into individual modules. Symlinks are used to keep the repository lightweight and promote reuse of schemas between modules.

An example single-column schema file can be found in ``schemas/ffpe_status/ffpe_status-1.0.yaml``, where as a multi-column schema file should look like the base schema, *i.e.* ``schemas/base/base-1.0.yaml``.

**Important:** Read the section below on :ref:`conditional-module-behaviour-dev` for an explanation on why you should avoid adding new wildcards beyond the standard ones described in :ref:`other-rules`.

.. _conditional-module-behaviour-dev:

Conditional Module Behaviour
----------------------------

One size doesn’t always fit all, so modules sometimes have to tailor their behaviour based on sample attributes. Snakemake offers more than one avenue to implement these conditional behaviours. The simplest approach is to create parallel rules, which will handle samples differently based on the file names, potentially using wildcard constraints. However, this approach has two major issues.

First, the resulting parallel rules are mostly identical except for a few, often minor differences (*e.g.* a single command-line argument). This redundancy violates the `DRY Principle`_, making the module harder to maintain and more vulnerable to bugs. This pitfall can be avoided by merging the two rules and using the :ref:`switch-on-wildcard-value` function from :py:mod:`oncopipe` described below.

Second, it requires the module developer to encode the sample attributes in the file names. While this is not a severe limitation on its own, it complicates the task of connecting modules together because the file names in downstream modules will need to include every wildcard from upstream modules. This would not only lead to unsustainably long file names, but the file names of a module shouldn’t depend on which modules are upstream to ensure modularity. The accumulation of module-specific wildcards can be avoided using the :ref:`switch-on-sample-metadata` function from :py:mod:`oncopipe` described below.

   To give a specific example, let’s say the ``salmon`` module requires the strandedness of the RNA-seq samples, so this information is encoded in the file name, *e.g.* ``{sample_id}.{strandedness}.quant``. Once we have quantified gene expression in all RNA-seq samples, we wish to perform cohort-wide correction for library size. Unfortunately, we need to pull the information about strandedness from the sample metadata in order to find the ``salmon`` output files because it’s part of the file names, even though that information isn’t relevant to our library size correction module.

**Important:** The :py:func:`oncopipe.switch_on_wildcard` and :py:func:`oncopipe.switch_on_column` functions do not currently support :ref:`directory-placeholders-dev`. This `issue <https://github.com/LCR-BCCRC/lcr-modules/issues/27>`__ will track the implementation.

.. _switch-on-wildcard-value:

Switch on Wildcard Value
~~~~~~~~~~~~~~~~~~~~~~~~

You can use the :py:func:`oncopipe.switch_on_wildcard` function to dynamically set the value of an input file or parameter for a snakemake rule based on the value of a wildcard. The first argument (``wildcard``) is the name of the wildcard, and the second argument (``options``) is a dictionary mapping possible values for the wildcard to the corresponding values that should be returned.

This dictionary can make use of special keys. The most important one to note is the ``"_default"`` special key, whose associated value is selected if the wildcard value isn’t among the other keys. You should check out :py:func:`oncopipe.switch_on_wildcard` to find out about the other special keys. (:ref:`faq-underscore`)

By default, the :py:func:`oncopipe.switch_on_wildcard` will replace any placeholders (using the same format as the ``shell`` directive; *e.g.* ``{wildcards.seq_type}``) with the actual values. This beheviour can be tweaked with the ``format`` (default = ``True``) and ``strict`` (default = ``False``) optional arguments. See the function docstring for more information on these optional arguments.

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

``CFG["switches"]["manta_config"]`` contains the dictionary representation of the ``manta_config`` section from the YAML file shown above. You can see how the ``"_default"`` special key is being used here (see :ref:`switch-on-wildcard-value` for more info) as well as the ``{MODSDIR}`` placeholder for the module subdirectory (see :ref:`directory-placeholders-dev` for more info).

.. code:: python

   # This is the dictionary stored in `CFG["switches"]["manta_config"]`
   {
       '_default': '{MODSDIR}/etc/manta_config.default.ini',
       'mrna': '{MODSDIR}/etc/manta_config.high_sensitivity.ini',
       'capture': '{MODSDIR}/etc/manta_config.high_sensitivity.ini'
   }

.. _switch-on-sample-metadata:

Switch on Sample Metadata
~~~~~~~~~~~~~~~~~~~~~~~~~

As I mentioned in :ref:`conditional-module-behaviour-dev`, adding wildcards for conditional behaviour in a Snakefile is unsustainable and goes against the core principle of modularity. One workaround is to query the metadata for each sample (or each tumour-normal pair) and to update the tool command accordingly. The approach is similar to a :ref:`switch-on-wildcard-value`, but with a few notable differences.

The function to use is :py:func:`oncopipe.switch_on_column(` where the first argument (``column``) is the column name, the second argument (``samples``) is the samples data frame (typically ``CFG["samples"]``), and the third argument (``options``) is a dictionary mapping possible values in the column to the corresponding values that should be returned. This dictionary follows the same structure as the :ref:`switch-on-wildcard-value`. An additional albeit optional argument is called ``match_on``, which needs to be set to either ``"tumour"`` (default) or ``"normal"`` to determine whether the function uses the ``wildcards.tumour_id`` or ``wildcards.normal_id`` to look up a sample ID. The function will automatically use ``wildcards.seq_type`` to also filter on sequencing data type.

   At the moment, this function only works for tumour-based modules (*e.g.* paired variant calling). It should soon be generalized to also work with sample-based modules (*e.g.* STAR alignment). This issue is tracked `here <https://github.com/LCR-BCCRC/lcr-modules/issues/35>`__.

The code block below shows how we could achieve the same outcome using :py:func:`oncopipe.switch_on_column` for the example given in :ref:`switch-on-wildcard-value`. The only difference other than the function name is the addition of the ``samples`` argument before providing the same ``options`` dictionary. By default, the function will use ``wildcards.tumour_id`` (and ``wildcards.seq_type``) to look up the sample in ``CFG["samples"]``. In practice, you would simply use :py:func:`oncopipe.switch_on_wildcard` since ``seq_type`` is available as a wildcard.

.. code:: python

   rule _manta_configure:
       input:
           config = op.switch_on_column("seq_type", CFG["samples"], CFG["switches"]["manta_config"])

Switch on File Contents
-----------------------

The behaviour of some module depends on the contents (or existence) of input or intermediate files. The best way to address this is using `Snakemake Checkpoints`_. They are a bit complicated to implement, but you can look at the ``manta`` module (version 1.0) for an example. Do note that checkpoints can be slow because the function using the checkpoint is run sequentially for each sample.
