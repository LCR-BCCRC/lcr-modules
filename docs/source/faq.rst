.. _faq-conda-fail:

How do I handle a conda environment that fails to build?
========================================================

While conda brings us much closer to computational reproducibility, it isn’t perfect. Issues arise when conda packages are removed from `Anaconda Cloud <https://anaconda.org/>`__ or when the dependency resolution algorithm changes. We suggest you try the following steps in order:

1. Remove the build IDs from the conda environment YAML file, although this should already be the case for all environments in ``lcr-modules``.
2. Remove the versions for the offending package(s) (*i.e.* the one(s) mentioned in the error message).
3. Remove the offending packages altogether.
4. Remove the dependency packages, leaving only the “target packages”. This generally means subsetting to the core conda packages listed in a module’s README for the environment in question. While extreme, the hope is that the versions of the dependency packages are not crucial for maintaining scientific reproducibility.
5. Remove the versions for the target packages.
6. If you reach this point, it usually means that a target package is problematic. If possible, replace that package with the same (or similar) version from another Anaconda channel. Ideally, restore the YAML file first and cycle through the previous steps.
7. Install the software tools manually (ideally the versions specified in the YAML file) and ensure they are available in your ``PATH`` environment variable.

.. _faq-underscore:

What does the underscore prefix mean?
=====================================

The underscore prefix is mainly used to avoid name conflicts. This convention is borrowed from Python. For instance, ``collections.namedtuple`` has an ``_asdict()`` method, where the underscore helps prevent clashes with user-defined attributes for the ``namedtuple``. For more examples in Python, check out this `blog post <https://medium.com/python-features/naming-conventions-with-underscores-in-python-791251ac7097>`__.

In ``lcr-modules``, the underscore prefix is used in a few areas. First, the name of every rule or function defined in a module starts with an underscore followed by the module name (*e.g.* ``_manta``). This minimizes the risk for clashing with other rule/function names defined elsewhere by the user, which isn’t allowed by Snakemake. Second, the underscore prefix is used for dictionary keys with special behaviour, such as the ``"_default"`` key in the ``op.switch_on_wildcard()`` `function <#switch-on-wildcard-value>`__. Third, the shared ``lcr-modules`` configuration is stored under the ``_shared`` key, which is done to avoid clashing with a potential module called ``shared``.

.. _faq-symlink:

What is the difference between ``op.relative_symlink()`` and ``os.symlink()``?
==============================================================================

Behind the scenes, ``op.relative_symlink()`` uses ``os.symlink()`` while ensuring that the symlinks are relative and correct regardless of the current working directory. This is equivalent to the ``-r`` option on modern version of the ``ln`` command-line tool.

.. _faq-cfg-nameerror:

Why am I running into a ``NameError: name 'CFG' is not defined`` exception?
===========================================================================

Each module creates a ``CFG`` variable as a convenient but temporary pointer to the module configuration (*i.e.* ``config["lcr-modules"]["<module_name>"]``). Because each module uses this variable name, the ``op.cleanup_module()`` function deletes the variable to be safe. Hence, you will run into this ``NameError`` exception if some code tries to use ``CFG`` after it’s been deleted. If you use ``CFG`` in the rule directives that are evaluated when the module snakefile is parsed (*e.g.* ``input``, ``output``, ``log``, ``params``, etc.), it’s not an issue. However, if you use this variable in a function or ``run`` directive, *i.e.* code that is run after the ``op.cleanup_module()`` function is run, you will get the error above. You can fix this error by adding this line of code before using the ``CFG`` variable, which recreates the variable in a local scope:

.. code:: python

   # Replace <module_name> with the actual module name (e.g., `star`)
   CFG = config["lcr-modules"]["<module_name>"]

If you're using ``CFG`` in an anonymous ``lambda`` function, then you can just use the ``config`` object directly. For example:

.. code:: python

   lambda w: config["lcr-modules"]["<module_name>"][w.genome_build]["some_ref"]

.. _faq-memory-per-thread:

How do I specify the available memory per thread for a command-line tool?
=========================================================================

The ``mem_mb`` resource is meant to represent the total amount of memory used by all threads of a given process. Some tools have command-line arguments allowing the user to specify the amount of memory they can use, such as any Java-based application (*i.e.* using ``-Xmx``). In some cases, the tool expects the amount of memory per thread (*e.g.* ``samtools sort``), whereas ``resources.mem_mb`` represents the total amount of memory. `Arithmetic expansion <https://www.shell-tips.com/bash/performing-math-calculation-in-bash#using-arithmetic-expansion-with-or>`__ in Bash allows you to circumvent this issue as long as you are dealing with integers, which should be the case with ``threads`` and ``mem_mb``. For example, here’s how you would divide two integers and print the result: ``echo $((12000 / 12))``. We can leverage the same syntax within the ``shell`` directive of a Snakemake rule. The example below is taken from the ``samtools sort`` rule in the ``utils`` module.

.. code:: bash

   rule:
       ...
       shell:
           op.as_one_line("""
           samtools sort {params.opts} -@ {threads} -m $(({resources.mem_mb} / {threads}))M
           -T {params.prefix} -o {output.bam} {input.bam} > {log.stdout} 2> {log.stderr}
           """)
