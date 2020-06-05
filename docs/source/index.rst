.. include:: links.rst

Welcome to lcr-modules's Documentation!
=======================================

Getting Started With lcr-modules
--------------------------------

- **Users:** Check out this :ref:`getting-started-user` guide and the `Demo Project`_.

- **Contributors:** Check out this :ref:`getting-started-dev` guide and the `lcr-modules repository`_.

Motivation
----------

This project aims to become a collection of standard analytical modules for genomic and transcriptomic data. Too often do we copy-paste from each other’s pipelines, which has several pitfalls:

.. code::

   * Too much time spent on routine analyses           * Increased risk for hidden logical bugs
   * Duplicated effort within and between labs         * No consistently used pipelining tool
   * Inefficient dissemination of best practices       * Steep learning curve for new members

Fortunately, all of these problems can be solved with standardized analytical modules, and the benefits are many:

.. code::

   * Projects can ramp up faster                       * Consistent intermediate/output files
   * Streamline efforts between labs                   * More reproducible analyses
   * Define analytical best practices                  * Easier-to-write methods
   * Consolidate collective expertise                  * Automated logging and “paper trail”
   * Simplify member onboarding                        * Easier peer review of code

                              * And happier bioinformaticians!

.. _what-are-modules:

What Are Modules?
-----------------

Each module accomplishes a specific analysis, generally centered around a specific tool (*e.g.* Strelka2, Manta, MutSigCV). Analyses—and by extension, modules—can be organized into different levels. The figure below contains for examples for each level.

- **Level-1 Analyses:** They process raw sequencing data, generally producing BAM/FASTQ files.
- **Level-2 Analyses:** They perform sample-level analyses on level-1 output, such as variant calling and gene expression quantification. 
- **Level-3 Analyses:** They aggregate sample-specific level-2 output and perform cohort-wide analyses, such as the identification of sifgnificantly mutated genes.
- **Level-4 Analyses:** They are project-specific and are meant to ask specific questions of the data. These are the analyses you ideally want to spend your time on. 

.. figure:: ../../images/module_levels.png
   :alt: Module Levels
   :target: _images/module_levels.png

   Module Levels


.. toctree::
   :maxdepth: 4
   :hidden:
   :caption: For Users
   
   for_users.rst


.. toctree::
   :maxdepth: 4
   :hidden:
   :caption: For Developers
   
   for_developers.rst


.. toctree::
   :maxdepth: 4
   :hidden:
   :glob:
   :caption: Oncopipe
   
   oncopipe/*


.. toctree::
   :maxdepth: 4
   :hidden:
   :caption: Frequently Asked Questions
   
   faq.rst
