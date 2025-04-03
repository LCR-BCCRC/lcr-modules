# Changelog

All notable changes to the `hotmaps` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-12-13

This release was authored by Manuela Cruz.

- HotMAPS maps missense mutations to protein residues on PDB structures. Prior to running, PDB structures must be downloaded and their paths must be specified in the config under `["lcr-modules"]["hotmaps"]["pdb_structure_dirs"]`. More information on downloading PDB structures here: https://github.com/KarchinLab/HotMAPS_2016/wiki/Tutorial-(Exome-scale) <br />
Please note that the `rsync` command provided at the link above will download <ins>all</ins> PDB structures, which will result in a huge directory. You can download only the necessary structures (`biounit` and `structures` PDBs) using these paths:<br />
`rsync.rcsb.org::ftp_data/biounit/coordinates/all/`<br />
`rsync.rcsb.org::ftp_data/structures/all/pdb/`<br />
More information here:<br />
https://www.wwpdb.org/ftp/pdb-ftp-sites <br />
Take a look at the directory structure of the wwPDB database here:<br />
https://files.wwpdb.org/pub/pdb/data/ <br />
Finally, the example config provided in the HotMAPS repo can help if you are confused about which PDB structures are required and how to specify their paths in the config file: https://github.com/KarchinLab/HotMAPS_2016/blob/master/config.txt. The HotMAPS scripts have been updated such that PDB structures from either `ftp_data/structures/all/pdb/` and `ftp_data/structures/divided/pdb/` are accepted.

- HotMAPS depends on a MySQL database created by the Karchin Lab that maps genomic coordinates to PDB residues. This database must be downloaded before running HotMAPS and the database name, host, user, and password information must be specified in the config. More information here: https://github.com/KarchinLab/HotMAPS_2016/wiki/MySQL-database

- The HotMAPS module relies on a grch37 "master maf" file that contains the mutation data for all samples, and an accompanying "sample_set" file that maps unique sample IDs listed in the first column to sample subsets using a binary approach.

- Multiple input "master maf" files can be specified and will be combined into one file.

- Paths to variant blacklists can be provided in the config. They must contain a column `chrpos` that contains mutations in `{chromosome}:{position}` format. Chromosomes should not include the "chr" prefix.

- HotMAPS can improperly map DNPs if they span two residues. The module includes rules that separate these DNPs from the analysis and reannotates them as two SNPs.

- HotMAPS reformats input mutations by loading them into the `mupit_modbase` MySQL database and then retrieving them into a file named `mutations.txt`. Upon loading of the mutations into the database, pre-existing mutations in that database are deleted. This prevents the `mutations.txt` file from including mutations from other sample subsets and prevents the file from growing too large as more subsets are run. Only one job for this rule should occur at any given time in order to ensure that mutations from one subset are not deleted if other subsets are simultaneously loading their mutations. To ensure this, the `mysql` resource for  `rule _hotmaps_get_mutations` is set to 1 in the config.

- If more than one sample subset is being run, ensure that only one database loading job of `rule _hotmaps_get_mutations` can run at a time by specifying `--resources mysql=1` on the command line.

## [1.1] - 2024-08-22

This update was authored by Sierra Gillis.

- Added an output directory for the prepare step
- Was tested with an update made to `lcr-scripts/generate_smg_inputs` that uses threads to read in maf files faster, and only reads in relevant columns