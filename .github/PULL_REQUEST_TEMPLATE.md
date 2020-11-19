# Pull Request Checklists

**Important:** When opening a pull request, keep only the applicable checklist and delete all other sections.

## Checklist for New Module

### Required

- [ ] I used the cookiecutter template and updated the placeholder rules.

- [ ] The snakemake rules follow the [design guidelines](https://lcr-modules.readthedocs.io/en/latest/for_developers.html#module-rules). 

  - [ ] All references to the `rules` object (_e.g._ for input files) are wrapped with `str()`.

- [ ] Every rule in the module is either listed under `localrules` or has the `threads` and `resources` directives.

- [ ] Input and output files are being symlinked into the `CFG["inputs"]` and `CFG["outputs"]` subdirectories, respectively.

- [ ] I updated the final target rule (`*_all`) to include every output rule.

- [ ] I explained important module design decisions in `CHANGELOG.md`.

- [ ] I tested the module on real data for all supported `seq_type` values.

- [ ] I updated the `default.yaml` configuration file to provide default values for each rule in the module snakefile.

- [ ] I did not set any global wildcard constraints. Any/all wildcard constraints are set on a per-rule basis. 

- [ ] I ensured that all symbolic links are relative and self-contained (_i.e._ do not point outside of the repository).

- [ ] I replaced every value that should (or might need to) be updated in the default configuration file with `__UPDATE__`.

- [ ] I recursively searched for all comments containing `TODO` to ensure none were left. For example:

  ```bash
  grep -r TODO modules/<module_name>/1.0
  ```

### If applicable

- [ ] I added more granular output subdirectories.

- [ ] I added rules to the `reference_files` workflow to generate any new reference files.

- [ ] I added subdirectories with large intermediate files to the list of `scratch_subdirectories` in the `default.yaml` configuration file.

- [ ] I updated the list of available wildcards for the input files in the `default.yaml` configuration file.

## Checklist for Updated Module

To be completed.
