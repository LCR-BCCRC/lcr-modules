# Lcr-modules continuous integration (CI) via Github actions

## About
This document briefly outlines the automatic unit tests configured for lcr-modules. 
This should reduce the likelihood that a module is merged into master in a broken state, 
and improve user experience by ensuring that a functional test setup is always available.

## How it works
The continuous integration (CI) tests are implemented via Github Actions, and are triggered every time a push is made to the lcr-modules repository. 
These tests are executed by a “runner” (i.e. a computer/server), which clones the repository and associated dependencies, sets up the environment, 
and perfoms the specified tests. These tests include:
1)	Execute a dry run of the capture demo workflow
2)	Run the capture demo workflow

To minimize runtime, a synthetic Tumour-Normal pair is used representing ultra-deep targeted sequencing.


NOTE: The pathseq module is currently not tested, as the bacterial database takes an exceptionally long time to prepare. 
Also note that the test environment state is NOT cleaned between runs (i.e. python, snakemake, and oncopipe installs will persist, 
and packages used via snakemake conda environments will be cached)

## Checking the test results

The result of all tests can be found under Actions > lcr-modules on Github: 
https://github.com/LCR-BCCRC/lcr-modules/actions/workflows/lcr-modules-test-actions.yml. 
Click on a given test for more details about what was executed, and what steps failed.

## Configuring the tests
### Configuring the runner
The CI tests are executed by a self-hosted runner (see https://docs.github.com/en/actions/hosting-your-own-runners/about-self-hosted-runners). 
This is both free and allows for improved flexibility compared to cloud-based runners. However, due to security implications, it is highly recommended
that the runner be installed inside a secure docker container. Conda and python should also be installed in this container, due to MANY painful 
limitations of Github Actions

Currently, the runner is hosted on Thanos (thanos.mbb.sfu.ca), in the container `thanos-githubrunner`. To prevent the `lcr-modules-reference`
directory from being re-generated during each run, it should be mounted to the container. For instance:
`sudo docker run -e RUNNER_ALLOW_RUNASROOT="1" --mount type=bind,source=/reference/lcr-modules-reference/,target=/reference/lcr-modules-reference/ thanos-githubrunner:1.5 /root/actions-runner/run.sh`

To configure the runner, run `actions-runner/config.sh` (this only needs to be done once), then start the runner using `actions-runner/run.sh`. 
This has already been performed for the docker container.

### The configuration
The steps the runner executes are configured via Github Actions, which is extensively documented (https://docs.github.com/en/actions). 
The corresponding actions file for these tests can be found under `.github/workflows/lcr-modules-test-actions.yml`. 

Note the ”runs-on” tags should match those assigned to the runner during runner configuration.

### Demo data
Demo data, including sample BAM files and capture regions, can be found in https://github.com/LCR-BCCRC/lcr-modules-demo. The sample config and snakemake
config for these tests files are found under `demo/capture_config_ci.yaml` and `demo/data/samples_ci.tsv`.

### The results
All tests are executed under `~/actions-runner/runner-work/lcr-modules/lcr-modules/`, which includes clones of all three repositories 
(lcr-modules, lcr-scripts, and lcr-modules-demo). If the tests fails, the corresponding logs can be found under 
`/home/github-runner/actions-runner/runner-work/lcr-modules/lcr-modules/lcr-modules/demo/.snakemake`. If the tests are successful, the `results` folder is
cleared.
