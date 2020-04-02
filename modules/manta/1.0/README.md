# Manta Module

## Software Requirements

_Please read the main `lcr-modules` `README.md` for information on the general software requirements for using this repository._

The required conda environments are available in `envs/`. If conda fails to build the environment (_e.g._ during dependency resolution), you can subset the package dependencies to the core packages listed below. While the environment will not be exactly as intended, using the intended versions for the core packages should be close to perfectly reproducible.

### Core Packages

```plain
tabix
manta
svtools
pyvcf
```
