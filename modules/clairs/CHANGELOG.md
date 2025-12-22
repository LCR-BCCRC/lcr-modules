# Changelog

All notable changes to the `clairs` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2025-10-14

This release was authored by Giuliano Banco.

This module is designed to only work with a sample table that includes the columns 'chemistry' and 'platform'. These are required by schemas.

This module only works with hg38, as of version 1.0.

This module uses singularity/apptainer containers. These containers require binding/mounting directories with real version of all files inside the container. You cannot bind symlinks, all binding should be handled by the module.