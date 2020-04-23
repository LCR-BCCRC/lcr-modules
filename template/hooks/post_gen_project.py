#!/usr/bin/python3

# Load modules
import os

# Symlink schema to avoid duplicating files
os.symlink("../../../../schemas/base.yaml", "1.0/schemas/base.yaml")

# Symlink samtools conda environment to avoid duplicating files
os.symlink("../../../../envs/samtools/samtools-1.9.yaml", "1.0/envs/samtools-1.9.yaml")

# Delete hidden files serving only to have the directories be visible to Git
os.remove("1.0/envs/.gitkeep")
os.remove("1.0/etc/.gitkeep")
os.remove("1.0/schemas/.gitkeep")
