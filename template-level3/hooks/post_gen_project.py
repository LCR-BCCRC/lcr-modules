#!/usr/bin/python3

# Load modules
import os

# Symlink schema to avoid duplicating files
os.symlink("../../../../schemas/base/base-1.0.yaml", "1.0/schemas/base-1.0.yaml")

# Delete hidden files serving only to have the directories be visible to Git
os.remove("1.0/envs/.gitkeep")
os.remove("1.0/etc/.gitkeep")
os.remove("1.0/schemas/.gitkeep")
