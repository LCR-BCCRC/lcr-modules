#!/usr/bin/python3

# Load modules
import re

# Get cookiecutter variables
module_name = "{{ cookiecutter.module_name }}"
module_author = "{{ cookiecutter.module_author }}"
original_author = "{{ cookiecutter.original_author }}"
input_file_type = "{{ cookiecutter.input_file_type }}"
output_file_type = "{{ cookiecutter.input_file_type }}"

# Ensure that module name conforms to requirements
assert re.fullmatch(r"[a-z0-9_]+", module_name), (
    f"`module_name` ('{module_name}') should only consist of lowercase "
    "alphanumerical characters or underscores (_i.e._ no spaces)."
)

# Ensure that output file type conforms to requirements
assert re.fullmatch(r"[a-z0-9_]+", output_file_type), (
    f"`output_file_type` ('{output_file_type}') should only consist of lowercase "
    "alphanumerical characters or underscores (_i.e._ no spaces)."
)
