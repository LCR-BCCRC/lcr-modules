#!/usr/bin/python3

# Load modules
import re

# Get cookiecutter variables
module_name = "{{ cookiecutter.module_name }}"
module_author = "{{ cookiecutter.module_author }}"
original_author = "{{ cookiecutter.original_author }}"
input_file_type = "{{ cookiecutter.input_file_type }}"
output_file_type = "{{ cookiecutter.input_file_type }}"
module_run_per = "{{ cookiecutter.module_run_per }}"

# Get seq_type options
seq_types = dict()
{%- for seq_type, mode in cookiecutter.items() if seq_type.startswith("seq_type.") %}
seq_types["{{ seq_type }}".replace("seq_type.", "", 1)] = "{{ mode }}"
{%- endfor %}

# Ensure that module name conforms to requirements
assert re.fullmatch(r"[a-z0-9_]+", module_name), (
    f"`module_name` ('{module_name}') should only consist of lowercase "
    "alphanumerical characters or underscores (_i.e._ no spaces)."
)

# Ensure that input file type conforms to requirements
assert re.fullmatch(r"[a-z0-9_]+", input_file_type), (
    f"`input_file_type` ('{input_file_type}') should only consist of lowercase "
    "alphanumerical characters or underscores (_i.e._ no spaces)."
)

# Ensure that output file type conforms to requirements
assert re.fullmatch(r"[a-z0-9_]+", output_file_type), (
    f"`output_file_type` ('{output_file_type}') should only consist of lowercase "
    "alphanumerical characters or underscores (_i.e._ no spaces)."
)

# Ensure that a sample-based module isn't run with a paired seq_type
if module_run_per == "sample":
    for seq_type, mode in seq_types.items():
        assert mode not in ["matched_only", "allow_unmatched"], (
            f"`module_run_per` set to 'sample', but `seq_type.{seq_type}` set to "
            f"'{mode}'. For paired analyses, select 'tumour' for `module_run_per`."
            f"For unpaired analyses, select 'unpaired' for `seq_type.{seq_type}`."
        )
