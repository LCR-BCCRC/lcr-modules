#!/usr/bin/env python3
"""Parse an LCR-modules runtime_config.yaml and print Snakemake flag assignments.

Usage: eval "$(python3 parse_runtime_config.py [runtime_config.yaml])"

Prints shell variable assignments suitable for eval in bash:
  SNAKEMAKE_CONDA_PREFIX   -- value for --conda-prefix, empty if not configured
  SNAKEMAKE_CONTAINER_FLAG -- --use-apptainer or --use-singularity (version-aware), empty if no apptainer section
  SNAKEMAKE_PREFIX_FLAG    -- --apptainer-prefix or --singularity-prefix (version-aware)
  SNAKEMAKE_ARGS_FLAG      -- --apptainer-args or --singularity-args (version-aware)
  SNAKEMAKE_SIF_PREFIX     -- value for the prefix flag, empty if not configured
  SNAKEMAKE_BIND_PATHS     -- value for the args flag (e.g. "-B /p1 -B /p2"), empty if not configured
"""

import shlex
import sys


def detect_container_flags():
    """Return (use_flag, prefix_flag, args_flag) appropriate for the installed Snakemake version."""
    try:
        import subprocess
        help_output = subprocess.check_output(
            ["snakemake", "--help"], stderr=subprocess.STDOUT, timeout=10
        ).decode()
        if "--use-apptainer" in help_output:
            return "--use-apptainer", "--apptainer-prefix", "--apptainer-args"
    except Exception:
        pass
    return "--use-singularity", "--singularity-prefix", "--singularity-args"


def parse(config_path):
    import yaml
    with open(config_path) as f:
        c = yaml.safe_load(f) or {}

    conda_prefix = (c.get("conda") or {}).get("prefix") or ""

    container_flag = ""
    prefix_flag = ""
    args_flag = ""
    sif_prefix = ""
    bind_paths = ""

    apptainer = c.get("apptainer")
    if apptainer:
        container_flag, prefix_flag, args_flag = detect_container_flags()
        sif_prefix = str(apptainer.get("prefix") or "")
        raw_paths = [str(p) for p in (apptainer.get("bind_paths") or []) if p]
        bind_args = " ".join("-B " + p for p in raw_paths)
        extra_args = str(apptainer.get("extra_args") or "").strip()
        bind_paths = " ".join(filter(None, [bind_args, extra_args]))

    return conda_prefix, container_flag, prefix_flag, args_flag, sif_prefix, bind_paths


def emit_empty():
    for var in ("SNAKEMAKE_CONDA_PREFIX", "SNAKEMAKE_CONTAINER_FLAG",
                "SNAKEMAKE_PREFIX_FLAG", "SNAKEMAKE_ARGS_FLAG",
                "SNAKEMAKE_SIF_PREFIX", "SNAKEMAKE_BIND_PATHS"):
        print(f"{var}=''")


def main():
    if len(sys.argv) < 2 or not sys.argv[1]:
        emit_empty()
        return

    try:
        conda_prefix, container_flag, prefix_flag, args_flag, sif_prefix, bind_paths = parse(sys.argv[1])
    except Exception as e:
        print(f"Error reading {sys.argv[1]}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"SNAKEMAKE_CONDA_PREFIX={shlex.quote(conda_prefix)}")
    print(f"SNAKEMAKE_CONTAINER_FLAG={shlex.quote(container_flag)}")
    print(f"SNAKEMAKE_PREFIX_FLAG={shlex.quote(prefix_flag)}")
    print(f"SNAKEMAKE_ARGS_FLAG={shlex.quote(args_flag)}")
    print(f"SNAKEMAKE_SIF_PREFIX={shlex.quote(sif_prefix)}")
    print(f"SNAKEMAKE_BIND_PATHS={shlex.quote(bind_paths)}")


if __name__ == "__main__":
    main()
