#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import os

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["ega_download"]`
CFG = op.setup_module(
    name = "ega_download",
    version = "2.0",
    subdirectories = ["inputs", "ega_download", "decrypt", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _ega_output_files,
    _ega_all


##### RULES #####


# Download the crypt4gh-encrypted file from the EGA Live Outbox via SFTP.
# SFTP authentication uses the user's registered SSH public key (id_rsa.pub or
# id_ed25519.pub), which must be shared with EGA separately from the crypt4gh key.
# Files are stored as {egad}/{file_name}.c4gh on the outbox SFTP server.
rule _ega_download_file:
    output:
        encrypted = CFG["dirs"]["ega_download"] + "{seq_type}/{file_name}.c4gh"
    log:
        stdout = CFG["logs"]["ega_download"] + "{seq_type}/{file_name}_download.stdout.log",
        stderr = CFG["logs"]["ega_download"] + "{seq_type}/{file_name}_download.stderr.log"
    params:
        sftp_server = CFG["sftp_server"],
        remote_path = lambda wc: (
            f"{CFG['sftp_outbox_path']}/{CFG['egad']}/{wc.file_name}.c4gh"
            if CFG["sftp_outbox_path"]
            else f"{CFG['egad']}/{wc.file_name}.c4gh"
        )
    threads:
        CFG["threads"]["ega_file_download"]
    resources:
        **CFG["resources"]["ega_file_download"]
    shell:
        op.as_one_line("""
        printf 'get {params.remote_path} {output.encrypted}\\n'
        | sftp -b - {params.sftp_server}
        > {log.stdout}
        2> {log.stderr}
        """)


# Decrypt the crypt4gh-encrypted file using the user's crypt4gh private key.
rule _ega_decrypt_file:
    input:
        encrypted = CFG["dirs"]["ega_download"] + "{seq_type}/{file_name}.c4gh"
    output:
        decrypted = CFG["dirs"]["decrypt"] + "{seq_type}/{file_name}"
    log:
        stderr = CFG["logs"]["decrypt"] + "{seq_type}/{file_name}_decrypt.stderr.log"
    params:
        private_key = CFG["private_key"]
    conda:
        CFG["conda_envs"]["crypt4gh"]
    threads:
        CFG["threads"]["ega_file_decrypt"]
    resources:
        **CFG["resources"]["ega_file_decrypt"]
    shell:
        op.as_one_line("""
        crypt4gh decrypt
        --sk {params.private_key}
        < {input.encrypted}
        > {output.decrypted}
        2> {log.stderr}
        """)


# Symlink decrypted files into 99-outputs/ preserving the EGA file name
rule _ega_output_files:
    input:
        ega_file = CFG["dirs"]["decrypt"] + "{seq_type}/{file_name}"
    output:
        ega_file = CFG["dirs"]["outputs"] + "{seq_type}/{file_name}"
    run:
        op.relative_symlink(input.ega_file, output.ega_file, in_module=True)


# Generates the target files for each run
rule _ega_all:
    input:
        expand(
            [str(rules._ega_output_files.output.ega_file)],
            zip,
            seq_type  = CFG["samples"]["seq_type"],
            file_name = CFG["samples"]["file_name"],
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
