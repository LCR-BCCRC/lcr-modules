# MODULES

import os
import copy
import inspect
import functools
import itertools
import subprocess
import collections.abc
from datetime import datetime
from collections import defaultdict, namedtuple

import yaml
import pandas as pd
import snakemake as smk
from snakemake.logging import logger


# CONSTANTS

DEFAULT_PAIRING_CONFIG = {
    "genome": {
        "run_unpaired_tumours_with": "unmatched_normal",
        "run_paired_tumours": True,
        "run_paired_tumours_as_unpaired": False,
    },
    "capture": {
        "run_unpaired_tumours_with": "unmatched_normal",
        "run_paired_tumours": True,
        "run_paired_tumours_as_unpaired": False,
    },
    "mrna": {
        "run_unpaired_tumours_with": "no_normal",
        "run_paired_tumours": False,
        "run_paired_tumours_as_unpaired": True,
    },
    "mirna": {
        "run_unpaired_tumours_with": "no_normal",
        "run_paired_tumours": False,
        "run_paired_tumours_as_unpaired": True,
    },
}

DOCS = {
    "update_config": "https://lcr-modules.readthedocs.io/en/latest/for_users.html#updating-configuration-values"
}


# SESSION


class _Session:
    """Session for storing Snakemake config internally."""

    def __init__(self):
        timestamp_key = "LCR_MODULES_LAUNCH_TIMESTAMP"
        if timestamp_key not in os.environ:
            log_path = logger.get_logfile()
            if log_path is None:
                logger.warning(
                    "Warning: The oncopipe package was imported outside of a "
                    "snakefile. Most functions are designed to work within a "
                    "snakefile. Some unexpected behaviour/errors might occur."
                )
                now = datetime.now()
                log_time_fmt = now.strftime("launched-%Y-%m-%d-at-%H-%M-%S")
            else:
                log_file = os.path.basename(log_path)
                log_time = datetime.strptime(
                    log_file, "%Y-%m-%dT%H%M%S.%f.snakemake.log"
                )
                log_time_fmt = log_time.strftime("launched-%Y-%m-%d-at-%H-%M-%S")
            os.environ[timestamp_key] = log_time_fmt
        self.config = None
        self.launched_fmt = os.environ[timestamp_key]

    def set_config(self, config):
        self.config = config


_session = _Session()


def enable_set_functions(config):
    """Enable the `set_*` oncopipe convenience functions.

    Parameters
    ----------
    config : dict
        The Snakemake configuration nested dictionary.
    """
    _session.set_config(config)


# CONVENIENCE FUNCTIONS


def set_input(module, name, value):
    """Use given value for an input file in a module.

    Parameters
    ----------
    module : str
        The module name.
    name : str
        The name of input file field. This is usually taken from the
        module's configuration YAML file.
    value : str or function
        The value to provide for the named input file. In most cases,
        this value will be a plain string, but you can also provide
        an input file function as per the Snakemake documentation,
        where the function would return strings. In all cases, the
        strings can make use of the wildcards that are usually listed
        in the configuration file.
    """
    config = _session.config
    new_value = {"lcr-modules": {module: {"inputs": {name: value}}}}
    smk.utils.update_config(config, new_value)


def set_samples(module, *samples):
    """Use given samples for a module.

    Parameters
    ----------
    module : str
        The module name. This can also be ``"_shared"`` for a value
        that should be inherited by all modules.
    *samples : list of pandas.DataFrame
        One or more pandas data frames that will be concatenated
        before being used by the module. These data frames should
        contain sample tables as described in the documentation.
    """
    config = _session.config
    samples_concat = pd.concat(samples)
    new_value = {"lcr-modules": {module: {"samples": samples_concat}}}
    smk.utils.update_config(config, new_value)


def set_value(value, *keys):
    """Update lcr-modules configuration using simpler syntax.

    This function will automatically create dictionaries if
    accessing a key that doesn't exist and notify the user.

    Parameters
    ----------
    value : anything
        The value to be set at the location specified by ``*keys``.
    *keys : list of str
        All subsequent arguments will be collected into a list of
        strings, which specify the location where to set ``value``.
        You do not need to include the ``"lcr-modules"`` key; it is
        assumed that you are accessing keys therein.
    """
    keys = ["lcr-modules"] + keys
    current_location = _session.config
    for key in keys:
        pass
    # ...
    smk.utils.update_config(current_location, value)


# UTILITIES


def relative_symlink(src, dest, overwrite=True):
    """Creates a relative symlink from any working directory.

    Parameters
    ----------
    src : str
        The source file or directory path.
    dest : str
        The destination file path. This can also be a destination
        directory, and the destination symlink name will be identical
        to the source file name (unless directory).
    overwrite : boolean
        Whether to overwrite the destination file if it exists.
    """

    # Coerce length-1 NamedList instances to strings
    def coerce_namedlist_to_string(obj):
        if isinstance(obj, smk.io.Namedlist) and len(obj) == 1:
            obj = str(obj)
        elif isinstance(obj, smk.io.Namedlist) and len(obj) != 1:
            raise AssertionError(
                f"Got the following Namedlist: {obj!r}. This function "
                "only supports Namedlists of length 1."
            )
        return obj

    src = coerce_namedlist_to_string(src)
    dest = coerce_namedlist_to_string(dest)

    # Check whether arguments are strings
    assert isinstance(src, str), "Source file path must be a string."
    assert isinstance(dest, str), "Destination file path must be a string."

    # Here, you're symlinking a file into a directory (same name)
    dest = dest.rstrip(os.path.sep)
    if not os.path.isdir(src) and os.path.isdir(dest):
        dest_dir = dest
        dest_file = os.path.split(src)[1]
    # Here, you can symlinking a file to a specific location
    # Or you are symlinking a directory to a specific location
    else:
        dest_dir, dest_file = os.path.split(dest)
    os.makedirs(dest_dir, exist_ok=True)

    dest = os.path.join(dest_dir, dest_file)
    if os.path.lexists(dest) and os.path.islink(dest):
        if os.path.realpath(src) == os.path.realpath(dest):
            return
        elif overwrite:
            os.remove(dest)
    assert not os.path.exists(dest), (
        "Symbolic link already exists but points elsewhere: \n"
        f"    Current: {dest} -> {os.path.realpath(dest)} \n"
        f"    Attempted: {dest} -> {os.path.realpath(src)}"
    )

    # Make `src` relative to destination parent directory
    if not os.path.isabs(src):
        dest_dir = os.path.realpath(dest_dir)
        src = os.path.relpath(src, dest_dir)
    os.symlink(src, dest)


def get_from_dict(dictionary, list_of_keys):
    """Access nested index/key in dictionary."""
    return functools.reduce(dict.get, list_of_keys, dictionary)


def as_one_line(text):
    """Collapses a triple-quoted string to one line.

    Line endings do not need to be escaped like in a shell script.
    Spaces and tabs are stripped from each side of each line to
    remove the indentation included in triple-quoted strings.

    This function is useful for long shell commands in a Snakefile,
    especially if it contains quotes that would need to be escaped
    (e.g., in an awk command).

    Returns
    -------
    str
        A single line (i.e., without line endings) of text.
    """
    lines = text.strip().split("\n")
    lines_dedented = [x.strip(" \t") for x in lines]
    return " ".join(lines_dedented)


def list_files(directory, file_ext):
    """Searches directory for all files with given extension.

    The search is performed recursively. The function first tries
    to use the faster `find` UNIX tool before falling back on a
    slower Python implementation.

    Parameters
    ----------
    directory : str
        The directory to search in.
    file_ext : str
        The file extension (excluding the period).

    Returns
    -------
    list of str
        The list of matching files.
    """

    files_all = []

    try:
        # Get list of BAM files using `find` UNIX tool (fast)
        command = ["find", directory, "-name", f"*.{file_ext}"]
        files = subprocess.check_output(command, text=True)
        files_split = files.rstrip("\n").split("\n")
        files_all.extend(files_split)
    except subprocess.CalledProcessError:
        # Slower fallback in Python
        for root, _subdirs, files in os.walk(directory):
            files = [f for f in files if f.endswith(f".{file_ext}")]
            files = [os.path.join(root, f) for f in files]
            files_all.extend(files)

    return files_all


# SNAKEMAKE INPUT/PARAM FUNCTIONS


def create_formatter(wildcards, input, output, threads, resources, strict):
    """Create formatter function based on rule variables."""

    variables = {
        "wildcards": wildcards,
        "input": input,
        "output": output,
        "threads": threads,
        "resources": resources,
    }

    def check_for_clashes(text, wildcards, strict):
        """Check for clashes between wildcard names and rule variables."""
        reserved = ["wildcards", "input", "output", "threads", "resources"]
        potential_clashes = [k for k in wildcards.keys() if k in reserved]
        actual_clashes = [c for c in potential_clashes if "{" + c + "}" in text]
        if len(actual_clashes) > 0 and strict is False:
            raise ValueError(
                f"Some wildcards {actual_clashes} cannot be unambiguously "
                "resolved. Consider setting `strict_mode` to True and "
                f"accessing the wildcards using `{{wildcards.<name>}}.`"
            )

    def format_str(text):
        """Customized format function for strings."""
        assert isinstance(text, str), "Can only format `str` objects."
        is_input_function = output is None
        if is_input_function:
            text_fmt = smk.utils.format(text, **wildcards)
        else:
            check_for_clashes(text, wildcards, strict)
            if not strict:
                variables.update(wildcards)
            text_fmt = smk.utils.format(text, **variables)
        return text_fmt

    def format_dict(dictionary):
        """Customized format function for dictionaries."""
        for key, value in dictionary.items():
            if isinstance(value, str):
                dictionary[key] = format_str(value)
            elif isinstance(value, list):
                dictionary[key] = [format_str(v) for v in value]
            else:
                raise ValueError(
                    "Can only format `str` or list of `str` objects "
                    "within a dictionary."
                )
        return dictionary

    def formatter(obj):
        """Customized formatter function."""
        if isinstance(obj, str):
            obj_fmt = format_str(obj)
        elif isinstance(obj, dict):
            obj_fmt = format_dict(obj)
        else:
            raise ValueError("Can only format `str` or `dict` objects.")
        return obj_fmt

    return formatter


def switch_on_column(
    column, samples, options, match_on="tumour", format=True, strict=False
):
    """Pick an option based on the value of a column for a sample.

    The function finds the relevant row in `samples` for either
    the tumour (the default) or normal sample, which is determined
    by the `match_on` argument. To find the row, the `seq_type`
    and `tumour_id` (or `normal_id`) wildcards are required.

    The following special keys are available:

    _default
        If you provide a value under the key '_default' in `options`,
        this value will be used if the column value is not among
        the other keys in `options` (instead of defaulting to "").
    _prefix, _suffix
        If you provide values for the '_prefix' and/or '_suffix' keys
        in `options`, these values will be prepended and/or appended,
        respectively, to the selected value (including '_default') as
        long as the selected value is a string (not a dictionary).

    Parameters
    ----------
    column : str
        The column name whose value determines the option to pick.
    samples : pandas.DataFrame
        The samples data frame for the current module.
    options : dict
        The mapping between the possible values in `column` and the
        corresponding options to be returned. Special key-value
        pairs can also be included (see above).
    match_on : {"tumour", "normal"}
        Whether to match on the `sample_id` column in `samples`
        using `wildcard.tumour_id` or `wildcard.normal_id`.
    format : boolean
        Whether to format the option using the rule variables.
    strict : boolean
        Whether to include the bare wildcards in formatting.
        For example, if you have a wildcards called 'seq_type',
        without strict mode, you can access it with `{seq_type}`
        or `{wildcards.seq_type}`, whereas in strict mode, only
        the latter option is possible. This mode is useful if a
        wildcard has the same name as a rule variable, namely
        wildcards, input, output, threads, resources.

    Returns
    -------
    function
        A Snakemake-compatible input file or parameter function.
    """

    assert isinstance(options, dict), "`options` must be a `dict` object."
    assert column in samples, "`column` must be a column name in `samples`."

    def _switch_on_column(
        wildcards, input=None, output=None, threads=None, resources=None
    ):
        """Customized Snakemake input file or param function."""
        if match_on == "tumour":
            sample_id = wildcards.tumour_id
        elif match_on == "normal":
            sample_id = wildcards.normal_id
        else:
            raise ValueError("Invalid value for `match_on`.")
        subset = samples.loc[samples["seq_type"] == wildcards.seq_type]
        row = subset.loc[subset["sample_id"] == sample_id]
        assert len(row) == 1, (
            f"More than one row (or no row) matched the given seq_type "
            f"({wildcards.seq_type}) and given sample ID ({sample_id})."
        )
        series = row.squeeze()
        column_value = series[column]
        default_option = options.get("_default", "")
        selected_option = options.get(column_value, default_option)
        if isinstance(selected_option, str):
            prefix = options.get("_prefix", "")
            suffix = options.get("_suffix", "")
            combined = [prefix, selected_option, suffix]
            selected_option = "".join(x for x in combined if x != "")
        if format:
            formatter = create_formatter(
                wildcards, input, output, threads, resources, strict
            )
            selected_option = formatter(selected_option)
        return selected_option

    return _switch_on_column


def switch_on_wildcard(wildcard, options, format=True, strict=False):
    """Pick an option based on the value of a wildcard for a run.

    The following special keys are available:

    _default
        If you provide a value under the key '_default' in `options`,
        this value will be used if the column value is not among
        the other keys in `options` (instead of defaulting to "").
    _prefix, _suffix
        If you provide values for the '_prefix' and/or '_suffix' keys
        in `options`, these values will be prepended and/or appended,
        respectively, to the selected value (including '_default') as
        long as the selected value is a string (not a dictionary).

    Parameters
    ----------
    wildcard : str
        The wildcard name whose value determines the option to pick.
    options : dict
        The mapping between the possible values in `column` and the
        corresponding options to be returned. Special key-value
        pairs can also be included (see above).
    format : boolean
        Whether to format the option using the rule variables.
    strict : boolean
        Whether to include the bare wildcards in formatting.
        For example, if you have a wildcards called 'seq_type',
        without strict mode, you can access it with `{seq_type}`
        or `{wildcards.seq_type}`, whereas in strict mode, only
        the latter option is possible. This mode is useful if a
        wildcard has the same name as a rule variable, namely
        wildcards, input, output, threads, resources.

    Returns
    -------
    function
        A Snakemake-compatible input file or parameter function.
    """

    assert isinstance(options, dict), "`options` must be a `dict` object."

    def _switch_on_wildcard(
        wildcards, input=None, output=None, threads=None, resources=None
    ):
        """Customized Snakemake input file or param function."""
        wildcard_value = wildcards.get(wildcard)
        default_option = options.get("_default", "")
        selected_option = options.get(wildcard_value, default_option)
        if isinstance(selected_option, str):
            prefix = options.get("_prefix", "")
            suffix = options.get("_suffix", "")
            combined = [prefix, selected_option, suffix]
            selected_option = "".join(x for x in combined if x != "")
        if format:
            formatter = create_formatter(
                wildcards, input, output, threads, resources, strict
            )
            selected_option = formatter(selected_option)
        return selected_option

    return _switch_on_wildcard


def locate_bam(
    bam_directory=None,
    sample_keys=("sample_id", "tumour_id", "normal_id"),
    sample_bams=("sample_bam", "tumour_bam", "normal_bam"),
):
    """Locates BAM file for a given sample ID in a directory.

    This function actually configures another function, which is
    returned to be used by Snakemake.

    Parameters
    ----------
    bam_directory : str, optional
        The directory containing all BAM files. If None is provided,
        then the default value of 'data/' will be used.
    sample_keys : list of str, optional
        The possible wildcards that contain identifiers for samples
        with BAM files.
    sample_bams : list of str, optional
        The respective names for the BAM file located for each sample
        in `sample_keys` in the dictionary returned by the input file
        function. For example, the BAM file for the sample specified
        in 'sample_id' wildcard will be stored under the key
        'sample_bam' in the returned dictionary.

    Returns
    -------
    function
        A Snakemake-compatible input file function taking wildcards as
        its only argument. This function will return a dictionary of
        BAM files for any wildcards appearing in `sample_keys` under
        the corresponding keys specified in `sample_bams`.
    """

    bam_directory = "data/" if bam_directory is None else bam_directory

    assert len(sample_keys) == len(sample_bams)
    key_to_bam = dict(zip(sample_keys, sample_bams))

    # Use [None] to differentiate from the case where no BAM files are found
    bam_files = [None]
    seq_type_memo = dict()

    def locate_bam_custom(wildcards):

        # Retrieve list of BAM files (if not already done)
        if bam_files == [None]:
            del bam_files[0]
            bam_files.extend(list_files(bam_directory, "bam"))

        # Retrieve (or create) BAM files by seq_type
        seq_type = wildcards.seq_type
        if seq_type not in seq_type_memo:
            seq_type_memo[seq_type] = [b for b in bam_files if seq_type in b]
        seqtype_bam_files = seq_type_memo[seq_type]

        # Create dictionary meant for unpacking in Snakemake
        bams = dict()
        for sample_key, sample_id in wildcards.items():
            if sample_key not in sample_keys:
                continue
            matches = [b for b in seqtype_bam_files if sample_id in b]
            assert len(matches) == 1, (
                f"The given sample ID ({sample_id}) and seq_type "
                f"({wildcards.seq_type}) failed to identify a unique "
                f"BAM file in the given directory ({bam_directory}). "
                f"Instead, {len(matches)} matching files were found: "
                f"{', '.join(matches)}"
            )
            bam_name = key_to_bam[sample_key]
            bams[bam_name] = matches[0]

        return bams

    return locate_bam_custom


def check_reference(module_config, reference_key=None):
    """Ensure that a required reference config (and file) is available.

    If there is no 'genome_build' column in `module_samples` and
    there is only one loaded reference, this function will assume
    that the loaded reference is the reference to be used.

    Parameters
    ----------
    module_config : dict
        The module-specific configuration, corresponding to
        `config['lcr-modules']['<module-name>']`.
    reference_key : str, optional
        The key for a required reference file.

    Returns
    -------
    None
    """

    error_message = (
        "Load the appropriate reference YAML file in the `lcr-modules` repository. "
    )
    assert "reference" in module_config, error_message

    module_samples = module_config["samples"]
    references_available = list(module_config["reference"].keys())

    if len(references_available) == 0:
        raise AssertionError(error_message)

    if len(references_available) == 1 and "genome_build" not in module_samples:
        reference = references_available[0]
        log_message = f"Defaulting to the only loaded reference, {reference}"
        logger.warning(log_message)
        module_samples["genome_build"] = reference

    if len(references_available) > 1 and "genome_build" not in module_samples:
        raise AssertionError(
            "More than one loaded reference {references_available}, yet no "
            "`genome_build` column in the samples table. Thus, the genome build "
            "cannot be inferred. Add a `genome_build` column to your samples table."
        )

    genome_builds = module_samples["genome_build"].unique()

    for genome_build in genome_builds:
        assert genome_build in module_config["reference"], error_message
        ref_config = module_config["reference"][genome_build]
        if reference_key is not None:
            error_message += f"And ensure that `{reference_key}` is available."
            assert reference_key in ref_config, error_message


def get_reference(module_config, reference_key):
    check_reference(module_config, reference_key)

    def get_reference_custom(wildcards):
        return module_config["reference"][wildcards.genome_build][reference_key]

    return get_reference_custom


# SAMPLE PROCESSING


def load_samples(
    file_path, sep="\t", to_lowercase=("tissue_status",), renamer=None, **maps
):
    """Loads samples metadata with some light processing.

    The advantage of using this function over `pandas.read_table()`
    directly is that this function processes the data frame as follows:

        1) Can convert columns to lowercase.
        2) Can rename columns using either a renamer function or
           a set of key-value pairs where the values are the
           original names and the keys are the desired names.

    If a renamer function is provided in addition to a set of key-value
    pairs, the renamer function will be used first.

    Parameters
    ----------
    file_path : str
        The path to the tabular file containing the sample metadata
        (including any required columns).
    sep : str, optional
        The column separator.
    to_lowercase : list of str, optional
        The columns to be converted to lowercase.
    renamer : function or dict-like, optional
        A function that transforms each column name or a dict-like
        object that maps the original names (keys) to the desired
        names (values).
    **maps : key-value pairs, optional
        Pairs that specify the actual names (values) of the expected
        columns (keys). For example, if you had a 'sample' column
        while `lcr-modules` expects 'sample_id', you can use:

        load_samples(..., sample_id = "sample")

    Returns
    -------
    pandas.DataFrame
    """
    samples = pd.read_table(file_path, sep=sep)

    if renamer:
        samples.rename(columns=renamer, inplace=True)
    if maps:
        maps_rev = {v: k for k, v in maps.items()}
        samples.rename(columns=maps_rev, inplace=True)
    for col in to_lowercase:
        if col in samples.columns:
            samples[col] = samples[col].str.lower()
    return samples


def filter_samples(samples, invert=False, **filters):
    """Subsets for rows with certain values in the given columns.

    Parameters
    ----------
    samples : pandas.DataFrame
        The samples.
    invert : boolean
        Whether to keep or discard samples that match the filters.
    **filters : key-value pairs
        Columns (keys) and the values they need to contain (values).
        Values can be any value or a list of values.

    Returns
    -------
    pandas.DataFrame
        A subset of rows from the input data frame.
    """
    samples = samples.copy()
    for column, value in filters.items():
        if column not in samples:
            logger.warning(f"Column '{column}' not in sample table. Skipping.")
            continue
        if not isinstance(value, (list, tuple)):
            value = [value]
        if invert:
            samples = samples[~samples[column].isin(value)]
        else:
            samples = samples[samples[column].isin(value)]
    return samples


def keep_samples(samples, **filters):
    """Convenience wrapper around ``filter_samples``."""
    return filter_samples(samples, invert=False, **filters)


def discard_samples(samples, **filters):
    """Convenience wrapper around ``filter_samples``."""
    return filter_samples(samples, invert=True, **filters)


def group_samples(samples, subgroups):
    """Organizes samples into nested dictionary.

    Parameters
    ----------
    samples : pandas.DataFrame
        The samples.
    subgroups : list of str
        Columns of `samples` by which to organize the samples.
        The order determines the nesting order.

    Returns
    -------
    nested dict
        The number of levels is determined by the list of subgroups.
        The number of 'splits' at each level is based on the number of
        different values in the samples data frame for that column.
        The 'terminal' values are lists of samples, which are stored
        as named tuples containing all metadata for that row.
    """
    assert len(subgroups) > 0, "Need to provide at least one subgroup."
    # Iterate over each row
    samples_dict = dict()
    Sample = None
    for _index, row in samples.iterrows():
        values = []
        # Initialize intermediate subgroups with dictionaries
        parent = samples_dict
        for subgroup in subgroups[:-1]:
            value = row[subgroup]
            if value not in parent:
                parent[value] = {}
            values.append(value)
            parent = get_from_dict(samples_dict, values)
        # Initialize "terminal" subgroups with sets
        value = row[subgroups[-1]]
        values.append(value)
        if value not in parent:
            parent[value] = list()
        # Add sample ID to the "terminal" subgroup
        parent = parent[value]
        if Sample is None:
            Sample = namedtuple("Sample", row.index.tolist())
        sample = Sample(*row)
        assert (
            sample not in parent
        ), f"`{sample}` not unique for these subgroups ({values})."
        parent.append(sample)
    return samples_dict


def generate_runs_for_patient(
    patient_samples,
    run_paired_tumours,
    run_unpaired_tumours_with,
    unmatched_normal=None,
    unmatched_normals=None,
    run_paired_tumours_as_unpaired=False,
    **kwargs,
):
    """Generates a run for every tumour with and/or without a paired normal.

    Note that 'unpaired tumours' in the argument names and documentation
    refers to tumours without a matched normal sample.

    Parameters
    ----------
    patient_samples : dict
        Lists of sample IDs (str) organized by tissue_status (tumour vs
        normal) for a given patient. The order of the samples in each
        list is irrelevant.
    run_paired_tumours : boolean
        Whether to run paired tumours. Setting this to False is useful
        for naturally unpaired analyses (e.g., for RNA-seq).
    run_unpaired_tumours_with : { None, 'no_normal', 'unmatched_normal' }
        What to pair with unpaired tumours. This cannot be set to None if
        `run_paired_tumours_as_unpaired` is True. Provide value for
        `unmatched_normal` argument if this is set to 'unmatched_normal'.
    unmatched_normal : namedtuple, optional
        The normal sample to be used with unpaired tumours when
        `run_unpaired_tumours_with` is set to 'unmatched_normal'.
    unmatched_normals : dict, optional
        The normal samples to be used with unpaired tumours when
        `run_unpaired_tumours_with` is set to 'unmatched_normal'.
        Unlike `unmatched_normal`, this parameter expects a mapping
        from "{seq_type}--{genome_build}" to Sample namedtuples.
        If this option is provided, it will take precedence over
        `unmatched_normal`.
    run_paired_tumours_as_unpaired : boolean, optional
        Whether paired tumours should also be run as unpaired
        (i.e., separate from their matched normal sample).
        This is useful for benchmarking purposes or preventing
        unwanted paired analyses (e.g., in RNA-seq analyses
        intended to be tumour-only)
    **kwargs : key-value pairs
        Any additional unused arguments (e.g, `unmatched_normal_id`).

    Returns
    -------
    dict
        Lists of sample features prefixed with `tumour_` and `normal_`
        for all tumours for the given patient. Depending on the argument
        values, tumour-normal pairs may not be matching, and normal
        samples may not be included. The 'pair_status' column specifies
        whether a tumour is paired with a matched normal sample.
    """

    # Check that `run_unpaired_tumours_with` is among possible values
    run_unpaired_tumours_with_options = (None, "no_normal", "unmatched_normal")
    assert run_unpaired_tumours_with in run_unpaired_tumours_with_options, (
        "`run_unpaired_tumours_with` must be one of the values below "
        f"(not `{run_unpaired_tumours_with!r}`): \n"
        f"{run_unpaired_tumours_with_options}"
    )

    run_unpaired_tumour = run_unpaired_tumours_with is not None

    # Require `run_unpaired_tumour` if `run_paired_tumours_as_unpaired` is True
    assert run_unpaired_tumour or not run_paired_tumours_as_unpaired, (
        "`run_paired_tumours_as_unpaired` was True whereas "
        "`run_unpaired_tumours_with` was None. Please set "
        "`run_unpaired_tumours_with` to 'unmatched_normal' "
        "or 'no_normal'."
    )

    # Retrieve tumour and normal samples
    runs = defaultdict(list)
    tumour_samples = patient_samples.get("tumour", [])
    tumour_samples += patient_samples.get("tumor", [])
    normal_samples = patient_samples.get("normal", [None])

    # Add an unpaired normal is there isn't one
    if run_paired_tumours_as_unpaired and None not in normal_samples:
        normal_samples.append(None)

    for tumour, normal in itertools.product(tumour_samples, normal_samples):
        # Check for paired samples
        paired = normal is not None
        if paired and run_paired_tumours is False:
            continue
        # Check for unpaired samples
        unpaired = normal is None
        if unpaired and run_unpaired_tumour is False:
            continue
        # Compile features
        tumour = tumour._asdict()
        if normal is None and run_unpaired_tumours_with == "unmatched_normal":
            # Check that `unmatched_normal` or `unmatched_normals` is given
            seq_type = tumour["seq_type"]
            genome_build = tumour["genome_build"]
            assert unmatched_normal is not None or unmatched_normals is not None, (
                "`run_unpaired_tumours_with` was set to 'unmatched_normal' "
                f"whereas `unmatched_normal` and `unmatched_normals` were both "
                "None. For {seq_type!r}, provide an unmatched normal sample ID. "
                "See README for format."
            )
            if unmatched_normals is not None:
                normal = unmatched_normals[f"{seq_type}--{genome_build}"]._asdict()
            else:
                normal = unmatched_normal._asdict()
            runs["pair_status"].append("unmatched")
        elif normal is None and run_unpaired_tumours_with == "no_normal":
            normal = {key: None for key in tumour.keys()}
            runs["pair_status"].append("no_normal")
        else:
            normal = normal._asdict()
            runs["pair_status"].append("matched")
        for field in tumour.keys():
            runs["tumour_" + field].append(tumour[field])
            runs["normal_" + field].append(normal[field])

    return dict(runs)


def generate_runs_for_patient_wrapper(patient_samples, pairing_config):
    """Runs generate_runs_for_patient for the current seq_type/genome_build.

    This function is meant as a wrapper for `generate_runs_for_patient()`,
    whose parameters depend on the sequencing data type (seq_type) and
    genome_build of the samples at hand. It assumes that all samples for
    the given patient share the same seq_type and genome_build.

    Parameters
    ----------
    patient_samples : dict
        Same as `generate_runs_for_patient()`.
    pairing_config : nested dict
        The top level is sequencing data types (seq_type; keys) mapped
        to dictionaries (values) specifying argument values meant for
        `generate_runs_for_patient()`. For example:

            {'genome': {'run_unpaired_tumours_with': 'unmatched_normal',
                        'unmatched_normal': Sample(...)},
            'mrna': {'run_paired_tumour': False,
                    'run_unpaired_tumours_with': 'no_normal'}}

    Returns
    -------
    dict
        Same as `generate_runs_for_patient()`.
    """

    seq_type_set = set()
    for samples_list in patient_samples.values():
        seq_type_set.update(s.seq_type for s in samples_list)

    assert len(seq_type_set) == 1, (
        "This function is only meant to be run on groups of samples for a "
        "given patient and a given sequencing data type. The current group "
        f"of samples has the following seq_types: \n    {seq_type_set}"
    )

    seq_type = seq_type_set.pop()
    return generate_runs_for_patient(patient_samples, **pairing_config[seq_type])


def combine_lists(dictionary, as_dataframe=False):
    """Merges lists for matching keys in nested dictionary.

    Parameters
    ----------
    dictionary : dict
        Nested dictionaries where the key names match up.

        ::

            {'genome': {'field1': [1, 2, 3],
                        'field2': [4, 5, 6]},
            'mrna': {'field1': [11, 12, 13],
                    'field2': [14, 15, 16]}}

    as_dataframe : boolean, optional
        Whether the return value is coerced to pandas.DataFrame.

    Returns
    -------
    dict or pandas.DataFrame
        The type of the return value depends on `as_dataframe`.
        If `as_dataframe` is False, the output will look like:

        ::

            {'field1': [1, 2, 3, 11, 12, 13],
            'field2': [4, 5, 6, 14, 15, 16]}

        If `as_dataframe` is True, the output will look like:

        ::

                field1  field2
            0       1       4
            1       2       5
            2       3       6
            3      11      14
            4      12      15
            5      13      16
    """
    combined = defaultdict(list)
    for d in dictionary.values():
        for k, v in d.items():
            combined[k].extend(v)
    combined = dict(combined)
    if as_dataframe:
        combined = pd.DataFrame(combined)
    return combined


def walk_through_dict(
    dictionary, end_fn, max_depth=None, _trace=None, _result=None, **kwargs
):
    """Runs a function at a given level in a nested dictionary.

    If `max_depth` is unspecified, `end_fn()` will be run whenever
    the recursion encounters an object other than a dictionary.

    Parameters
    ----------
    dictionary : foo
        The dictionary to be recursively walked through.
    end_fn : function
        THe function to be run once recursion ends, either at
        `max_depth` or when a non-dictionary is encountered.
    max_depth : int, optional
        How far deep the recursion is allowed to go. By default, the
        recursion is allowed to go as deep as possible (i.e., until
        it encounters something other than a dictionary).
    _trace : tuple, optional
        List of dictionary keys used internally to track nested position.
    _result : dict
        Used internally to pass new dictionaries and avoid changing the
        input dictionary.
    **kwargs : key-value pairs
        Argument values that are passed to `end_fn()`.

    Returns
    -------
    dict
        A processed dictionary. The input dictionary remains unchanged.
    """

    # Define default values
    if max_depth is None:
        max_depth = float("inf")
    if _trace is None:
        _trace = tuple()
    if _result is None:
        _result = dict()

    # If the max_depth is zero, simply run `end_fn()` right away
    if max_depth <= 0:
        return end_fn(dictionary, **kwargs)

    # Iterate over every dictionary key and run `end_fn()` if the value
    # isn't a dictionary and the end depth isn't met. Otherwise, walk
    # through nested dictionary recursively.
    for k, v in dictionary.items():
        # Track nested position
        current_trace = _trace + (k,)
        if isinstance(v, dict) and len(current_trace) < max_depth:
            _result[k] = dict()
            walk_through_dict(v, end_fn, max_depth, current_trace, _result[k], **kwargs)
        else:
            _result[k] = end_fn(v, **kwargs)

    return _result


def generate_runs(
    samples,
    pairing_config=None,
    unmatched_normal_ids=None,
    subgroups=("seq_type", "genome_build", "patient_id", "tissue_status"),
):
    """Produces a data frame of tumour runs from a data frame of samples.

    Here, a 'tumour run' can consist of a tumour-only run or
    a paired run. In the case of a paired run, it can either
    be with a matched or unmatched normal sample.

    Parameters
    ----------
    samples : pandas.DataFrame
        The samples.
    pairing_config : dict, optional
        Same as `generate_runs_for_patient_wrapper()`. If left unset
        (or None is provided), this function will fallback on a
        default value (see `oncopipe.DEFAULT_PAIRING_CONFIG`).
    unmatched_normal_ids : dict, optional
        The mapping from seq_type and genome_build to the unmatched
        normal sample IDs that should be used for unmatched analyses.
        The keys must take the form of '{seq_type}--{genome_build}'.
    subgroups : list of str, optional
        Same as `group_samples()`.

    Returns
    -------
    pandas.DataFrame
        The generated runs with columns matching the keys of the
        return value for `generate_runs_for_patient()`.
    """

    # Set default values and create copy
    if pairing_config is None:
        pairing_config = DEFAULT_PAIRING_CONFIG
    pairing_config = copy.deepcopy(pairing_config)

    # Generate Sample instances for unmatched normal samples from sample IDs
    Sample = namedtuple("Sample", samples.columns.tolist())
    sample_genome_builds = samples["genome_build"].unique()
    for seq_type, args_dict in pairing_config.items():
        if (
            "run_unpaired_tumours_with" in args_dict
            and args_dict["run_unpaired_tumours_with"] == "unmatched_normal"
            and unmatched_normal_ids is not None
        ):
            unmatched_normals = dict()
            for key, normal_id in unmatched_normal_ids.items():
                _, genome_build = key.split("--", 1)
                if (
                    not key.startswith(f"{seq_type}--")
                    or genome_build not in sample_genome_builds
                ):
                    continue
                normal_row = samples[
                    (samples.sample_id == normal_id) & (samples.seq_type == seq_type)
                ]
                num_matches = len(normal_row)
                assert num_matches == 1, (
                    f"There are {num_matches} {seq_type} samples matching "
                    f"the normal ID {normal_id} (instead of just one)."
                )
                unmatched_normals[key] = Sample(*normal_row.squeeze())
            args_dict["unmatched_normals"] = unmatched_normals
        elif (
            "run_unpaired_tumours_with" in args_dict
            and args_dict["run_unpaired_tumours_with"] == "unmatched_normal"
            and "unmatched_normal_id" in args_dict
        ):
            normal_id = args_dict["unmatched_normal_id"]
            normal_row = samples[
                (samples.sample_id == normal_id) & (samples.seq_type == seq_type)
            ]
            num_matches = len(normal_row)
            assert num_matches == 1, (
                f"There are {num_matches} {seq_type} samples matching "
                f"the normal ID {normal_id} (instead of just one)."
            )
            args_dict["unmatched_normal"] = Sample(*normal_row.squeeze())

    # Organize samples by patient and tissue status (tumour vs. normal)
    patients = group_samples(samples, subgroups)

    # Find every possible tumour-normal pair for each patient
    end_depth = len(subgroups) - 1
    runs = walk_through_dict(
        patients,
        generate_runs_for_patient_wrapper,
        end_depth,
        pairing_config=pairing_config,
    )
    while end_depth > 0:
        runs = walk_through_dict(runs, combine_lists, end_depth - 1, as_dataframe=True)
        end_depth -= 1

    # Warn if runs have duplicates
    if any(runs.duplicated()):
        logger.warning("Duplicate runs exist. This probably shouldn't happen.")

    # Fix column names if data frame is empty
    if runs.empty:
        column_names = (
            ["pair_status"]
            + ("tumour_" + samples.columns).to_list()
            + ("normal_" + samples.columns).to_list()
        )
        runs = pd.DataFrame(columns=column_names)

    return runs


def generate_pairs(samples, unmatched_normal_ids=None, **seq_types):
    """Generate tumour-normal pairs using sensible defaults.

    Each sequencing data type (``seq_type``) is provided as
    separate arguments with a specified "pairing mode". This
    mode determines how the samples for that ``seq_type``
    are paired. Only the listed ``seq_type`` values will be
    included in the output. The available pairing modes are:

    1. ``matched_only``: Only tumour samples with matched
       normal samples will be returned. In other words,
       unpaired tumour or normal samples will be omitted.

       .. code:: python

          generate_pairs(SAMPLES, genome='matched_only')

    2. ``allow_unmatched``: All tumour samples will be returned
       whether they are paired with a matched normal sample
       or not. If they are not paired, they will be returned
       with an unmatched normal sample specified by the user.
       This mode must be specified alongside the ID for the
       sample to be paired with unpaired tumours as a tuple.
       This sample must be present in the ``samples`` table.

       .. code:: python

          generate_pairs(SAMPLES, genome=('allow_unmatched', 'PT003-N'))

    3. ``no_normal``: All tumour samples will be returned
       without a paired normal sample. This is simply a
       shortcut for filtering for tumour samples, but this
       ensures that the column names will be consistent
       with other calls to ``generate_pairs()``.

       .. code:: python

          generate_pairs(SAMPLES, mrna='no_normal')

    Parameters
    ----------
    samples : pandas.DataFrame
        The sample table. This data frame must include the
        following columns: ``sample_id``, ``patient_id``,
        ``seq_type``, and ``tissue_status`` ('normal' or
        'tumour'/'tumor'). If ``genome_build`` is included,
        no tumour-normal pairs will be made between different
        genome builds.
    unmatched_normal_ids : dict, optional
        The mapping from seq_type and genome_build to the unmatched
        normal sample IDs that should be used for unmatched analyses.
        The keys must take the form of '{seq_type}--{genome_build}'.
    **seq_types : {'matched_only', 'allow_unmatched', 'no_normal'}
        A mapping between values of ``seq_type`` and
        pairing modes. See above for description of each
        pairing mode.

    Returns
    -------
    pandas.DataFrame
        The tumour-normal pairs (one pair per row), but
        the normal sample is omitted if the ``no_normal``
        pairing mode is used. Every column in the input
        ``samples`` data frame will appear twice in the
        output, once for the tumour sample and once for
        the normal sample, prefixed by ``tumour_`` and
        ``normal_``, respectively. An additional column
        called ``pair_status`` will indicate whether the
        tumour-normal samples in the row are matched or
        unmatched. If the normal sample is omitted due
        to the ``no_normal`` mode, this column will be
        set to ``no_normal``.

    Examples
    --------
    Among the samples in the ``SAMPLES`` data frame, the
    ``genome`` tumour samples will be paired with a matched
    normal samples if one exists or with the given unmatched
    normal sample (``PT003-N``) if no matched normal samples
    are present; the ``capture`` tumour samples will only be
    paired with matched normal samples; and the ``mrna``
    tumour samples will be returned without matched or
    unmatched normal samples.

    >>> PAIRS = generate_pairs(SAMPLES, genome=('allow_unmatched', 'PT003-N'),
    >>>                        capture='matched_only', mrna='no_normal')
    """

    # Define pairing modes
    pairing_modes = {
        "matched_only": {
            "run_paired_tumours": True,
            "run_unpaired_tumours_with": None,
            "run_paired_tumours_as_unpaired": False,
        },
        "allow_unmatched": {
            "run_paired_tumours": True,
            "run_unpaired_tumours_with": "unmatched_normal",
            "run_paired_tumours_as_unpaired": False,
            # unmatched_normal_id must be added
        },
        "no_normal": {
            "run_paired_tumours": False,
            "run_unpaired_tumours_with": "no_normal",
            "run_paired_tumours_as_unpaired": True,
        },
    }

    # Iterate over seq_types
    pairing_config = dict()
    available_pairing_modes = list(pairing_modes.keys())
    for seq_type, mode in seq_types.items():

        # Check if mode was provided as a two-element iterable (list or tuple)
        unmatched_normal_id = None
        if len(mode) == 2 and mode[0] == "allow_unmatched":
            unmatched_normal_id = mode[1]
            mode = "allow_unmatched"

        # Make sure mode is a string and among the available options
        assert isinstance(mode, str) and mode in available_pairing_modes, (
            f"The pairing mode specified for {seq_type!r} isn't valid. "
            f"The available modes are: {available_pairing_modes}."
        )

        # Make sure that the `allow_unmatched` mode is provided
        # with unmatched_normal_id or unmatched_normal_ids
        assert mode != "allow_unmatched" or (
            unmatched_normal_id is not None or unmatched_normal_ids is not None
        ), (
            "The 'allow_unmatched' mode must be provided with the "
            "`unmatched_normal_ids` parameter or paired with a normal "
            "sample ID, such as:\n    "
            "generate_pairs(SAMPLES, genome=('allow_unmatched', 'PT003-N'))\n"
        )

        pairing_config[seq_type] = pairing_modes[mode]

        if mode == "allow_unmatched" and unmatched_normal_id is not None:
            pairing_config[seq_type]["unmatched_normal_id"] = unmatched_normal_id

    # Subgroup using `genome_build` if available
    if "genome_build" in samples:
        subgroups = ("seq_type", "genome_build", "patient_id", "tissue_status")
    else:
        subgroups = ("seq_type", "patient_id", "tissue_status")

    # Make sure all of the required columns are present
    required_columns = list(subgroups) + ["sample_id"]
    assert all(column in samples for column in required_columns), (
        "The sample table doesn't include all of the "
        f"expected columns, namely {required_columns}."
    )

    # Generate the runs using the generated pairing configuration
    samples = filter_samples(samples, seq_type=list(seq_types.keys()))
    runs = generate_runs(samples, pairing_config, unmatched_normal_ids, subgroups)

    return runs


# MODULE SETUP/CLEANUP


def check_for_none_strings(config, name):
    """Warn the user if 'None'/'null' strings are found in config."""

    def check_for_none_strings_(obj):
        if isinstance(obj, str):
            if obj in ["None", "null"]:
                logger.warning(
                    f"Found the value `{obj!r}` (string) in the configuration for "
                    f"the {name} module. You probably intended to use the value "
                    "`null` (without quotes) in the YAML configuration files."
                )
        return obj

    walk_through_dict(config, check_for_none_strings_)


def check_for_update_strings(config, name):
    """Warn the user if '__UPDATE__' strings are found in config."""

    def check_for_update_strings_(obj):
        if isinstance(obj, str):
            assert "__UPDATE__" not in obj, (
                "Found the value '__UPDATE__' in the configuration for the "
                f"{name} module. This usually means some values from the "
                "module's default configuration file haven't been updated. "
                f"For more info, check out {DOCS['update_config']}."
            )
        return obj

    walk_through_dict(config, check_for_update_strings_)


def setup_module(name, version, subdirectories):
    """Prepares and validates configuration for the given module.

    This function performs a number of convenient tasks:

        1) It ensures that the `CFG` variable doesn't exist. This is
           intended as a safeguard since the modules use `CFG` as a
           convenient shorthand.
        2) It ensures that Snakemake meets the required version.
        3) It ensures that the required configuration is loaded.
        4) It initializes the module configuration with the `_shared`
           configuration, but recursively overwrites values from the
           module-specific configuration. In other words, the
           specific overrides the general.
        5) It ensures that the module configuration has the expected
           fields to avoid errors downstream.
        6) It's updates any strings containing placeholders such as
           `{REPODIR}`, `{MODSDIR}`, and `{SCRIPTSDIR}` with the
           actual values.
        7) It validates the samples table using all of the schema
           YAML files in the module's `schemas/` folder.
        8) It configures, numbers, and creates the output and
           log subdirectories.
        9) It generates a table of runs consisting of tumour-
           normal pairs in case that's useful.
        10) It will automatically filter the samples for those
            whose `seq_type` appear in `pairing_config`.

    Parameters
    ----------
    name : str
        The name of the module.
    version : str
        The semantic version of the module.
    subdirectories : list of str
        The subdirectories of the module output directory where the
        results will be produced. They will be numbered incrementally
        and created on disk. This should include 'inputs' and 'outputs'.

    Returns
    -------
    dict
        The module-specific configuration, including any shared
        configuration from `config['lcr-modules']['_shared']`.
    """

    # Get namespace where module is being set up
    module_frame = inspect.currentframe().f_back
    module_globals = module_frame.f_globals
    config = module_globals["config"]

    # Make sure the `CFG` variable doesn't exist yet
    assert "CFG" not in module_globals, "`CFG` is a reserved variable for lcr-modules."

    # Ensure minimum version of Snakemake
    smk.utils.min_version("5.4.0")

    # Ensure that the lcr-modules _shared config is loaded
    assert "lcr-modules" in config and "_shared" in config["lcr-modules"], (
        "Shared lcr-modules configuration is not loaded. "
        "See README.md in lcr-modules for more information."
    )

    # Ensure that this module's config is loaded
    assert name in config["lcr-modules"], (
        f"The configuration for the {name!r} module is not loaded. "
        "It should be loaded before the module Snakefile (.smk) is "
        "included. See README.md for more information."
    )

    # Get configuration for the given module and create samples shorthand
    mconfig = copy.deepcopy(config["lcr-modules"]["_shared"])
    smk.utils.update_config(mconfig, config["lcr-modules"][name])
    msamples = mconfig["samples"].copy()

    # Check whether there are "None" strings
    check_for_none_strings(mconfig, name)

    # Check whether there are "__UPDATE__" strings
    check_for_update_strings(mconfig, name)

    # Drop samples whose seq_types do not appear in pairing_config
    assert "pairing_config" in mconfig, "`pairing_config` missing from module config."
    sample_seq_types = msamples["seq_type"].unique()
    pairing_config = mconfig["pairing_config"]
    supported_seq_types = [
        k for k, v in pairing_config.items() if "run_paired_tumours" in v
    ]
    unsupported_seq_types = set(sample_seq_types) - set(supported_seq_types)
    if len(unsupported_seq_types) > 0:
        logger.warning(
            f"Some samples have seq_types {unsupported_seq_types} that are "
            f"not configured in the pairing config for the {name} module. "
            "They will be excluded from the analysis."
        )
    msamples = msamples[msamples["seq_type"].isin(supported_seq_types)].copy()
    mconfig["samples"] = msamples

    # Set module name and version
    mconfig["name"] = name
    mconfig["version"] = version

    # Ensure that common module sub-fields are present
    subfields = ["inputs", "dirs", "conda_envs", "options", "threads", "mem_mb"]
    for subfield in subfields:
        if subfield not in mconfig:
            mconfig[subfield] = dict()

    # Check reference
    assert (
        "genome_build" in msamples
    ), "Add a `genome_build` column to your samples data frame."

    # Update placeholders in any string in the module-specific config
    def update_placeholders(obj, **placeholders):
        if isinstance(obj, str):
            result = obj
            for placeholder, value in placeholders.items():
                result = result.replace("{" + placeholder + "}", value)
        else:
            result = obj
        return result

    # Find repository and module directories
    repodir = os.path.normpath(mconfig["lcr-modules"])
    modsdir = os.path.join(repodir, "modules", name, version)
    scriptsdir = os.path.normpath(mconfig["lcr-scripts"])

    placeholders = {
        "REPODIR": repodir,
        "MODSDIR": modsdir,
        "SCRIPTSDIR": scriptsdir,
    }
    mconfig = walk_through_dict(mconfig, update_placeholders, **placeholders)

    # Validate samples data frame
    schemas_dir = os.path.join(modsdir, "schemas")
    schemas = os.listdir(schemas_dir)
    for schema in schemas:
        smk.utils.validate(msamples, schema=os.path.join(schemas_dir, schema))

    # Configure output directory if not specified and create it
    if mconfig["dirs"].get("_parent") is None:
        root_output_dir = mconfig.get("root_output_dir", "results")
        output_dir = os.path.join(root_output_dir, f"{name}-{version}")
        mconfig["dirs"]["_parent"] = output_dir
    mconfig["dirs"]["_parent"] = mconfig["dirs"]["_parent"].rstrip("/") + "/"
    os.makedirs(mconfig["dirs"]["_parent"], exist_ok=True)

    # Update paths to conda environments to be relative to the module directory
    for env_name, env_val in mconfig["conda_envs"].items():
        if env_val is not None:
            mconfig["conda_envs"][env_name] = os.path.relpath(env_val, modsdir)

    # Setup output sub-directories
    scratch_subdirs = mconfig.get("scratch_subdirectories", [])
    mconfig = setup_subdirs(mconfig, subdirectories, scratch_subdirs)

    # Setup log sub-directories
    mconfig["logs"] = dict()
    parent_dir = mconfig["dirs"]["_parent"]
    launched_fmt = _session.launched_fmt
    logs_parent_dir = os.path.join(parent_dir, "logs", launched_fmt)
    logs_parent_dir = logs_parent_dir.rstrip("/") + "/"
    mconfig["logs"]["_parent"] = logs_parent_dir
    os.makedirs(logs_parent_dir, exist_ok=True)
    for subdir, value in mconfig["dirs"].items():
        if subdir == "_parent":
            continue
        mconfig["logs"][subdir] = value.replace(parent_dir, logs_parent_dir)

    # Generate runs
    assert "pairing_config" in mconfig, "Module config must have 'pairing_config'."
    runs = generate_runs(
        msamples, mconfig["pairing_config"], mconfig.get("unmatched_normal_ids")
    )

    # Split runs based on pair_status
    mconfig["runs"] = runs
    mconfig["paired_runs"] = runs[runs.pair_status != "no_normal"]
    mconfig["unpaired_runs"] = runs[runs.pair_status == "no_normal"]

    # Return module-specific configuration
    return mconfig


def setup_subdirs(module_config, subdirectories, scratch_subdirs=()):
    """Numbers and creates module output subdirectories.

    Parameters
    ----------
    module_config : dict
        The module-specific configuration.
    subdirectories : list of str
        The names (without numbering) of the output subdirectories.
    scratch_subdirs : list of str, optional
        A subset of `subdirectories` that should be symlinked into the given
        scratch directory, specified under:

            `config["lcr_modules"]["_shared"]["scratch_directory"]`

        This should not include 'inputs' and 'outputs', which only
        contain symlinks.

    Returns
    -------
    dict
        The updated module-specific configuration with the paths
        to the numbered output subdirectories.
    """

    # Check for any issues with the subdirectory names
    assert "_parent" not in subdirectories, "You cannot have a '_parent' sub-directory."
    assert subdirectories[0] == "inputs", "The first subdirectory must be 'inputs'."
    assert subdirectories[-1] == "outputs", "The last subdirectory must be 'outputs'."

    # Check if `scratch_directory` is needed
    if len(scratch_subdirs) > 0 and "scratch_directory" not in module_config:
        raise AssertionError(
            "`scratch_directory` is not specified in the `_shared` config. If "
            "you don't want to use a scratch directory, set the following "
            "(note the lack of quotes around null): \nscratch_directory: null"
        )

    # If `scratch_directory` is None, then don't worry about `scratch_subdirs`
    scratch_directory = module_config.get("scratch_directory")
    if scratch_directory is None:
        scratch_subdirs = ()

    # Check it 'inputs' or 'outputs' are among the `scratch_subdirs`
    msg = "'inputs' and 'outputs' cannot be `scratch_subdirs`."
    assert "inputs" not in scratch_subdirs, msg
    assert "outputs" not in scratch_subdirs, msg

    # Generate zero-padded numbers
    numbers = [f"{x:02}" for x in (*range(0, len(subdirectories) - 1), 99)]

    # Join numbers and names, and create subdirectory
    name = module_config["name"]
    version = module_config["version"]
    parent_dir = module_config["dirs"]["_parent"]
    for num, subdir in zip(numbers, subdirectories):
        subdir_full = os.path.join(parent_dir, f"{num}-{subdir}/")
        module_config["dirs"][subdir] = subdir_full
        if subdir in scratch_subdirs:
            scratch_parent_dir = os.path.join(scratch_directory, f"{name}-{version}")
            scratch_subdir_full = os.path.join(scratch_parent_dir, f"{num}-{subdir}/")
            os.makedirs(scratch_subdir_full, exist_ok=True)
            relative_symlink(scratch_subdir_full, subdir_full, overwrite=False)
        else:
            os.makedirs(subdir_full, exist_ok=True)

    return module_config


def cleanup_module(module_config):
    """Save module-specific configuration, sample, and runs to disk."""

    # Get namespace where module is being cleaned up
    module_frame = inspect.currentframe().f_back
    module_globals = module_frame.f_globals

    # Delete `CFG` from module namespace to avoid conflicts
    del module_globals["CFG"]

    # Define useful variables
    parent_dir = module_config["logs"]["_parent"]

    # Define fields to be output as TSV files
    tsv_fields = {
        "samples": None,
        "runs": None,
        "paired_runs": None,
        "unpaired_runs": None,
    }
    tsv_fields_skip = ["paired_runs", "unpaired_runs"]
    for field in tsv_fields.keys():
        if field not in module_config:
            continue
        tsv_fields[field] = module_config.pop(field)
        if field not in tsv_fields_skip:
            output_file = os.path.join(parent_dir, f"{field}.tsv")
            tsv_fields[field].to_csv(output_file, sep="\t", index=False)

    # Output current configuration for future reference
    config_file = os.path.join(parent_dir, "config.yaml")
    with open(config_file, "w") as config_file_handler:
        yaml.dump(module_config, config_file_handler)

    # Add back the TSV fields
    for field in tsv_fields.keys():
        module_config[field] = tsv_fields[field]
