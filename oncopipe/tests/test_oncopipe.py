import oncopipe as op
import pandas as pd
import snakemake as smk
from snakemake.io import Namedlist
import yaml
import copy as cp

##### INPUTS #####
tsv = "sample_test.tsv"
csv = "sample_test.csv"

##### VARIABLES #####
#wildcards = Namedlist(fromdict = {"seq_type": "capture"})


##### LOAD GENERAL AND MODULE (MANTA) CONFIG FILES #####
DF = op.load_samples(tsv, sep = "\t") # pd.read_table(tsv, sep = "\t")

with open(r"config.yaml") as file:
    config = yaml.full_load(file)

with open(r"manta.yaml") as file:
    mconfig = yaml.full_load(file)

smk.utils.update_config(config, mconfig)
op.enable_set_functions(config)
op.set_samples("_shared", DF)
op.set_samples("manta", DF)

##### UNIT TEST FOR FUNCTIONS IN ONCOPIPE #####

def test_load_samples():
    # func -> load_samples(file_path, sep="\t", lowercase_cols=LOWERCASE_COLS, renamer=None, **maps)
    r_tsv = cp.deepcopy(DF) #op.load_samples(tsv, sep = "\t")
    r_csv = op.load_samples(csv, sep = ",")
    r_lowercase = op.load_samples(tsv, sep = "\t", lowercase_cols = "ff_or_ffpe")
    assert type(r_tsv) and type(r_csv) is pd.DataFrame
    assert "FF" or "FFPE" not in r_lowercase.columns


def test_filter_samples(): 
    # func -> filter_samples(samples, **filters)
    result = op.filter_samples(DF, seq_type = ["mrna", "capture"])
    assert type(result) is pd.DataFrame
    assert "mrna" and "capture" in list(result["seq_type"])
    assert "genome" not in result["seq_type"]


def test_group_samples():
    # func -> group_samples(samples, subgroups)
    result = op.group_samples(DF, ["seq_type", "patient_id", "tissue_status"])
    key = list(DF.seq_type.unique())

    assert type(result) is pd.DataFrame
    assert type(result.get("capture").get("patient1").get("Tumour")) is list
    for i in key:
        assert i in result.keys()


def test_set_samples():
    # func -> set_samples(module, *samples)
    op.set_samples("manta", DF)
    assert type(config) is dict
    assert type(config["lcr-modules"]["manta"]["samples"]) is pd.DataFrame


def test_set_input():
    """
    func -> set_input(module, name, value)
    also tests get_from_dict(dictionary, list_of_keys) function
    """
    op.set_input("manta", "bam", "bwa_mem")
    assert type(config) is dict
    assert type(config["lcr-modules"]["manta"]["inputs"]) is dict
    assert "bam" in config["lcr-modules"]["manta"]["inputs"]
    assert "bam" in op.get_from_dict(config, ["lcr-modules", "manta", "inputs"])
    assert config["lcr-modules"]["manta"]["inputs"]["bam"] == "bwa_mem"
    assert op.get_from_dict(config, ["lcr-modules", "manta", "inputs", "bam"]) == "bwa_mem"

"""
## INCOMPLETE
def test_setup_module():
    # func -> setup_module(name, version, subdirectories)
    tconfig = copy.deepcopy(config)
    result = op.setup_module("manta", "1.0", ["inputs", "manta", "outputs"])
"""

def test_as_one_line():
    # func -> as_one_line(text)
    result = op.as_one_line("""
        line 1     
            line 2
        line 3
        """)
    assert result == "line 1 line 2 line 3"
    assert type(result) is str

"""
def test_get_reference():
    # get_reference(module_config, reference_key)
    wildcards = Namedlist(fromdict = {"genome_build": "grch37"})
    #result = op.get_reference(config["lcr-modules"]["manta"], "genome_fasta")
    assert op.get_reference(config["lcr-modules"]["manta"], "genome_fasta") == "grch37.fa"
    assert op.get_reference(config["lcr-modules"]["manta"], "intervals") == "grch37_intervals.txt"



##### TESTING COMPLICATED FUNCTIONS #####

def test_switch_on():
    #func -> switch_on_column(
    #column, samples, options, match_on="tumour", format=True, strict=False)
    # setup variables and results
    wildcards = Namedlist(fromdict = {"seq_type": "capture", "tumour_id": "A88333"})
    op.set_samples("manta", DF)

    r_col = op.switch_on_column("genome_build", config["lcr-modules"]["manta"]["samples"], {"_default": "no_genome", "grch37": "GRCh37_genome"})

    r_wc = op.switch_on_column("seq_type", config["lcr-modules"]["manta"]["samples"], config["lcr-modules"]["manta"]["options"]["configure"])

    # assertions
    #assert r_col == "GRCh37_genome"
    #assert r_wc == config["lcr-modules"]["manta"]["options"]["configure"]["capture"]

"""