import oncopipe as op
import pandas as pd

##### INPUTS #####
data = "sample_test.tsv"
DF = pd.read_table(dataframe, sep = "\t")

def test_load_samples():
    result = op.load_samples(data)

def test_set_input():
    # func -> set_input(module, name, value)
    result = op.set_input("manta", "bam_files", "bam/")
    assert result == {"lcr_modules": {manta: {"inputs": {"bam_files": "bam/"}}}}
    assert type(result) is dict

def test_set_samples():
    # func -> set_samples(module, *samples)
    result = op.set_samples("manta", )