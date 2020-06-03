#!/usr/bin/env python3

from setuptools import setup

setup(
    name="oncopipe",
    version="1.0",
    description="Functions for running Snakemake modules",
    url="https://github.com/LCR-BCCRC/lcr-modules",
    author="Bruno Grande",
    author_email="bgrande@sfu.ca",
    license="MIT",
    py_modules=["oncopipe"],
    install_requires=["pyyaml", "pandas", "snakemake>=5.4",],
    zip_safe=False,
    python_requires=">3.6.0",
)
