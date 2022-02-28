# Changelog

All notable changes to the `lymphgen` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-11-02

This release was authored by Chris "coolbeans" Rushton.

Initial release, adding LymphGen, LGenIC, and allowing for CNV data (if availible)
You can run LymphGen using just SNVs, or with CNV and SV data
Note if you provide CNV and SV files, you should specify the appropriate column names in the config file
All possible iterations of LymphGen will be run (i.e. if you provide both CNVs and SNVs, LymphGen will be run
with both CNVs and SNVs, as well as just with SNVs)
