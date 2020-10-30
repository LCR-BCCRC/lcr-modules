# Changelog

## [1.0] - 2020-10-30

This release was authored by Ryan Morin.
I am using a modified Starfish Python script because the one in github was buggy. This tool is packaged with the modified version. Ideally we will shift to the up to date version eventually to avoid this. This version is rough but meets the needs of the developer. 
The conversion of VCFs to bed and merging of indels is tailored to the goal of running the Mutect module on candidate indels from Strelka while only running on the SNVs commont to LoFreq and Strelka.