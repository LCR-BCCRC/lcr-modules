# Changelog

## [1.0] - 2020-10-30

This release was authored by Ryan Morin.
I am using a modified Starfish Python script because the one in github was buggy. This tool is packaged with the modified version. Ideally we will shift to the up to date version eventually to avoid this. This version is rough but meets the needs of the developer. 
The conversion of VCFs to bed and merging of indels is tailored to the goal of running the Mutect module on candidate indels from Strelka while only running on the SNVs commont to LoFreq and Strelka.

## [2.0] - 2020-11-17

This release was authored by Laura Hilton. 

 - This version is compatible with up to 6 input VCF files, and generates a union VCF file from the output. 
 - The included src/starfish.py script is modified from the developer verson to throw an error if the user expects a venn diagrm output but the requirements aren't satisfied. 
 **WARNING:** Despite the use of a target sentinel file, Starfish still occasionally outputs empty VCF files without throwing any errors. Check your results carefully. 