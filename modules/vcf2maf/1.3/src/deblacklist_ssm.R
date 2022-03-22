#! /usr/bin/env Rscript

# Uses GAMBLR to remove blacklisted genes from augmented MAFs

log <- file(snakemake@log[[1]], open="wt")
sink(stdout(), type = "message")
sink(log)

library(GAMBLR)
library(data.table)
library(tidyverse)
library(readr)
library(glue)

genome_build = ifelse(
    str_detect(snakemake@wildcards[["genome_build"]], "38"), 
    "hg38", 
    "grch37"
)

blacklist_template = snakemake@params[["blacklist_template"]]
project_base = snakemake@params[["project_base"]]
seq_type_param = snakemake@params[["seq_type_blacklist"]]
drop_threshold = snakemake@params[["threshold"]]

maf = snakemake@input[["maf"]]

seq_type = if (seq_type_param == "both") c("genome", "capture") else snakemake@wildcards[["seq_type"]]

fread_maf(maf) %>% 
    mutate(Chromosome = as.character(Chromosome)) %>% 
    annotate_ssm_blacklist(
        genome_build = genome_build, 
        blacklist_file_template = blacklist_template, 
        project_base = project_base, 
        seq_type = seq_type, 
        verbose = TRUE, 
        drop_threshold = drop_threshold
    ) %>% 
    write_tsv(snakemake@output[["maf"]])

sink()
