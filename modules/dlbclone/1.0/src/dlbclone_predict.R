#!/usr/bin/env Rscript
#
# Notes: 
#   This script is intended for use with the dlbclone-1.0 module in LCR-modules. 
#   It uses GAMBL.predict stored DLBCLone models OR custom DLBCLone models created
#   by users using rule: _dlbclone_build_model. 
#
#   This script predicts (GAMBLR.predict::DLBCLone_predict()) DLBCL molecular subgroups 
#   for one or more samples using a pre-trained DLBCLone model. The pretrained DLBCLone
#   model is loaded (GAMBLR.predict::DLBCLone_load_optimized()) and activated
#   (GAMBLR.predict::DLBCLone_activate()) prior to prediction. A predictions tsv file 
#   is then generated for user analysis of their testing data. 
#
#   The user must appropriately set the config file with a valid DLBCLone model they 
#   want to predict their test samples on. 

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressPackageStartupMessages({
  library(GAMBLR.predict)
  library(optparse)
  library(tibble)
  library(dplyr)
  library(readr)
})

# Determine arguments -----------------------------------------------------

# Define command line options
option_list <- list(
    make_option(c("-f", "--test_features"), type = "character", help = "Input test mutation matrix file path"),
    make_option(c("-o", "--output_dir"), type = "character", help = "Output file path for predictions"),
    make_option(c("-p", "--model_path"), type = "character", help = "Pre-trained model file path"),
    make_option(c("-n", "--model_prefix"), type = "character", help = "Pre-trained model name prefix"),
    make_option(c("-m", "--fill_missing"), action="store_true", help="if TRUE, any features present in the model but not seen in test_features are added and set to zero. If FALSE, missing model features cause an error"),
    make_option(c("-e", "--drop_extra"), action="store_true", help="if TRUE, any features present in test_features but not seen during model training are dropped")
)

opt <- parse_args(OptionParser(option_list=option_list))

test_features_path <- opt$test_features
test_features <- read_tsv(test_features_path) %>%
    column_to_rownames("sample_id")

model_path <- opt$model_path
if (model_path == "GAMBLR.predict"){
    model_path = system.file("extdata/models",package = "GAMBLR.predict")
}

output_dir <- opt$output_dir
model_prefix <- opt$model_prefix
fill_missing <- opt$fill_missing
drop_extra <- opt$drop_extra

# Processing stored DLBCLone model object ---------------------------------

# Run GAMBLR.predict::DLBCLone_load_optimized()
loaded_model <- DLBCLone_load_optimized( 
    path = model_path,
    name_prefix = model_prefix
)

# Run GAMBLR.predict::DLBCLone_activate()
active_model <- DLBCLone_activate( 
    loaded_model, 
    force = TRUE
)

# Run GAMBLR.predict::DLBCLone_predict() ----------------------------
predictions <- DLBCLone_predict(
    mutation_status = test_features, 
    optimized_model = active_model,
    fill_missing = fill_missing,
    drop_extra = drop_extra
)

# Output ------------------------------------------------------------

message("Writing DLBCLone predictions to file...")
write_tsv(predictions$prediction, output_dir)
