#!/usr/bin/env Rscript

library(tidyverse)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# run script in base Ce-PopGen directory

# args
# 1 - analysis type - name of folder corresponding to what admixure analysis

# example - 
# args <- c("GATK-STRELKA_Intersection_Complete")
args <- commandArgs(TRUE)

analysis_type <- args[1]

pr_mask_directory <- "Data/Divergent_Masks/Masks_By_Type/"

test <- enrichment_results[[4]]@result

head(test)
