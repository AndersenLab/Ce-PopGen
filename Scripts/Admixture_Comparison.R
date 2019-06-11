#!/usr/bin/env Rscript

library(tidyverse)
library(CCA)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

analysis_type_a <- "GATK-STRELKA_Intersection"
analysis_type_b <- "DIVERGENT-MASKED"

analysis_a <- data.table::fread(glue::glue("Processed_Data/Admixture/{analysis_type_a}_Ancestral_Pops.tsv"))
analysis_b <- data.table::fread(glue::glue("Processed_Data/Admixture/{analysis_type_b}_Ancestral_Pops.tsv"))

# if analysis returned different optimal K, filter data sets to be equal for comparion
analysis_a <- analysis_a %>% dplyr::filter(K %in% analysis_b$K)
analysis_b <- analysis_b %>% dplyr::filter(K %in% analysis_a$K)


for(kpop in 1:length(unique(analysis_a$K))) {
  
  temp_a <- analysis_a %>% 
    dplyr::filter(K == unique(analysis_a$K)[kpop]) %>%
    tidyr::spread(cluster, frac_cluster)
  temp_b <- analysis_b %>% 
    dplyr::filter(K == unique(analysis_b$K)[kpop]) %>%
    tidyr::spread(cluster, frac_cluster)

  a_matrix <- as.matrix(temp_a[5:ncol(temp_a)])
  b_matrix <- as.matrix(temp_b[5:ncol(temp_b)])
  
  correl <- matcor(a_matrix, b_matrix)
  img.matcor(correl, type = 2)
  
  res.cc <- cc(a_matrix, b_matrix)
  
  barplot(res.cc$cor, xlab = "Dimension", ylab = "Canonical correlations", names.arg = 1:unique(analysis_a$K)[kpop], ylim = c(0,1))
  
  # red is a
  # blue is b
  png(filename=glue::glue("Plots/Admixture/Ancestry_Comparison_{analysis_type_a}_{analysis_type_b}_K{unique(analysis_a$K)[kpop]}.png"), width = 8, height = 8, res = 300, units = "in")
  plt.cc(res.cc, var.label = T, ind.names = temp_a$samples, type = "v")
  dev.off()
  
  
}


