#!/usr/bin/env Rscript

library(pophelper)
library(tidyverse)
library(ggthemes)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# args
# 1 - analysis type - name of folder corresponding to what admixure analysis
# 2 - Outgroup

# example - 
# args <- c("DIVERGENT-MASKED", "XZ1516")
# args <- c("GATK-STRELKA_Intersection", "XZ1516")
# args <- c("DIVERGENT-MASKED_Complete", "XZ1516")
# args <- c("GATK-STRELKA_Intersection_Complete", "XZ1516")
args <- commandArgs(TRUE)
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  Generate Population Summary File for TREEMIX  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

source("scripts/PLOT_TREEMIX.R")

analysis_type <- args[1]
# define outgroup
outgroup_strain <- args[2]

data_dir <- glue::glue("Data/ADMIXTURE/{analysis_type}/")

dir.create(glue::glue("{data_dir}/Treemix"))
dir.create(glue::glue("Plots/Treemix/{analysis_type}"))

# strain names
# get sample information
sample_names <- data.table::fread(glue::glue("data/{analysis_type}_samples.txt"), header = F) %>% dplyr::pull(V1)



for(kpops in 1:length(grep(".Q", list.files(glue::glue("{data_dir}")), value = T))){ 
  
  K <- as.numeric(strsplit(grep(".Q", list.files(glue::glue("{data_dir}")), value = T)[kpops], split = "\\.")[[1]][4])
  
  # load P files
  pfile_name <- grep(pattern = glue::glue("{K}\\.P$"), value = T, x = list.files(glue::glue("{data_dir}")))
  pfile <- pophelper::readQ(files = paste0(glue::glue("{data_dir}"),pfile_name))[[1]]
  
  # label P file rownames and colnames
  colnames(pfile) <- LETTERS[1:K]
  
  
  treemix_input <- apply(pfile, MARGIN = c(1,2), function(x){
    f <- round(x*length(sample_names),digits = 0)
    paste(length(sample_names)-f, f, sep = ",")
  })
  
  # load Q files
  qfile_name <- grep(pattern = glue::glue("{K}\\.Q$"), value = T, x = list.files(glue::glue("{data_dir}")))
  qfile <- pophelper::readQ(files = paste0(glue::glue("{data_dir}"),qfile_name))[[1]]
  # add pop names
  colnames(qfile) <- LETTERS[1:K]
  qfile$strain <- sample_names
  
  # find outgroup population
  outgroup_population <- qfile%>%
    dplyr::filter(strain == outgroup_strain)%>%
    tidyr::gather(Population, Frequency, -strain)%>%
    dplyr::filter(Frequency == max(Frequency))
  
  write.table(x = treemix_input,
              file = glue::glue("{data_dir}Treemix/K-{K}_Outgroup={outgroup_population$Population[1]}=TREEMIX_input.txt"),
              quote = F,
              col.names = T,
              row.names = F)
  
  treemix_input_name <- strsplit(glue::glue("{data_dir}Treemix/K-{K}_Outgroup={outgroup_population$Population[1]}=TREEMIX_input.txt"),
                                 split = "/")[[1]][5]
  
  system(glue::glue("gzip {data_dir}Treemix/{treemix_input_name}"))
  
  outgroup_pop <- outgroup_population$Population[1]
  
  output_name <- strsplit(strsplit(glue::glue("{data_dir}Treemix/K-{K}_Outgroup={outgroup_population$Population[1]}=TREEMIX_input.txt"),
                                   split = "/")[[1]][5], split = "\\.txt")[[1]][1]
  
  system(glue::glue("treemix -i {data_dir}Treemix/{treemix_input_name}.gz -o {data_dir}Treemix/{output_name} -root {outgroup_pop} -k 500"))
  
  system(glue::glue("threepop -i {data_dir}Treemix/{treemix_input_name}.gz -k 500 | grep ';' > {data_dir}Treemix/{treemix_input_name}_threepop" ))
  
  threepop <- data.table::fread(glue::glue("{data_dir}Treemix/{treemix_input_name}_threepop"), sep = " ") %>%
    dplyr::select(tree = V1, F3 = V2, SE = V3, Zscore = V4)
  
  system(glue::glue("fourpop -i {data_dir}Treemix/{treemix_input_name}.gz -k 500 | grep ';' > {data_dir}Treemix/{treemix_input_name}_fourpop" ))
  
  fourpop <- data.table::fread(glue::glue("{data_dir}Treemix/{treemix_input_name}_fourpop"), sep = " ") %>%
    dplyr::select(tree = V1, F4 = V2, SE = V3, Zscore = V4)
  
  
  system(paste0(paste0("printf \"", paste(LETTERS[1:K], collapse = '\n')), '\n\" > ', glue::glue("{data_dir}Treemix/poporder")))
  
  treemix_stem <- glue::glue("{data_dir}Treemix/{output_name}")
  
  pdf(glue::glue("Plots/Treemix/{analysis_type}/TREEMIX_phylogeny_K{K}_NO_MIGRATION.pdf"))
  plot_tree(stem = treemix_stem)
  dev.off()
  
  pdf(glue::glue("Plots/Treemix/{analysis_type}/TREEMIX_phylogeny_K{K}_NO_MIGRATION_Residuals.pdf"))
  plot_resid(stem = treemix_stem, glue::glue("{data_dir}Treemix/poporder"))
  dev.off()
  
  # treemix with migration
  
  # run treemix with migration -
  dir.create(glue::glue("{data_dir}Treemix/Migrations"))
  
  for(m in 1:4){
    system(glue::glue("treemix -i {data_dir}Treemix/{treemix_input_name}.gz -o {data_dir}Treemix/Migrations/{output_name}_{m} -root {outgroup_pop} -k 500 -m {m} -se"))
  }
  
  # # # # #  PLOT TREEMIX WITH MIGRATIONS
  
  for(migration in gsub("\\.llik","",grep(glue::glue("K-{K}"), grep("llik",list.files(glue::glue("{data_dir}Treemix/Migrations")),value = T), value = T))){
    
    system(paste0(paste0("printf \"", paste(LETTERS[1:K], collapse = '\n')), '\n\" > ', glue::glue("{data_dir}Treemix/Migrations/poporder")))
    
    pdf(glue::glue("Plots/Treemix/{analysis_type}/TREEMIX_phylogeny_K{K}_{migration}_migrations.pdf"))
    plot_tree(stem =  glue::glue("{data_dir}Treemix/Migrations/{migration}"))
    dev.off()
    
    pdf(glue::glue("Plots/Treemix/{analysis_type}/TREEMIX_phylogeny_K{K}_{migration}_migrations_RESIDUALS.pdf"))
    plot_resid(stem = glue::glue("{data_dir}Treemix/Migrations/{migration}"), glue::glue("{data_dir}Treemix/poporder"))
    dev.off()
    # 
  }
  
}






