#!/usr/bin/env Rscript

library(pophelper)
library(tidyverse)
library(ggthemes)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# run script in base Ce-PopGen directory

# args
# 1 - analysis type - name of folder corresponding to what admixure analysis

# example - 
# args <- c("GATK-STRELKA_Intersection_Complete")
args <- commandArgs(TRUE)

analysis_type <- args[1]

data_dir <- glue::glue("Data/POPGENOME/{analysis_type}/SUBPOPs/GENE/")

# define a color pallette with maximal contrast
# Colors from - https://gist.github.com/ollieglass/f6ddd781eeae1d24e391265432297538
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "K"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")
# chnage names tp match FST df pop names
ancestry.colours2 <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                       "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                       "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "K"="black", 
                       "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                       "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")
names(ancestry.colours2) <- paste0("pop",c(1:length(ancestry.colours2)))

# extract chromosomes
chroms <- grep('[A-Z]', list.dirs(data_dir, full.names = F), value = T)


for(chrom in chroms){
  FST_files <- grep("FST",list.files(glue::glue("{data_dir}{chrom}")), value=T)
  
  for(k in FST_files) {
    
    load(glue::glue("{data_dir}{chrom}/{k}"))
    
    temp_fst <- pairFst %>%
      dplyr::mutate(CHROM = chrom,
                    K = strsplit(k, split = "_")[[1]][2])
    
    if(!exists("fst_df")) {
      fst_df <- temp_fst
    } else {
      fst_df <- dplyr::bind_rows(fst_df, temp_fst)
    }
  }
}

fst_df%>%
  dplyr::filter(K == "K14", Fst > 0.1)%>%
  ggplot()+
  aes(x = Gene_Start, y = Fst)+
  geom_point()+
  facet_grid(.~CHROM)

