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

data_dir <- glue::glue("Data/POPGENOME/{analysis_type}/SUBPOPs/WINDOW/")

# define a color pallette with maximal contrast
# Colors from - https://gist.github.com/ollieglass/f6ddd781eeae1d24e391265432297538
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "K"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")



load(glue::glue("{data_dir}Ce_Genome-wide_Neutrality_stats.Rda")) 

# Set k and load
kpop = 15
focal_pops <- c("B","D","E","F","L","O","N")


nd_df <- neutrality_df %>%
  dplyr::filter(K == kpop)


plot_st <- c("Fu.F_S")
nd_df %>%
  dplyr::filter(Population %in% c("D", "N"), statistic == plot_st) %>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = value, color = Population) +
  geom_line()+
  scale_color_manual(values = ancestry.colours) +
  facet_grid(CHROM~., scales = "free", space = "free")+
  theme_bw(18)+
  labs(x = "Genomic Position (Mb)")

plot_st <- c("Fu.F_S")
nd_df %>%
  dplyr::filter(Population %in% c("E", "O"), statistic == plot_st) %>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = value, color = Population) +
  geom_line()+
  scale_color_manual(values = ancestry.colours) +
  facet_grid(CHROM~statistic, scales = "free", space = "free")+
  theme_bw(18)+
  labs(x = "Genomic Position (Mb)")
