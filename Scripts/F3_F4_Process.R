#!/usr/bin/env Rscript

library(pophelper)
library(tidyverse)
library(ggthemes)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# run script in base Ce-PopGen directory

# args
# 1 - analysis type - name of folder corresponding to what admixure analysis

# example - 
# args <- c("GATK-STRELKA_Intersection_Complete", "14")
args <- commandArgs(TRUE)

analysis_type <- args[1]
analysis_k <- args[2]

data_dir <- glue::glue("Data/ADMIXTURE/{analysis_type}/Treemix")

# define a color pallette with maximal contrast
# Colors from - https://gist.github.com/ollieglass/f6ddd781eeae1d24e391265432297538
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "K"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")

list.files(glue::glue("{data_dir}"))
K-{analysis_k}*
grep("threepop", )

fthree <- data.table::fread(glue::glue("{data_dir}K-{analysis_k}"),sep = " ")%>%
  tidyr::separate(V1, into = c("Outgroup","otherpops"), sep = ";")%>%
  dplyr::rename(F3 = V2, SE = V3, Z = V4)%>%
  dplyr::group_by(Outgroup)%>%
  dplyr::mutate(mf3 = median(F3))%>%
  dplyr::ungroup()%>%
  dplyr::arrange(desc(mf3))%>%
  dplyr::mutate(Outgroup2 = factor(Outgroup,levels = unique(Outgroup), labels = unique(Outgroup)))

ggplot(fthree)+
  aes(x = Outgroup2, y = F3, fill = Outgroup)+
  geom_boxplot(alpha = 0.85)+
  theme_bw()+
  scale_fill_manual(values=ancestry.colours)+
  labs(x = "Outgroup", y = "F3 Statistic")+
  theme(axis.text.x = ggplot2::element_text(size = 16),
        axis.text.y = ggplot2::element_text(size = 16),
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3))
