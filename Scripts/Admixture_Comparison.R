#!/usr/bin/env Rscript
install.packages("remotes")
remotes::install_github("gavinsimpson/ggvegan")
library(tidyverse)
library(CCA)
library(ggvegan)

source("Scripts/Figure_Themes.R")
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# analysis_type_a <- "GATK-STRELKA_Intersection"
# analysis_type_b <- "DIVERGENT-MASKED"

analysis_type_a <- "GATK-STRELKA_Intersection_Complete"
analysis_type_b <- "DIVERGENT-MASKED_Complete"

analysis_a <- data.table::fread(glue::glue("Processed_Data/Admixture/{analysis_type_a}_Ancestral_Pops.tsv"))
analysis_b <- data.table::fread(glue::glue("Processed_Data/Admixture/{analysis_type_b}_Ancestral_Pops.tsv"))

# if analysis returned different optimal K, filter data sets to be equal for comparion
analysis_a <- analysis_a %>% dplyr::filter(K %in% analysis_b$K)
analysis_b <- analysis_b %>% dplyr::filter(K %in% analysis_a$K)

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

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
  a_cor <- data.frame(Pop = res.cc$names$Xnames,
                      vcf_type = analysis_type_a,
                      d1 = res.cc$scores$corr.X.xscores[,1],
                      d2 = res.cc$scores$corr.X.xscores[,2])
  
  b_cor <- data.frame(Pop = res.cc$names$Ynames,
                      vcf_type = analysis_type_b,
                      d1 = res.cc$scores$corr.Y.xscores[,1],
                      d2 = res.cc$scores$corr.Y.xscores[,2])
  
  pltdf <- dplyr::bind_rows(a_cor,b_cor) 
  
  ggplot(pltdf)+
    base_theme+
    geom_path(aes(x,y),data=circleFun(c(0,0),2,npoints = 100), color = "gray60") +
    geom_path(aes(x,y),data=circleFun(c(0,0),1,npoints = 100), color = "gray60") +
    geom_text(aes(x = d1, y = d2, color = vcf_type, label = Pop), size = 8, fontface = "bold") +
    scale_color_manual(values = c("blue","red"),
                       name="VCF",
                       breaks=c(analysis_type_a, analysis_type_b),
                       labels=c("Intersection", "Masked"))+
    labs(x= "Dimension 1", y = "Dimension 2") +
    theme(legend.position = "top")
  
  ggsave(glue::glue("Plots/Admixture/Ancestry_Comparison_{analysis_type_a}_{analysis_type_b}_K{unique(analysis_a$K)[kpop]}.pdf"),
         width = 8, height = 8)
}








  

