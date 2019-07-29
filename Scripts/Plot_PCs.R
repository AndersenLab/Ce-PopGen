#!/usr/bin/env Rscript

library(tidyverse)

options(scipen=999)
source("Scripts/Figure_Themes.R")
#setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# args
# 1 - analysis type A - name of folder corresponding to what admixure analysis
# 2 - analysis type B - name of folder corresponding to what admixure analysis

# example - 
# args <- c("DIVERGENT-MASKED", "GATK-STRELKA_Intersection")
# args <- c("DIVERGENT-MASKED_Complete", "GATK-STRELKA_Intersection_Complete")
# args <- c("DIVERGENT-MASKED", "DIVERGENT-MASKED_Complete")
# args <- c("GATK-STRELKA_Intersection", "GATK-STRELKA_Intersection_Complete")

args <- commandArgs(TRUE)


dir.create("Plots/PCA")

analysis_type_a <- args[1]
analysis_type_b <- args[2]

A_n_remove <- glue::glue("Data/EIGENSTRAT/{analysis_type_a}/NO_REMOVAL/eigenstrat_no_removal.evac")
A_y_remove <- glue::glue("Data/EIGENSTRAT/{analysis_type_a}/OUTLIER_REMOVAL/eigenstrat_outliers_removed.evac")

B_n_remove <- glue::glue("Data/EIGENSTRAT/{analysis_type_b}/NO_REMOVAL/eigenstrat_no_removal.evac")
B_y_remove <- glue::glue("Data/EIGENSTRAT/{analysis_type_b}/OUTLIER_REMOVAL/eigenstrat_outliers_removed.evac")

compare_pca <- function(masked, intersection, removal = F) {
  
  pca_masked <- data.table::fread(masked, skip = 1, header = F) %>%
    dplyr::arrange(V1)
  
  n_pcs <- ncol(pca_masked)-1
  pca_masked <- pca_masked[,c(1:n_pcs)]
  colnames(pca_masked) <- c("Strain", paste0("PC", 1:(n_pcs-1)))
  
  pca_intersect <- data.table::fread(intersection, skip = 1, header = F) %>%
    dplyr::arrange(V1)
  
  n_pcs <- ncol(pca_intersect)-1
  pca_intersect <- pca_intersect[,c(1:n_pcs)]
  
  colnames(pca_intersect) <- c("Strain", paste0("PC", 1:(n_pcs-1)))
  
  if(removal == F) {
    pc_cor <- c()
    for(pcc in 2:(ncol(pca_intersect))){
      pc_cor <- append(x = pc_cor, cor(pca_masked[,pcc],pca_intersect[,pcc], method = "spearman"))
    }
    
    pc_cor <- data.frame(PC_COR = pc_cor, PC = paste0("PC", 1:(n_pcs-1)))
  } else {
    
    pca_masked <- pca_masked %>% dplyr::filter(Strain%in%pca_intersect$Strain)
    pca_intersect <- pca_intersect %>% dplyr::filter(Strain%in%pca_masked$Strain)
    
    pc_cor <- c()
    for(pcc in 2:51){
      pc_cor <- append(x = pc_cor, cor(pca_masked[,pcc],pca_intersect[,pcc], method = "spearman"))
    }
    
    pc_cor <- data.frame(PC_COR = pc_cor, PC = paste0("PC", 1:(n_pcs-1)))
    
  }
  
  return(pc_cor)
}

no_removal_pca_cor <- compare_pca(A_n_remove, B_n_remove)

ggplot(no_removal_pca_cor)+
  aes(x = factor(PC, levels  = PC), y = abs(PC_COR))+
  geom_point() +
  base_theme+
  theme(axis.text.x = element_text(angle = 90),
        axis.line = element_line())+
  labs(y = "Spearman's Correlation", x = "Principal Component", title = glue::glue("{analysis_type_a} - v - {analysis_type_b} \n NO REMOVAL"))

ggsave(filename = glue::glue("Plots/PCA/PCA_COR_{analysis_type_a}-v-{analysis_type_b}-NOREMOVAL.pdf"), height = 4, width = 12)

removal_pca_cor <- compare_pca(A_y_remove, B_y_remove, removal = T)

ggplot(removal_pca_cor)+
  aes(x = factor(PC, levels  = PC), y = abs(PC_COR))+
  geom_point() +
  base_theme+
  theme(axis.text.x = element_text(angle = 90),
        axis.line = element_line())+
  labs(y = "Spearman's Correlation", x = "Principal Component", title = glue::glue("{analysis_type_a} - v - {analysis_type_b} \n REMOVAL"))

ggsave(filename = glue::glue("Plots/PCA/PCA_COR_{analysis_type_a}-v-{analysis_type_b}-REMOVAL.pdf"), height = 4, width = 12)
