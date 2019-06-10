library(tidyverse)
library(clusterProfiler)
library(org.Ce.eg.db)

options(scipen=999)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

masked_n_remove <- "Data/EIGENSTRAT/DIVERGENT-MASKED/NO_REMOVAL/eigenstrat_no_removal.evac"
masked_y_remove <- "Data/EIGENSTRAT/DIVERGENT-MASKED/OUTLIER_REMOVAL/eigenstrat_outliers_removed.evac"

intersect_n_remove <- "Data/EIGENSTRAT/GATK-STRELKA_Intersection/NO_REMOVAL/eigenstrat_no_removal.evac"
intersect_y_remove <- "Data/EIGENSTRAT/GATK-STRELKA_Intersection/OUTLIER_REMOVAL/eigenstrat_outliers_removed.evac"

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

no_removal_pca_cor <- compare_pca(masked_n_remove, intersect_n_remove)

ggplot(no_removal_pca_cor)+
  aes(x = factor(PC, levels  = PC), y = abs(PC_COR))+
  geom_point() +
  theme_classic(18) +
  theme(axis.text.x = element_text(angle = 90))+
  labs(y = "Spearman's Correlation", x = "Principal Component")


removal_pca_cor <- compare_pca(masked_y_remove, intersect_y_remove, removal = T)

ggplot(removal_pca_cor)+
  aes(x = factor(PC, levels  = PC), y = abs(PC_COR))+
  geom_point() +
  theme_classic(18) +
  theme(axis.text.x = element_text(angle = 90))+
  labs(y = "Spearman's Correlation", x = "Principal Component")


