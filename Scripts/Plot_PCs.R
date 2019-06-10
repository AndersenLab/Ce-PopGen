library(tidyverse)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))





pc_masked <- data.table::fread("~/transfer/popgen/Masked_eigenstrat_no_removal.evac",skip = 1, header = F) %>%
  dplyr::arrange(V1)

colnames(pc_masked) <- c("Strain", paste0("PC", 1:(ncol(pc_masked)-1))) 


ggplot(pc_masked)+
  aes(x = PC4, y = PC5)+
  geom_point()+
  theme_classic(18)


pcs <- data.table::fread("~/transfer/popgen/eigenstrat_no_removal.evac",skip = 1, header = F)%>%
  dplyr::arrange(V1)

colnames(pcs) <- c("Strain", paste0("PC", 1:(ncol(pcs)-1))) 


ggplot(pcs)+
  aes(x = PC4, y = PC5)+
  geom_point()+
  theme_classic(18)


pc_cor <- c()
for(pcc in 2:51){
  pc_cor <- append(x = pc_cor,cor(pcs[,pcc],pc_masked[,pcc], method = "spearman"))
  }

pc_cordf <- data.frame(pc_cor = pc_cor, pc = 1:length(pc_cor))

ggplot(pc_cordf)+
  aes(x = pc, y = abs(pc_cor))+
  geom_point()

pc_masked <- data.table::fread("~/transfer/popgen/Masked_eigenstrat_outliers_removed.evac",skip = 1, header = F) %>%
  dplyr::arrange(V1)

colnames(pc_masked) <- c("Strain", paste0("PC", 1:(ncol(pc_masked)-1))) 


ggplot(pc_masked)+
  aes(x = PC3, y = PC4)+
  geom_point()+
  theme_classic(18)


pcs <- data.table::fread("~/transfer/popgen/eigenstrat_outliers_removed.evac",skip = 1, header = F)%>%
  dplyr::arrange(V1)

colnames(pcs) <- c("Strain", paste0("PC", 1:(ncol(pcs)-1))) 


ggplot(pcs)+
  aes(x = PC3, y = PC4)+
  geom_point()+
  theme_classic(18)

pc_match <- pcs %>% dplyr::filter(Strain%in%pc_masked$Strain)
pcm_match <- pc_masked %>% dplyr::filter(Strain%in%pcs$Strain)

pc_cor <- c()
for(pcc in 2:51){
  pc_cor <- append(x = pc_cor,cor(pcm_match[,pcc],pc_match[,pcc], method = "spearman"))
}

pc_cordf <- data.frame(pc_cor = pc_cor, pc = 1:length(pc_cor))

ggplot(pc_cordf)+
  aes(x = pc, y = abs(pc_cor))+
  geom_point()
