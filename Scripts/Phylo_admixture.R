#!/usr/bin/env Rscript
library(tidyverse)
library(ggtree)
library(ggstance)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

k = 14 
admix_pops <- data.table::fread("Processed_Data/Admixture/GATK-STRELKA_Intersection_Complete_Ancestral_Pops.tsv") %>%
  dplyr::filter(K== k)%>%
  dplyr::rename(strain=samples) %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::distinct()

btree <- ape::read.tree(glue::glue("Data/Phylogeny/GATK-STRELKA_Intersection_Complete/Population_Phylogeny.raxml.bestTree"))



# define a color pallette with maximal contrast
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "K"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "gray51","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")

# make long and determin order of plotting


branch_strains <- list()
for(i in 1:length(LETTERS[1:(k+1)])){
  if(i > k){
    branch_strains[[i]] <- admix_pops %>% dplyr::filter(frac_cluster < 0.9) %>% dplyr::pull(strain)
  } else {
    branch_strains[[i]] <- admix_pops %>% dplyr::filter(cluster == LETTERS[i], frac_cluster > 0.9) %>% dplyr::pull(strain)
  }
  
}

names(branch_strains) <- LETTERS[1:(k+1)]

ptree <- ggtree(tree_pt_h,
       branch.length="rate", 
       aes(color=group)) + 
  scale_color_manual(values=ancestry.colours) + 
  theme(legend.position="right")+
  geom_tiplab(align = T)+
  theme_tree2() + 
  xlim_tree(0.45)


admix_pops <- data.table::fread("Processed_Data/Admixture/GATK-STRELKA_Intersection_Complete_Ancestral_Pops.tsv") %>%
  dplyr::filter(K== k)

pw <- facet_plot(ptree, 
                 panel = 'Admixture', 
                 data = admix_pops, 
                 geom = geom_barh, 
           mapping = aes(x = frac_cluster, fill = cluster),
           stat='identity', 
           color = NA) + 
  scale_fill_manual(values = ancestry.colours)
pw

ggsave(pw, filename = glue::glue("Plots/Phylogeny/Phylo_Admixture_K{k}.pdf"),height = 48, width = 24)

