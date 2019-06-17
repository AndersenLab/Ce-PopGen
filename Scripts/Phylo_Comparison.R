library(tidyverse)
library(ggtree)
library(ape)
library(phytools)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

analysis_type_a <- "GATK-STRELKA_Intersection"
analysis_type_b <- "DIVERGENT-MASKED"

tree_a <- ape::read.tree(glue::glue("Data/Phylogeny/{analysis_type_a}/Population_Phylogeny.raxml.bestTree"))
tree_b <- ape::read.tree(glue::glue("Data/Phylogeny/{analysis_type_b}/Population_Phylogeny.raxml.bestTree"))


plot_rep_strain_phylo <- function(tree, strains = c("JU2825", "JU3137", "MY23", "ECA928", "ECA746", "MY518", "NIC515","JU3127", "ECA347","N2","CB4856","ECA191","JU775","WN2001")) {
  # strains of interest
  strains_to_plot <- strains
  
  # define colors
  highlight_color <- "#D7263D"
  background_color <- "#000F08"
  
  # highlight branches for strains of interest
  branch_strains <- list(CONNECT = strains_to_plot)
  
  tree_pt_h <- ggtree::groupOTU(tree, branch_strains)
  
  highlight_tree <- ggtree(tree_pt_h,
         branch.length="rate", 
         aes(color=group)) + 
    geom_tiplab(align = T) +
    scale_color_manual(values=c(background_color, highlight_color)) + 
    theme(legend.position="right")+
    theme_tree2() 
  return(highlight_tree)
}

plot_rep_strain_phylo(tree_a) +xlim(c(0, 0.35))
ggsave(glue::glue("Plots/Phylogeny/{analysis_type_a}_Phylogeny.pdf"), height = 42, width = 12)
plot_rep_strain_phylo(tree_b) +xlim(c(0, 0.35))
ggsave(glue::glue("Plots/Phylogeny/{analysis_type_b}_Phylogeny.pdf"), height = 42, width = 12)


# ta <- ggtree(tree_a)
# tb <- ggtree(tree_b)
# 
# d1 <- ta$data
# d2 <- tb$data
# 
# d2$x <- max(d2$x) - d2$x + max(d1$x) + 1
# 
# pp <- ta + geom_tiplab() + geom_tree(data=d2) + geom_tiplab(data = d2, hjust=1)
# 
# dd <- bind_rows(d1, d2) %>% 
#   filter(!is.na(label))
# 
# pp + geom_line(aes(x, y, group=label), data=dd, color='grey')


