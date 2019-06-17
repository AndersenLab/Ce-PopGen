#!/usr/bin/env Rscript

library(pophelper)
library(tidyverse)
library(ggthemes)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# run script in base Ce-PopGen directory

# args
# 1 - analysis type - name of folder corresponding to what admixure analysis

# example - 
# args <- c("DIVERGENT-MASKED")
# args <- c("GATK-STRELKA_Intersection")
# args <- c("DIVERGENT-MASKED_Complete")
# args <- c("GATK-STRELKA_Intersection_Complete")
args <- commandArgs(TRUE)

analysis_type <- args[1]

data_dir <- glue::glue("Data/ADMIXTURE/{analysis_type}/")

strain_islands <- c("XZ1514" = "#E69F00", "XZ1516" = "#E69F00","XZ1513" = "#E69F00","ECA372" = "#E69F00","ECA701" = "#E69F00","XZ1515" = "#E69F00",
                    "CB4856" = "#56B4E9",
                    "ECA369" = "#009E73","ECA738" = "#009E73",
                    "QX1792" = "#0072B2", "QX1794" = "#0072B2", "QX1793" = "#0072B2", "QX1791" = "#0072B2", "ECA740" = "#0072B2", "ECA741" = "#0072B2", "ECA363" = "#0072B2", "ECA743" = "#0072B2", "ECA742" = "#0072B2",
                    "ECA760" = "#CC79A7","ECA768" = "#CC79A7","ECA777" = "#CC79A7","ECA706" = "#CC79A7","ECA705" = "#CC79A7","ECA703" = "#CC79A7","ECA807" = "#CC79A7","ECA778" = "#CC79A7",
                    "ECA812" = "#CC79A7","ECA710" = "#CC79A7","ECA744" = "#CC79A7","ECA745" = "#CC79A7","ECA732" = "#CC79A7","ECA733" = "#CC79A7","ECA746" = "#CC79A7","DL238" = "#CC79A7",
                    "ECA347" = "#CC79A7","ECA730" = "#CC79A7","ECA724" = "#CC79A7","ECA722" = "#CC79A7","ECA189" = "#CC79A7","ECA191" = "#CC79A7","ECA723" = "#CC79A7","ECA712" = "#CC79A7",
                    "ECA396" = "#CC79A7")

# define a color pallette with maximal contrast
# Colors from - https://gist.github.com/ollieglass/f6ddd781eeae1d24e391265432297538
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "WHOLE_POPULATION"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")

# strain names
# get sample information
sample_names <- data.table::fread(glue::glue("data/{analysis_type}_samples.txt"), header = F) %>% dplyr::pull(V1)

# plot panel A, K by CV summary
k_summary <- data.table::fread(glue::glue("{data_dir}CV_Summary/admix_replicates_CV.tsv") ,header = T) 

ksum_plot <- ggplot(k_summary)+
  aes(x = factor(K), y = CV)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = .1)+
  theme_bw()+
  labs(x = "K", title = analysis_type)


# generate K summary plot - Supplemental figure XX
admix_plots <- list()
admix_plots_with_names <- list()

for(kpops in 1:length(grep(".Q", list.files(glue::glue("{data_dir}")), value = T))){
  K <- as.numeric(strsplit(grep(".Q", list.files(glue::glue("{data_dir}")), value = T)[kpops], split = "\\.")[[1]][4])
  
  # load Q files
  qfile_name <- grep(pattern = glue::glue("{K}\\.Q$"), value = T, x = list.files(glue::glue("{data_dir}")))
  qfile <- pophelper::readQ(files = paste0(glue::glue("{data_dir}"),qfile_name))[[1]]
  # add pop names
  colnames(qfile) <- LETTERS[1:K]
  
  # make long and determin order of plotting
  long_admix_pops <- qfile %>%
    dplyr::mutate(samples = sample_names) %>%
    tidyr::gather(cluster, frac_cluster, -samples) %>%
    dplyr::group_by(samples) %>%
    dplyr::mutate(max_frac = max(frac_cluster)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)),
                  analysis_type = analysis_type,
                  K = K)
  
  # establish plot order of strains based on anc pop and max fraction
  plot_order <- long_admix_pops %>%
    dplyr::filter(frac_cluster == max_frac) %>%
    dplyr::arrange(cluster, -max_frac) %>%
    dplyr::mutate(samples = factor(samples, levels = unique(samples)))
  
  admix_plots[[kpops]] <- long_admix_pops %>%
    dplyr::mutate(ordered_samples = factor(samples, levels = plot_order$samples)) %>%
    ggplot() +
    geom_bar(stat = "identity", 
             aes(x = ordered_samples, 
                 y = frac_cluster, 
                 fill = cluster)) +
    scale_fill_manual(values = ancestry.colours) +
    labs(fill = "", x = "", y =  glue::glue("K = {K}")) +
    theme_bw() +
    theme(axis.text.x=element_blank(),    
          axis.text.y=element_blank(),
          axis.title.y = element_text(angle = 0, vjust = .5),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "cm"))
  
  admix_plots_with_names[[kpops]] <- long_admix_pops %>%
    dplyr::mutate(ordered_samples = factor(samples, levels = plot_order$samples)) %>%
    ggplot() +
    geom_bar(stat = "identity", 
             aes(x = ordered_samples, 
                 y = frac_cluster, 
                 fill = cluster)) +
    scale_fill_manual(values = ancestry.colours) +
    labs(fill = "", x = "", y =  glue::glue("K = {K}")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 90, hjust = 1),    
          axis.text.y=element_blank(),
          axis.title.y = element_text(angle = 0, vjust = .5),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          plot.margin = unit(c(0,0,0,0), units = "cm"))
  
  if(!exists("final_admix_groups")){
    final_admix_groups <- long_admix_pops
  } else{
    final_admix_groups <- dplyr::bind_rows(final_admix_groups, long_admix_pops)
  }
  
}

# make panel B
admixture_plots <- cowplot::plot_grid(admix_plots[[1]] + theme(legend.position = "none"),
                                      admix_plots[[2]] + theme(legend.position = "none"),
                                      admix_plots[[3]] + theme(legend.position = "none"),
                                      admix_plots[[4]] + theme(legend.position = "none"),
                                      admix_plots[[5]] + theme(legend.position = "none"),
                                      ncol = 1)

# make final figure
ksummary_plot <- cowplot::plot_grid(ksum_plot,
                                    admixture_plots,
                                    ncol = 2,
                                    labels = c("A", "B"),
                                    rel_widths = c(0.5, 1))

# get big legend
admix_legend <- cowplot::plot_grid(cowplot::get_legend(admix_plots[[5]]))

# save
ggsave(ksummary_plot, filename = glue::glue("Plots/Admixture/{analysis_type}_KSummary.pdf"), height = 8, width = 12)
ggsave(ksummary_plot, filename = glue::glue("Plots/Admixture/{analysis_type}_KSummary.png"), height = 8, width = 12,dpi=300)
ggsave(admix_legend, filename = glue::glue("Plots/Admixture/{analysis_type}_admix_legend.pdf"), height = 8, width = 12)

write.table(final_admix_groups, file = glue::glue("Processed_Data/Admixture/{analysis_type}_Ancestral_Pops.tsv"), 
            col.names = T, row.names = F, quote = F, sep = "\t")

for(pk in 1:length(unique(final_admix_groups$K))){
  
  ggsave(admix_plots_with_names[[pk]] + labs(title = analysis_type), 
         filename = glue::glue("Plots/Admixture/{analysis_type}_Ancestor_Frequencies_K{unique(final_admix_groups$K)[pk]}.pdf"),
         height = 6, width = 36)
  
}

