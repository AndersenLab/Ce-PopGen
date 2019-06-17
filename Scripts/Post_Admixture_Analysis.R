library(tidyverse)
library(clusterProfiler)
library(org.Ce.eg.db)
library(biomaRt)

options(scipen=999)

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

# load admixture pop assignments
for(f in list.files(glue::glue("Processed_Data/Admixture/"))){
  temp_res <- data.table::fread(glue::glue("Processed_Data/ADMIXTURE/{f}"))
  if(!exists("admix_res")){
    admix_res <- temp_res
  } else {
    admix_res <- dplyr::bind_rows(admix_res, temp_res)
  }
  
}

# load mask data
masked_directory <- "Data/Divergent_Masks/Processed_Masks/"
mask_files <- list.files(masked_directory)

# Extract masked regions and retain type of mask
for(sm in 1:length(mask_files)){
  temp_sm_mask <- data.table::fread(glue::glue("{masked_directory}{mask_files[sm]}")) %>%
    dplyr::select(CHROM, START_BIN, END_BIN, STRAIN, window_mask) %>%
    dplyr::filter(window_mask != "Pass")
  
  if(!exists("population_masks")){
    population_masks <- temp_sm_mask
  } else {
    population_masks <- dplyr::bind_rows(population_masks, temp_sm_mask)
  }
}

all_masks <- population_masks %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(frq = n()/n_sm) 

ggplot(all_masks) +
  aes(x = END_BIN/1e6, y = frq)+
  geom_point()+
  facet_grid(CHROM~.) +
  theme_bw(18) +
  labs(x = "Genomic Position (Mb)", y = "Frequency")

ggsave("Plots/Enrichment/Mask_Frequency.pdf", height = 10, width = 12)

k <- 7
atype <- "DIVERGENT-MASKED_Complete"

mask_by_pop <- admix_res %>%
  dplyr::filter(K==k) %>%
  dplyr::filter(frac_cluster == max_frac, analysis_type == atype) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(n_st = n()) %>%
  dplyr::ungroup() %>%
  dplyr::rename(STRAIN = samples) %>%
  dplyr::left_join(population_masks, ., by="STRAIN") %>%
  dplyr::group_by(cluster, analysis_type, CHROM, START_BIN, END_BIN, n_st) %>%
  dplyr::mutate(mask_ct = n()) %>%
  dplyr::mutate(frq = mask_ct/n_st) %>%
  dplyr::ungroup() %>%
  # dplyr::mutate(frq = ifelse(frq > 0.5, 1-frq, frq)) %>%
  dplyr::distinct(CHROM, START_BIN, END_BIN, cluster, frq)

# mask_by_pop%>%
#   dplyr::filter(cluster == "A") %>%
  ggplot(mask_by_pop) +
  aes(x = START_BIN, y = frq, fill = cluster) +
  geom_bar(stat='identity') +
  facet_grid(CHROM~.) +
  theme_bw()
