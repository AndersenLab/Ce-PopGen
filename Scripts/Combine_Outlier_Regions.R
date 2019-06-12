library(tidyverse)
library(clusterProfiler)
library(org.Ce.eg.db)


options(scipen=999)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

pr_mask_directory <- "Data/Divergent_Masks/Masks_By_Type/"

dir.create(pr_mask_directory)

masked_directory <- "Data/Divergent_Masks/Processed_Masks/"

mask_files <- list.files(masked_directory)

n_sm <- length(mask_files)

# mask frequency thresholds
mask_frq_thresh <- c(0.01, 0.05, 0.5)

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

# get mask frequency
all_masks <- population_masks %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(frq = n()/n_sm) %>%
  dplyr::mutate(frq = ifelse(frq > 0.5, 1-frq, frq))

# generate bed files for each frequency class
# low
low_freq_mask <- dplyr::filter(all_masks, frq < mask_frq_thresh[1]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(START_BIN = ifelse(START_BIN == 0, 1, START_BIN))

write.table(low_freq_mask, 
            file = glue::glue("{pr_mask_directory}Low_Freq_Masks_All_Classes.bed"), 
            col.names = F, row.names = F, quote = F, sep = "\t")

# intermediate
int_freq_mask <- dplyr::filter(all_masks, frq > mask_frq_thresh[1]) %>%
  dplyr::filter(frq < mask_frq_thresh[2]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(START_BIN = ifelse(START_BIN == 0, 1, START_BIN))

write.table(int_freq_mask, 
            file = glue::glue("{pr_mask_directory}Intermediate_Freq_Masks_All_Classes.bed"), 
            col.names = F, row.names = F, quote = F, sep = "\t")

# common
com_freq_mask <- dplyr::filter(all_masks, frq > mask_frq_thresh[2]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(START_BIN = ifelse(START_BIN == 0, 1, START_BIN))

write.table(com_freq_mask, 
            file = glue::glue("{pr_mask_directory}Common_Freq_Masks_All_Classes.bed"), 
            col.names = F, row.names = F, quote = F, sep = "\t")

# find genes in masked regions
system("Scripts/Masked_Genes.sh")

# load gene files
gene_files <- grep("Genes", list.files(pr_mask_directory), value = T)

for(freq_cut in 1:length(gene_files)){
  
  cutoff <- strsplit(gene_files[freq_cut], split = "_")[[1]][1]
  
  temp_genes <- data.table::fread(glue::glue("{pr_mask_directory}{gene_files[freq_cut]}")) %>%
    dplyr::select(WBGeneID = V8) %>%
    dplyr::mutate(FreqCutoff = cutoff)
  
  if(!exists("masked_genes")){
    masked_genes <- temp_genes
  } else {
    masked_genes <- dplyr::bind_rows(masked_genes, temp_genes)
  }
}

# perform enrichment analysis
enrich_go <- function(x = NULL){
  mf <- enrichGO(gene          = x,
           OrgDb         = org.Ce.eg.db,
           ont           = "MF",
           pAdjustMethod = "BH",
           keyType       = 'ENSEMBL',
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05)
  
  bp <- enrichGO(gene          = x,
           OrgDb         = org.Ce.eg.db,
           ont           = "BP",
           pAdjustMethod = "BH",
           keyType       = 'ENSEMBL',
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05)
  return(list(mf, bp))
}

common_cds_gene_MF <- enrich_go(x = dplyr::filter(masked_genes, FreqCutoff == "Common") %>% dplyr::pull(WBGeneID))
intermediate_cds_gene_MF <- enrich_go(x = dplyr::filter(masked_genes, FreqCutoff == "Intermediate") %>% dplyr::pull(WBGeneID))
low_cds_gene_MF <- enrich_go(x = dplyr::filter(masked_genes, FreqCutoff == "Low") %>% dplyr::pull(WBGeneID))

barplot(common_cds_gene_MF[[1]])
barplot(common_cds_gene_MF[[2]])

goplot(common_cds_gene_MF[[1]])
goplot(common_cds_gene_MF[[2]])

barplot(intermediate_cds_gene_MF[[1]])
barplot(intermediate_cds_gene_MF[[2]])

barplot(low_cds_gene_MF[[1]])
barplot(low_cds_gene_MF[[2]])



goplot(intermediate_cds_gene_MF[[1]])
goplot(intermediate_cds_gene_MF[[2]])

goplot(low_cds_gene_MF[[1]])
goplot(low_cds_gene_MF[[2]])

