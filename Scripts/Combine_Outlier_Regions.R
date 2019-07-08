library(tidyverse)
library(clusterProfiler)
library(org.Ce.eg.db)
library(biomaRt)

options(scipen=999)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

pr_mask_directory <- "Data/Divergent_Masks/Masks_By_Type/"
dir.create(pr_mask_directory)

masked_directory <- "Data/Divergent_Masks/Processed_Masks/"

join_masks <- function(mask_file = NULL, jump_tolerance = 1000){
  common_cluster <- NA
  common_cluster_start <- NA
  
  for (i in c(1:nrow(mask_file)-1)) {
    
    if(i==1){
      
      common_cluster[i] <- ifelse(mask_file$END_BIN[i] == mask_file$START_BIN[i+1], 'yes', 'no')
      common_cluster_start[i] <- mask_file$START_BIN[i]
      
    } else {
      
      common_cluster[i] <- ifelse(mask_file$END_BIN[i-1] > mask_file$START_BIN[i] - jump_tolerance | mask_file$END_BIN[i] > mask_file$START_BIN[i+1] - jump_tolerance, 'yes', 'no')
      common_cluster_start[i] <- ifelse(mask_file$CHROM[i] == mask_file$CHROM[i-1] & common_cluster[i] == 'yes', 
                                        ifelse(mask_file$END_BIN[i-1] > mask_file$START_BIN[i] - jump_tolerance | mask_file$END_BIN[i] > mask_file$START_BIN[i+1] - jump_tolerance, 
                                               ifelse(common_cluster[i-1] == "no", mask_file$START_BIN[i], 
                                                      ifelse(mask_file$END_BIN[i-1] == mask_file$START_BIN[i], common_cluster_start[i-1], mask_file$START_BIN[i])), mask_file$START_BIN[i]), mask_file$START_BIN[i])
      
    }
  }
  
  common_cluster[nrow(mask_file)] <- NA
  common_cluster_start[nrow(mask_file)] <- NA
  
  mask_file_cluster <- data.frame(mask_file, common_cluster, common_cluster_start) %>%
    dplyr::group_by(CHROM, common_cluster_start) %>%
    dplyr::mutate(common_cluster_size=ifelse(is.na(common_cluster_start), 1000, n()*1000)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(common_cluster_end=common_cluster_start+common_cluster_size)
  
  
}


mask_files <- list.files(masked_directory)

n_sm <- length(mask_files)

# mask frequency thresholds
mask_frq_thresh <- c(0.01, 0.05)

joined_masks <- list()
# Extract masked regions and retain type of mask
for(sm in 1:length(mask_files)){
  temp_sm_mask <- data.table::fread(glue::glue("{masked_directory}{mask_files[sm]}")) %>%
    dplyr::select(CHROM, START_BIN, END_BIN, STRAIN, window_mask) %>%
    dplyr::filter(window_mask != "Pass")
  
  
  joined_masks[[sm]] <- join_masks(temp_sm_mask)
  
  if(!exists("population_masks")){
    population_masks <- temp_sm_mask
  } else {
    population_masks <- dplyr::bind_rows(population_masks, temp_sm_mask)
  }
}

joined_masks <- dplyr::bind_rows(joined_masks)

# get mask frequency
all_masks <- population_masks %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(frq = n()/n_sm) %>%
  dplyr::mutate(frq = ifelse(frq > 0.5, 1-frq, frq))

population_masks %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(frq = n()/n_sm) %>%
  ggplot() +
    aes(x = END_BIN/1e6, y = frq)+
    geom_point()+
    facet_grid(.~CHROM, scales = "free", space = "free") +
  theme_bw(18) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))+
  labs(x = "Genomic Position (Mb)", y = "Frequency")

ggsave("Plots/Diversity/Masked_region_frequency_genome_wide.pdf", height = 6, width = 18)

population_masks %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(frq = n()/n_sm) %>%
  ggplot() +
  aes(x = frq)+
  geom_histogram(binwidth = 0.01)+
  theme_bw(18) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))+
  labs(x = "Frequency", y = "Count")

ggsave("Plots/Diversity/Masked_region_frequency_histogram.pdf", height = 6, width = 10)



population_masks %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(frq = n()/n_sm) %>%
  dplyr::filter(frq < 0.004) %>% nrow()
# [1] 8045
population_masks %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(frq = n()/n_sm) %>%
  dplyr::filter(frq < 0.01) %>% nrow()
# [1] 13760
population_masks %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(frq = n()/n_sm) %>%
  dplyr::filter(frq < 0.05) %>% nrow()
# [1] 24816
population_masks %>%
  dplyr::group_by(CHROM, START_BIN, END_BIN) %>%
  dplyr::summarise(frq = n()/n_sm) %>%
  nrow()
# [1] 47014

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
    dplyr::select(WBGeneID = V8, TranscriptID=V9) %>%
    dplyr::mutate(FreqCutoff = cutoff)
  
  if(!exists("masked_genes")){
    masked_genes <- temp_genes
  } else {
    masked_genes <- dplyr::bind_rows(masked_genes, temp_genes)
  }
}

# perform enrichment analysis
enrich_go <- function(x = NULL, freq_cut){
  
  # GO term enrichment  
  wb_ids <- dplyr::filter(x, FreqCutoff == freq_cut) %>% dplyr::pull(WBGeneID)
  
  mf <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  bp <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  # KEGG enrichment
  trancript_ids <- dplyr::filter(x, FreqCutoff == freq_cut) %>% dplyr::pull(TranscriptID)
  
  kk <- enrichKEGG(gene         = trancript_ids,
                   organism     = 'cel',
                   pvalueCutoff = 0.05)
  
  entrez_names <- getBM(attributes = c('wormbase_gene', 'entrezgene', 'uniprot_gn'), mart = ensembl) %>%
    dplyr::distinct(wormbase_gene, .keep_all = T)
  
  entrez_search_names <- dplyr::filter(entrez_names, wormbase_gene %in% c(dplyr::filter(x, FreqCutoff == freq_cut) %>% dplyr::pull(WBGeneID))) %>% 
    dplyr::distinct(uniprot_gn) %>%
    dplyr::pull(uniprot_gn)
  
  entr_enrich <- rev(c(1:length(entrez_search_names)))
  names(entr_enrich) <- entrez_search_names
  
  kk2 <- gseKEGG(geneList     = entr_enrich,
                 organism     = 'cel',
                 nPerm        = 1000,
                 keyType = "uniprot",
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  
  
  return(list(mf, bp, kk, kk2))
}

dir.create("Processed_Data/Enrichment")
dir.create("Plots/Enrichment")

# names of enrichment function output, used for saving below
etyps <- c(1:4)
names(etyps) <- c("GO_Molecular_Function", "GO_Biological_Process","Kegg_Enrichment", "Kegg_GeneSet")

# set up data base to extract gene ids
ensembl = useMart("ensembl",dataset="celegans_gene_ensembl")

for(frequency_cutoff in c("Low","Intermediate","Common")) {
  
  enrichment_results <- enrich_go(x = masked_genes, freq_cut = frequency_cutoff)
  
  save(enrichment_results, file = glue::glue("Processed_Data/Enrichment/{frequency_cutoff}_enrichment_results.rda"))
  
  for(enrich_type in 1:4){
    
    dp <- tryCatch({
      dotplot(enrichment_results[[enrich_type]], showCategory = 20)
    }, error = function(error_condition) {
      x <- c("NO ENRICHMENT")
      return(x)
    })
    
    if(!is.character(dp)){
      ggsave(dp, file = glue::glue("Plots/Enrichment/{frequency_cutoff}_{names(etyps)[enrich_type]}_Dotplot.pdf"),
             height = 6, width = 18)
    }
    
    cnp <- tryCatch({
      clusterProfiler::cnetplot(enrichment_results[[enrich_type]])
    }, error = function(error_condition) {
      x <- c("NO ENRICHMENT")
      return(x)
    }) 
    
    if(!is.character(cnp)){
      ggsave(cnp, file = glue::glue("Plots/Enrichment/{frequency_cutoff}_{names(etyps)[enrich_type]}_Cnet.pdf"),
             height = 12, width = 12)
    }
    
    hp <- tryCatch({
      heatplot(enrichment_results[[enrich_type]])
    }, error = function(error_condition) {
      x <- c("NO ENRICHMENT")
      return(x)
    }) 
    
    if(!is.character(hp)){
      tryCatch({
        ggsave(hp, file = glue::glue("Plots/Enrichment/{frequency_cutoff}_{names(etyps)[enrich_type]}_Heatmap.pdf"),
               height = 6, width = 18)
      }, error = function(error_condition) {
        print("NO ENRICHMENT")
      }) 
    }
    
    emap <- tryCatch({
      emapplot(enrichment_results[[enrich_type]])
    }, error = function(error_condition) {
      x <- c("NO ENRICHMENT")
      return(x)
    })
    
    if(!is.character(emap)){
      ggsave(emap, file = glue::glue("Plots/Enrichment/{frequency_cutoff}_{names(etyps)[enrich_type]}_EnrichmentMap.pdf"),
             height = 10, width = 12)
    }
    
    if(enrich_type<3){
      gp <- tryCatch({
        goplot(enrichment_results[[enrich_type]])
      }, error = function(error_condition) {
        x <- c("NO ENRICHMENT")
        return(x)
      })
      
      if(!is.character(gp)){
        ggsave(gp, file = glue::glue("Plots/Enrichment/{frequency_cutoff}_{names(etyps)[enrich_type]}_GOplot.pdf"),
               height = 16, width = 24)
      }
    }
  }
}




