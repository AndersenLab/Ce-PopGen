#!/usr/bin/env Rscript

library(tidyverse)
library(clusterProfiler)
library(org.Ce.eg.db)
library(biomaRt)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# run script in base Ce-PopGen directory

# args
# 1 - analysis type - name of folder corresponding to what admixure analysis

# example - 
# args <- c("GATK-STRELKA_Intersection_Complete")
args <- commandArgs(TRUE)

analysis_type <- args[1]

ensembl = useMart("ensembl",dataset="celegans_gene_ensembl")
data_dir <- glue::glue("Data/POPGENOME/{analysis_type}/WHOLE_GENOME/")
pr_mask_directory <- "Data/Divergent_Masks/Masks_By_Type/"

p_st <- "Tajima.D"
st_cut <- 2

load(glue::glue("{data_dir}Ce_Genome-wide_Neutrality_stats.Rda")) 


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
  
  entrez_names <- getBM(attributes = c('wormbase_gene', 'entrezgene_id', 'uniprotsptrembl'), mart = ensembl) %>%
    dplyr::distinct(wormbase_gene, .keep_all = T)
  
  entrez_search_names <- dplyr::filter(entrez_names, wormbase_gene %in% c(dplyr::filter(x, FreqCutoff == freq_cut) %>% dplyr::pull(WBGeneID))) %>% 
    dplyr::distinct(uniprotsptrembl) %>%
    dplyr::pull(uniprotsptrembl)
  
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



test_join <- neutrality_df %>%
  dplyr::filter(statistic == p_st, value > st_cut) %>% 
  dplyr::mutate(left_r = WindowPosition-10000,
                right_r = WindowPosition+10000) %>%
  dplyr::select(CHROM, left_r, right_r, value)



write.table(test_join, 
            file = glue::glue("{pr_mask_directory}Window_Balancing_selection.bed"), 
            col.names = F, row.names = F, quote = F, sep = "\t")

system("Scripts/Balancing_selection_popgenome_window.sh")

atype <- "GATK-STRELKA_Intersection_Complete_POPPGENOME_window"

gene_files <- data.table::fread(glue::glue("{pr_mask_directory}Window_Balancing_Selection_Genes.bed")) %>%
  dplyr::select(WBGeneID = V8, TranscriptID=V9) %>%
  dplyr::mutate(FreqCutoff = atype)

# names of enrichment function output, used for saving below
etyps <- c(1:4)
names(etyps) <- c("GO_Molecular_Function", "GO_Biological_Process","Kegg_Enrichment", "Kegg_GeneSet")


for(frequency_cutoff in c(atype)) {
  
  enrichment_results <- enrich_go(x = gene_files, freq_cut = frequency_cutoff)
  
  save(enrichment_results, file = glue::glue("Processed_Data/Enrichment/{frequency_cutoff}_enrichment_results.rda"))
  
  for(enrich_type in 1:4){
    
    dp <- tryCatch({
      dotplot(enrichment_results[[enrich_type]])
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



