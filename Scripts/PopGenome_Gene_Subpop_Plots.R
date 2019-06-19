#!/usr/bin/env Rscript

library(pophelper)
library(tidyverse)
library(ggthemes)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# run script in base Ce-PopGen directory

# args
# 1 - analysis type - name of folder corresponding to what admixure analysis

# example - 
# args <- c("GATK-STRELKA_Intersection_Complete")
args <- commandArgs(TRUE)

analysis_type <- args[1]

data_dir <- glue::glue("Data/POPGENOME/{analysis_type}/SUBPOPs/GENE/")

# define a color pallette with maximal contrast
# Colors from - https://gist.github.com/ollieglass/f6ddd781eeae1d24e391265432297538
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "K"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")
# chnage names tp match FST df pop names
ancestry.colours2 <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                       "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                       "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "K"="black", 
                       "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                       "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")
names(ancestry.colours2) <- paste0("pop",c(1:length(ancestry.colours2)))

# extract chromosomes
chroms <- grep('[A-Z]', list.dirs(data_dir, full.names = F), value = T)


for(chrom in chroms){
  popgenome_files <- list.files(glue::glue("{data_dir}{chrom}"))
  
  for(k in popgenome_files) {
    
    ksize = as.numeric(gsub("K", "", as.character(strsplit(k, split = "_")[[1]][2])))
    
    load(glue::glue("{data_dir}{chrom}/{k}"))
    
    genes <- data.frame(Gene_region = GENOME_OBJECT@region.names) %>%
      tidyr::separate(Gene_region, into = c("Gene_Start", "Gene_End"), sep = " - ", convert = T)
    
    temp_fst <- data.frame(genes,
                           PopGenome::get.F_ST(GENOME_OBJECT, mode="nucleotide", pairwise = T)[[1]]) %>%
      tidyr::gather(Pops, Fst, -Gene_Start, -Gene_End) %>%
      dplyr::mutate(CHROM = chrom,
                    K = ksize) %>%
      dplyr::select(CHROM, Gene_Start, Gene_End, K, Pops, Fst) %>%
      na.omit()
    
    neutrality_ls <- list()
    for(popnumber in 1:ksize){
      popname <- LETTERS[[popnumber]]
      neutrality_ls[[popnumber]] <- data.frame(PopGenome::get.neutrality(GENOME_OBJECT, theta = T, stats = T)[[popnumber]]) %>%
        dplyr::mutate(Population = popname,
                      Gene_Start = genes$Gene_Start,
                      Gene_End = genes$Gene_End,
                      CHROM = chrom) 
    }
    
    temp_nd <- dplyr::bind_rows(neutrality_ls) %>%
      tidyr::gather(Statistic, Value, -(Population:CHROM)) %>%
      dplyr::select(CHROM, Gene_Start, Gene_End, Population, Statistic, Value) %>%
      dplyr::mutate(K = ksize) %>%
      na.omit()
    
    if(!exists("fst_df")) {
      fst_df <- temp_fst
      nd_df <- temp_nd
    } else {
      fst_df <- dplyr::bind_rows(fst_df, temp_fst)
      nd_df <- dplyr::bind_rows(nd_df, temp_nd)
    }
  }
}

# fix pop names to match other analyses
clean_fst <- fst_df%>%
  tidyr::separate(Pops, into = c("PopA" ,"PopB"), sep = "\\.") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(PopA = names(ancestry.colours[ancestry.colours==ancestry.colours2[PopA]]),
                PopB = names(ancestry.colours[ancestry.colours==ancestry.colours2[PopB]]))



dir.create(glue::glue("Processed_Data/Popgenome/{analysis_type}/GENE"),recursive = T)
write_tsv(nd_df, path = glue::glue("Processed_Data/Popgenome/{analysis_type}/GENE/SUBPOP_Gene_ND.tsv.gz"))

for(kpop in unique(clean_fst$K)){
  save_fst <- clean_fst %>% dplyr::filter(K == kpop)
  write_tsv(save_fst, path = glue::glue("Processed_Data/Popgenome/{analysis_type}/GENE/SUBPOP_Gene_Fst_K{kpop}.tsv.gz"))
}

kpop = 15
fst_df <- readr::read_tsv(glue::glue("Processed_Data/Popgenome/{analysis_type}/GENE/SUBPOP_Gene_Fst_K{kpop}.tsv.gz"))

fst_df %>%
  dplyr::filter(PopA == "H", PopB == "L") %>%
  ggplot()+
  aes(x = Gene_Start/1e6, y = Fst) +
  geom_point()+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw(18)+
  ylim(0,1)+
  labs(x = "Genomic Position (Mb)")


nd_df <- readr::read_tsv(glue::glue("Processed_Data/Popgenome/{analysis_type}/GENE/SUBPOP_Gene_ND.tsv.gz"))

plot_st <- "Tajima.D"
nd_df %>%
  dplyr::filter(Population %in% c("B","D","E","F","L","O"), K == kpop, Statistic == plot_st) %>%
  ggplot()+
  aes(x = Gene_Start/1e6, y = Value, color = Population) +
  geom_point()+
  scale_color_manual(values = ancestry.colours) +
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw(18)+
  labs(x = "Genomic Position (Mb)")
