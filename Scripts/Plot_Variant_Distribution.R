#!/usr/bin/env Rscript

library(tidyverse)


# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# run script in base Ce-PopGen directory

# args
# 1 - small variant counts bed file
# 2 - small variant effect file
# 3 - structural variant bed file

# example - 
# args <- c("WI_GATK_variant_counts.txt", "WI.HARD-FILTERED.bed.gz","combinedSV.bed")
# args <- c("GATK-STRELKA_Intersection")
# args <- c("DIVERGENT-MASKED_Complete")
# args <- c("GATK-STRELKA_Intersection_Complete")
args <- commandArgs(TRUE)

data_dir <- "Data/Diversity/"

dir.create("Plots/Diversity")


small_v <- data.table::fread(glue::glue("{data_dir}{args[1]}"), col.names = c("CHROM","w_start","w_end","count"))
large_v <- data.table::fread(glue::glue("{data_dir}{args[3]}"))
small_v_eff <- readr::read_tsv(glue::glue("{data_dir}{args[2]}"), col_names = F)


# small variant count per bin
ggplot(small_v)+
  aes(x = w_start/1e6, y = count)+
  geom_point(alpha = 0.5)+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  theme_bw(18) +
  labs(x = "Genomic Position (Mb)", y = "Variant Count") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("Plots/Diversity/Small_Variant_Count.pdf", height = 6, width = 12)

test <- small_v_eff %>%
  transform(X4 = strsplit(X4, ";")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(freq_miss = sum(grepl("\\./\\.",X4)/330),
                freq_alt = sum(grepl("1/1",X4)/330),
                freq_ref = sum(grepl("0/0",X4))/330) %>%
  dplyr::select(CHROM = X1, POS = X2, X5:freq_ref)

write_tsv(test, glue::glue("{data_dir}Processed_Variant_Effects.tsv.gz"),col_names = F)

hi_eff <- test %>%
  dplyr::group_by(CHROM) %>%
  dplyr::ungroup() %>%
  dplyr::filter(X7 == "HIGH", X6 != "bidirectional_gene_fusion") %>%
  dplyr::distinct(X9, X6, .keep_all = T)

# small variant count per bin
hi_eff %>%
  dplyr::group_by(X6) %>%
  dplyr::summarise(ct = n()) %>%
  dplyr::arrange(desc(ct)) %>%
  dplyr::filter(ct > 200) %>%
  dplyr::mutate(X6 = gsub("&","\n",X6)) %>%
  dplyr::mutate(eff_type = factor(X6,levels = X6))%>%
  ggplot()+
  aes(x = eff_type, y = ct)+
  geom_bar(stat = "identity") +
  theme_bw(18) +
  theme(axis.text.x = element_text(size =8))+
  labs(x = "SnpEff Prediction", y = "Count")

ggsave("Plots/Diversity/Small_Variant_SnpEff_High.pdf", height = 6, width = 18)

# Distribution of high impact variation
hi_eff %>%
  dplyr::group_by(X6) %>%
  dplyr::mutate(ct = n()) %>%
  dplyr::arrange(desc(ct)) %>%
  dplyr::filter(ct > 200, CHROM!="MtDNA") %>%
  dplyr::ungroup()%>%
  dplyr::distinct(CHROM, POS, .keep_all = T) %>%
  dplyr::mutate(X6 = gsub("&","\n",X6)) %>%
  dplyr::mutate(eff_type = factor(X6,levels = unique(X6)))%>%
  ggplot()+
  aes(x = POS/1e6, y = freq_alt)+
  geom_point(alpha = 0.5)+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  theme_bw(18) +
  labs(x = "Genomic Position (Mb)", y = "Count")+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("Plots/Diversity/Small_Variant_SnpEff_High_Distribution.pdf", height = 6, width = 18)

# unique genes with predicted LOF
length(unique(hi_eff$X9))
#[1] 10484

large_v %>%
  dplyr::distinct(CHROM, START, END, SAMPLE, .keep_all = T) %>%
  dplyr::group_by(CHROM, START, END, SVTYPE_CLEAN) %>%
  dplyr::mutate(st_ct = n()/330) %>%
  dplyr::distinct(st_ct) %>%
  ggplot()+
  aes(x = st_ct)+
  geom_histogram(binwidth = 0.01)+
  theme_bw(18)+
  labs(x = "Structural Variant Frequency", y = "Count")

ggsave("Plots/Diversity/Structural_Variant_Frequency.pdf", height = 4, width = 8)

large_v %>%
  dplyr::mutate(SIZE = as.numeric(SIZE)) %>%
  dplyr::filter(SIZE < 1e5, CHROM != "MtDNA") %>%
  dplyr::mutate(SUPPORT = as.numeric(SUPPORT)) %>%
  dplyr::filter(SUPPORT > 2 || SVTYPE == "INS") %>%
  ggplot() +
  aes(x = as.numeric(START)/1e6, xend = as.numeric(END)/1e6, y=SIZE, yend=SIZE, color = SVTYPE_CLEAN) +
  geom_segment(size =1, arrow = arrow(length = unit(0.1,"cm")))+
  facet_grid(CHROM~HIGH_EFF, scales = "free")+
  theme_bw(18) +
  theme(strip.background = element_blank(),
        legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray70"),
        panel.grid = element_blank())+
  scale_color_manual(values =c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(x = "Genomic Position (Mb)", color = "SVTYPE", y = "Structural Variant Size") 

ggsave("Plots/Diversity/Structural_Variant_Distribution.pdf", height = 10, width = 18)

