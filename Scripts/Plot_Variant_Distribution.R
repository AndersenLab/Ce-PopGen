#!/usr/bin/env Rscript

library(tidyverse)
library(maps)
library(ggthemes)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# run script in base Ce-PopGen directory

# args
# 1 - small variant counts bed file
# 2 - small variant effect file
# 3 - structural variant bed file

# example - 
# args <- c("INTERSECTION_variant_snp_counts.txt","INTERSECTION_variant_indel_counts.txt", "INTERSECTION_snpeff.bed.gz","combinedSV.bed")
# args <- c("GATK-STRELKA_Intersection")
# args <- c("DIVERGENT-MASKED_Complete")
# args <- c("GATK-STRELKA_Intersection_Complete")
args <- commandArgs(TRUE)

data_dir <- "Data/Diversity/"

dir.create("Plots/Diversity")


snp_v <- data.table::fread(glue::glue("{data_dir}{args[1]}"), col.names = c("CHROM","w_start","w_end","count"))
indel_v <- data.table::fread(glue::glue("{data_dir}{args[2]}"), col.names = c("CHROM","w_start","w_end","count"))
large_v <- data.table::fread(glue::glue("{data_dir}{args[4]}"))
small_v_eff <- readr::read_tsv(glue::glue("{data_dir}{args[3]}"), col_names = F)

sample_names <- data.table::fread(glue::glue("data/GATK-STRELKA_Intersection_Complete_samples.txt"), header = F) %>% dplyr::pull(V1)

# 
# > sum(indel_v$count)
# [1] 390833
# > sum(snp_v$count)
# [1] 2369213

# small variant count per bin
snp_v %>%
  dplyr::mutate(norm_count = count/max(count)) %>%
  ggplot()+
  aes(x = w_start/1e6, y = norm_count)+
  # geom_point(alpha = 0.5)+
  geom_smooth()+
  geom_smooth(color = "red", data = indel_v %>%
                dplyr::mutate(norm_count = count/max(count)))+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  theme_bw(18) +
  labs(x = "Genomic Position (Mb)", y = "Normalized Variant Count") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("Plots/Diversity/Small_Variant_Count.pdf", height = 6, width = 12)

snp_v %>%
  ggplot()+
  aes(x = w_start/1e6, y = count)+
  geom_point(alpha = 0.5)+
  facet_grid(.~CHROM, space = "free", scales = "free")+
  theme_bw(18) +
  labs(x = "Genomic Position (Mb)", y = "SNV Count") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("Plots/Diversity/SNV_Variant_Count_point.pdf", height = 6, width = 12)

test <- small_v_eff %>%
  transform(X6 = strsplit(X6, ";")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(freq_miss = sum(grepl("\\./\\.",X6)/330),
                freq_alt = sum(grepl("1/1",X6)/330),
                freq_ref = sum(grepl("0/0",X6))/330,
                miss_st = gsub("=\\./\\.=|=\\./\\.",",",paste(grep("\\./\\.",X6, value = T), collapse = "=")),
                alt_st = gsub("=1/1=|=1/1",",",paste(grep("1/1",X6, value = T), collapse = "=")),
                ref_st = gsub("=0/0=|=0/0",",",paste(grep("0/0",X6, value = T), collapse = "="))) %>%
  dplyr::select(CHROM = X1, POS = X2, REF = X4, ALT = X5, X7:ref_st) %>%
  dplyr::filter(freq_alt != 0 )

write_tsv(test, glue::glue("{data_dir}Processed_Variant_Effects.tsv.gz"),col_names = F)

for(sm in sample_names){
  alt_ct <- length(grep(glue::glue("{sm},"), test$alt_st))
  temp_alt_ct <- data.frame(sample = sm, alt_ct = alt_ct)
  if(!exists("sm_alt_ct")){
    sm_alt_ct <- temp_alt_ct
  } else {
    sm_alt_ct <- dplyr::bind_rows(sm_alt_ct, temp_alt_ct)
  }
}

isolation_info <- googlesheets::gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
  googlesheets::gs_read("WI C. elegans") %>%
  dplyr::filter(reference_strain == 1)%>%
  dplyr::select(sample = isotype, long = longitude, lat = latitude, state, country)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::distinct(sample, long, lat, .keep_all = TRUE) %>%
  dplyr::left_join(sm_alt_ct, ., by = "sample") %>%
  dplyr::arrange((alt_ct)) %>%
  dplyr::mutate(size_point = scale(1:n())) %>%
  dplyr::mutate(size_point = ifelse(size_point <= quantile(size_point, probs = 0.90), 1, 2),
                lat = as.numeric(lat),
                long = as.numeric(long))

world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

write_tsv(isolation_info, glue::glue("{data_dir}SM_ALT_counts.tsv"),col_names = F)


ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="black", fill="white", size=0.5)+
  theme_map(18)+
  geom_point(aes(long, lat, fill = alt_ct, size = factor(size_point)),data=isolation_info, color = "black", alpha = 0.5, shape = 21)+
  scale_fill_viridis_c(option = "A")+
  labs(fill = "ALT Count", size = "90th Quantile")

ggsave("Plots/Diversity/ALT_Count_World_Map.pdf", height = 10, width = 20)

ggplot(isolation_info)+
  aes(x = (alt_ct)) +
  geom_histogram()+
  theme_bw(18) +
  labs(x = "ALT Count", y = "Number of strains")

ggsave("Plots/Diversity/ALT_Count_histogram.pdf", height = 6, width = 8)

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

# test %>%
#   dplyr::group_by(CHROM) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(X7 = ifelse(is.na(X7), "INTERGENIC", X7)) %>%
#   ggplot()+
#   aes(x = freq_alt,fill = X7)+
#   geom_density(binwidth = 0.01, alpha = 0.5)+
#   theme_bw(18) +
#   # facet_grid(X7~.)+
#   theme(axis.text.x = element_text(size =8))+
#   labs(x = "SnpEff Prediction", y = "Count")


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
#[1] 8776

# filter svs found by two variants
multi_support_sv <- large_v %>%
  dplyr::filter(SVTYPE_CLEAN!="SVTYPE_CLEAN") %>%
  dplyr::mutate(SIZE = as.numeric(SIZE),
                SUPPORT = as.numeric(SUPPORT)) %>%
  dplyr::filter(SIZE < 1e5, CHROM != "MtDNA") %>%
  dplyr::filter(SUPPORT > 2 | SVTYPE == "INS") 

# total svs
multi_support_sv %>%
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  dplyr::group_by(SVTYPE_CLEAN) %>%
  dplyr::summarise(sv_ct = n(),
                   med_size = median(as.numeric(SIZE), na.rm = T))%>%
  dplyr::ungroup()%>%
  dplyr::summarise(ssv = sum(sv_ct))

# genomic distribution
multi_support_sv %>%
  dplyr::distinct(CHROM, START, END,SVTYPE_CLEAN, .keep_all = T) %>%
  ggplot() +
  aes(x = as.numeric(START)/1e6, fill = SVTYPE_CLEAN) +
  geom_histogram(binwidth = 0.15)+
  facet_grid(.~CHROM, scales = "free")+
  theme_bw(18) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))+
  scale_fill_manual(values =c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(x = "Genomic Position (Mb)", fill = "SVTYPE", y = "Count") 

ggsave("Plots/Diversity/Structural_Variant_Distribution_by_type.pdf", height = 6, width = 18)

# counts and size by type
multi_support_sv %>%
dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  dplyr::mutate(SIZE = as.numeric(SIZE)) %>%
  dplyr::group_by(SVTYPE_CLEAN) %>%
  dplyr::summarise(sv_ct = n(),
                   med_size = median(as.numeric(SIZE), na.rm = T))

multi_support_sv %>%
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  dplyr::mutate(SIZE = as.numeric(SIZE)) %>%
  dplyr::filter(SVTYPE_CLEAN!="INS")%>%
  ggplot()+
  aes(x = SVTYPE_CLEAN, y = log(SIZE), fill = SVTYPE_CLEAN)+
  geom_boxplot() +
  scale_fill_manual(values =c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  theme_bw(18)+
  theme(axis.title.x = element_blank())+
  labs(fill = "SVTYPE")
ggsave("Plots/Diversity/Structural_size_boxplot.pdf", height = 4, width = 8)

# counts and size by type and chrom
multi_support_sv %>%
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  dplyr::mutate(SIZE = as.numeric(SIZE)) %>%
  dplyr::group_by(SVTYPE_CLEAN, CHROM) %>%
  dplyr::summarise(sv_ct = n(),
                   med_size = median(as.numeric(SIZE), na.rm = T)) %>%
  dplyr::mutate(frac_by_chrom = sv_ct/sum(sv_ct)) %>%
  dplyr::group_by(CHROM) %>%
  dplyr::mutate(avg_by_chrom = mean(frac_by_chrom)) %>%
  dplyr::arrange(SVTYPE_CLEAN,desc(frac_by_chrom))%>%
  View()

# frequency
multi_support_sv %>%
  dplyr::distinct(CHROM, START, END, SAMPLE, .keep_all = T) %>%
  dplyr::group_by(CHROM, START, END, SVTYPE_CLEAN) %>%
  dplyr::mutate(st_ct = n()/330) %>%
  dplyr::ungroup()%>% 
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  # dplyr::group_by(SVTYPE_CLEAN) %>%
  dplyr::summarise(med_f = median(st_ct),
                   meanf = mean(st_ct))

multi_support_sv %>%
  dplyr::filter(SVTYPE_CLEAN!="SVTYPE_CLEAN") %>%
  dplyr::distinct(CHROM, START, END, SAMPLE, .keep_all = T) %>%
  dplyr::group_by(CHROM, START, END, SVTYPE_CLEAN) %>%
  dplyr::mutate(st_ct = n()/330) %>%
  dplyr::distinct(st_ct,.keep_all=T) %>%
  ggplot()+
  aes(x = st_ct)+
  geom_histogram(binwidth = 0.01)+
  theme_bw(18)+
  labs(x = "Structural Variant Frequency", y = "Count")

ggsave("Plots/Diversity/Structural_Variant_Frequency.pdf", height = 4, width = 8)

# frequency by type
multi_support_sv %>%
  dplyr::distinct(CHROM, START, END, SAMPLE, .keep_all = T) %>%
  dplyr::group_by(CHROM, START, END, SVTYPE_CLEAN) %>%
  dplyr::mutate(st_ct = n()/330) %>%
  dplyr::ungroup()%>% 
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  dplyr::group_by(SVTYPE_CLEAN) %>%
  dplyr::summarise(med_f = median(st_ct),
                   meanf = mean(st_ct))

# sv by genomic region
multi_support_sv %>%
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  dplyr::mutate(clean_region = ifelse(SNPEFF_PRED == "intron_variant", "intron", 
                                      ifelse(TRANSCRIPT != "Intergenic", "genic","Intergenic")))%>%
  dplyr::ungroup()%>% 
  dplyr::group_by(clean_region) %>%
  dplyr::summarise(region_ct = n())

# sv by genomic region - genic onlu
multi_support_sv %>%
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  dplyr::mutate(clean_region = ifelse(SNPEFF_PRED == "intron_variant", "intron", 
                                      ifelse(TRANSCRIPT != "Intergenic", "genic","Intergenic")))%>%
  dplyr::filter(clean_region == "genic")%>%
  dplyr::ungroup()%>% 
  dplyr::group_by(SNPEFF_EFF) %>%
  dplyr::summarise(ef_ct = n())

# sv by type by genomic region - genic onlu
multi_support_sv %>%
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  dplyr::mutate(clean_region = ifelse(SNPEFF_PRED == "intron_variant", "intron", 
                                      ifelse(TRANSCRIPT != "Intergenic", "genic","Intergenic")))%>%
  dplyr::filter(clean_region == "genic")%>%
  dplyr::ungroup()%>% 
  dplyr::group_by(SNPEFF_EFF,SVTYPE_CLEAN) %>%
  dplyr::summarise(ef_ct = n())

# sv by type by genomic region - genic onlu - high only
multi_support_sv %>%
  dplyr::distinct(CHROM, START, END, .keep_all = T) %>%
  dplyr::mutate(clean_region = ifelse(SNPEFF_PRED == "intron_variant", "intron", 
                                      ifelse(TRANSCRIPT != "Intergenic", "genic","Intergenic")))%>%
  dplyr::filter(clean_region == "genic",SNPEFF_EFF=="HIGH")%>%
  dplyr::ungroup()%>%
  dplyr::mutate(u_genes = length(unique(TRANSCRIPT)))


large_v %>%
  dplyr::mutate(SIZE = as.numeric(SIZE)) %>%
  dplyr::filter(SIZE < 1e5, CHROM != "MtDNA") %>%
  dplyr::mutate(SUPPORT = as.numeric(SUPPORT)) %>%
  dplyr::filter(SUPPORT > 2 || SVTYPE == "INS") %>%
  dplyr::distinct(CHROM, START, END,SVTYPE_CLEAN, .keep_all = T) %>%
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

ggsave("Plots/Diversity/Structural_Variant_Distribution_by_type_size.pdf", height = 10, width = 18)

