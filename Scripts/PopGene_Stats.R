library(tidyverse)


setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

analysis_type_a <- "GATK-STRELKA_Intersection_Complete"
analysis_type_b <- "DIVERGENT-MASKED_Complete"

load(glue::glue("Data/POPGENOME/{analysis_type_a}/WHOLE_GENOME/Ce_Genome-wide_Neutrality_stats.Rda"))
pg_a <- neutrality_df
load(glue::glue("Data/POPGENOME/{analysis_type_b}/WHOLE_GENOME/Ce_Genome-wide_Neutrality_stats.Rda"))
pg_b <- neutrality_df
rm(neutrality_df)

unique(pg_a$statistic)
p_st <- "Fu.Li.F"
# Pi
ggplot() +
  aes(x = WindowPosition/1e6, y = value)+
  facet_grid(CHROM~., scales = "free")+
  geom_point( color = "red", data = pg_a %>% dplyr::filter(statistic == p_st) )+
  geom_point(color = "blue", data = pg_b %>% dplyr::filter(statistic == p_st) )+
  geom_line( color = "pink", data = pg_a %>% dplyr::filter(statistic == p_st) )+
  geom_line(color = "cyan", data = pg_b %>% dplyr::filter(statistic == p_st) )+
  theme_bw()+
  theme(legend.box = NULL,
        strip.background = element_blank(),
        legend.background = element_rect(colour = NA),
        legend.key = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = p_st, color = "Population", fill = "Population")

