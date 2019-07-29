library(maps)
library(ggthemes)
library(ggrepel)
# devtools::install_github('kevinblighe/PCAtools')
library(PCAtools)
library(tidyverse)
library(GGally)
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))
pr_mask_directory <- "Data/Divergent_Masks/Masks_By_Type/"
joined_masks <- read_tsv( glue::glue("{pr_mask_directory}Combine_neighboring_masks.tsv.gz"), col_names = T)

source("Scripts/Figure_Themes.R")

joined_masks %>%
  dplyr::filter(common_cluster_size>1000) %>%
  ggplot()+
  aes(x = common_cluster_size) +
  geom_histogram(binwidth = 1000)

joined_masks %>%
  dplyr::filter(common_cluster_size>1000) %>%
  dplyr::distinct(CHROM, common_cluster_start, common_cluster_end, STRAIN, .keep_all = T) %>%
  dplyr::group_by(CHROM, common_cluster_start, common_cluster_end) %>%
  dplyr::mutate(ct = n())%>%
  dplyr::ungroup() %>%
  ggplot()+
  aes(x = common_cluster_size, y=ct) +
  geom_point()

joined_masks %>%
  dplyr::filter(common_cluster_size>1000) %>%
  dplyr::distinct(CHROM, common_cluster_start, common_cluster_end, STRAIN, .keep_all = T) %>%
  dplyr::group_by(CHROM, common_cluster_start, common_cluster_end) %>%
  dplyr::mutate(ct = n()/330)%>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN, CHROM) %>%
  dplyr::mutate(total_divergent = sum(common_cluster_size)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(CHROM, STRAIN, total_divergent) %>%
  dplyr::mutate(total_divergent = (total_divergent)/1e6) %>%
  tidyr::spread(CHROM, total_divergent) %>%
  dplyr::mutate(n_strain = 1:n())%>%
  ggpairs(., columns=2:7,  
          upper = list(
            continuous = wrap('cor', method = "spearman")
          )) + 
  ggtitle("Total Divergent Regions") +
  base_theme 


joined_masks %>%
  dplyr::filter(common_cluster_size>1000) %>%
  dplyr::distinct(CHROM, common_cluster_start, common_cluster_end, STRAIN, .keep_all = T) %>%
  dplyr::group_by(CHROM, common_cluster_start, common_cluster_end) %>%
  dplyr::mutate(ct = n()/330)%>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN, CHROM) %>%
  dplyr::mutate(total_divergent = sum(common_cluster_size)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(total_divergent)) %>%
  dplyr::mutate(f_strain = factor(STRAIN, levels = unique(STRAIN)))%>%
  ggplot()+
  geom_segment(aes(x = common_cluster_start/1e6, y = f_strain, xend = common_cluster_end/1e6, yend = f_strain, colour = ct))+
  facet_grid(.~CHROM, space="free", scales = "free") +
  theme(axis.text.y = element_blank())+
  scale_color_viridis_c(direction = 1, option = "A")


joined_masks %>%
  dplyr::filter(common_cluster_size>5000) %>%
  dplyr::distinct(CHROM, common_cluster_start, common_cluster_end, STRAIN, .keep_all = T) %>%
  dplyr::group_by(CHROM, common_cluster_start, common_cluster_end) %>%
  dplyr::mutate(ct = n())%>%
  ggplot()+
  geom_segment(aes(x = common_cluster_start/1e6, y = ct, xend = common_cluster_end/1e6, yend = ct, colour = common_cluster_size))+
  facet_grid(.~CHROM, space="free")+
  scale_color_viridis_c(direction = 1, option = "A")

check_substrate_enrichment <- joined_masks %>%
  dplyr::filter(common_cluster_size>5000) %>%
  dplyr::distinct(CHROM, common_cluster_start, common_cluster_end, STRAIN, .keep_all = T) %>%
  dplyr::group_by(CHROM, common_cluster_start, common_cluster_end) %>%
  dplyr::mutate(ct = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(ct > 5 & common_cluster_size > 50000)

write.table(check_substrate_enrichment, file = "~/transfer/daehan_enrichment.tsv", sep = "\t", quote = F, col.names = T,row.names = F)

world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

loc_info <- googlesheets::gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
  googlesheets::gs_read("WI C. elegans") %>%
  dplyr::filter(reference_strain == 1)%>%
  dplyr::select(STRAIN = isotype, long = longitude, lat = latitude, state, country, substrate, landscape)%>%
  dplyr::filter(lat != "None") %>%
  dplyr::mutate(lat = as.numeric(lat),
                long = as.numeric(long)) %>%
  dplyr::pull(STRAIN)


check_substrate_enrichment <- joined_masks %>%
  dplyr::filter(STRAIN %in% loc_info) %>%
  dplyr::filter(common_cluster_size>5000) %>%
  dplyr::distinct(CHROM, common_cluster_start, common_cluster_end, STRAIN, .keep_all = T) %>%
  dplyr::group_by(CHROM, common_cluster_start, common_cluster_end) %>%
  dplyr::mutate(ct = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(ct > 5 & common_cluster_size > 50000)

check_substrate_enrichment$group_id <- check_substrate_enrichment %>% group_indices(CHROM, common_cluster_start, common_cluster_end) 

to_merge <- check_substrate_enrichment %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(common_cluster_size))


isolation_info <- googlesheets::gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
  googlesheets::gs_read("WI C. elegans") %>%
  dplyr::filter(reference_strain == 1)%>%
  dplyr::select(STRAIN = isotype, long = longitude, lat = latitude, state, country, substrate, landscape)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(.,to_merge,by="STRAIN")%>%
  dplyr::distinct(STRAIN, long, lat,group_id, .keep_all = TRUE) %>%
  dplyr::mutate(lat = as.numeric(lat),
                long = as.numeric(long)) %>%
  dplyr::arrange(desc(common_cluster_size)) %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(n_sub = length(unique(substrate)),
                n_land = length(unique(landscape)),
                n_country = length(unique(country)),
                n_state = length(unique(state)))


ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="black", fill="#7f7f7f", size=0.5)+
  theme_map()+
  geom_label_repel(aes(long, lat, label = STRAIN, fill = substrate),
                   data=dplyr::filter(isolation_info, group_id == 28),
                   segment.color = 'cyan')

# in close proximity
131
138
# hawaii
143
98
1
# largest cluster
28


comb_mask_matrix <- joined_masks %>%
  tidyr::unite(mask_region, CHROM, common_cluster_start, common_cluster_end, sep="_") %>%
  dplyr::select(mask_region, STRAIN) %>%
  dplyr::mutate(GT = 1) %>%
  dplyr::distinct(mask_region, STRAIN, .keep_all=T) %>%
  tidyr::spread(STRAIN, GT)


comb_mask_matrix <- population_masks %>%
  tidyr::unite(mask_region, CHROM, START_BIN, END_BIN, sep="_") %>%
  dplyr::select(mask_region, STRAIN) %>%
  dplyr::mutate(GT = 1) %>%
  dplyr::distinct(mask_region, STRAIN, .keep_all=T) %>%
  tidyr::spread(STRAIN, GT)

comb_mask_matrix[is.na(comb_mask_matrix)] <- 0

isolation_info <- googlesheets::gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
  googlesheets::gs_read("WI C. elegans") %>%
  dplyr::filter(reference_strain == 1)%>%
  dplyr::select(STRAIN = isotype, long = longitude, lat = latitude, state, country, substrate, landscape)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::mutate(lat = as.numeric(lat),
                long = as.numeric(long))


metadata <- data.frame(isolation_info[,1:7]) %>%
  dplyr::arrange(STRAIN) %>%
  na.omit() %>%
  dplyr::mutate(n_sub = as.numeric(factor(substrate)),
                n_land = as.numeric(factor(landscape))) %>%
  dplyr::group_by(substrate) %>%
  dplyr::mutate(ct_sub = n()) %>%
  dplyr::filter(ct_sub > 5) %>%
  dplyr::group_by(landscape) %>%
  dplyr::mutate(ct_sub = n()) %>%
  dplyr::filter(ct_sub > 5) %>%
  dplyr::select(-ct_sub) %>%
  dplyr::ungroup()

row.names(metadata) <- metadata$STRAIN

x <- comb_mask_matrix[colnames(comb_mask_matrix) %in% metadata$STRAIN]

row.names(x)  <- comb_mask_matrix$mask_region

p <- pca(x, metadata = metadata)

screeplot(p,
          components = getComponents(p, 1:20),
          hline = 80, vline = 27) +
  geom_text(aes(20, 80, label = '80% explained variation', vjust = -1))

eigencorplot(p, metavars = c('long','lat',"n_sub","n_land"),components = getComponents(p, seq(20)))

pairsplot(p)
pairsplot(p,colby = "lat", getComponents(p, c(1,3,4,5,7,8,10)), legendPosition = "left")

plotloadings(p,components = getComponents(p, c(1,3,4,5,7,8,10)))


# substrate
ldplot <- plotloadings(p,components = getComponents(p, c(7,8,10,18,20)))

pairsplot(p,colby = "substrate", getComponents(p, c(7,8,10,18,20)), legendPosition = "left")

loadings_cor_subs <- ldplot$data



#
library(pegas)
require(adegenet)
library(ape)

x <- as.character(woodmouse)
x[, 1:20]
str(as.alignment(x))

#chrom II - eca701 and xz1516 look very different

comb_mask_matrix <- joined_masks %>%
  # dplyr::filter(common_cluster_size>5000) %>%
  dplyr::distinct(CHROM, common_cluster_start, common_cluster_end, STRAIN, .keep_all = T) %>%
  dplyr::group_by(CHROM, common_cluster_start, common_cluster_end) %>%
  dplyr::mutate(ct = n()/330)%>%
  dplyr::filter(ct > 0.05) %>%
  tidyr::unite(mask_region, CHROM, common_cluster_start, common_cluster_end, sep="_") %>%
  dplyr::select(mask_region, STRAIN) %>%
  dplyr::mutate(GT = "A") %>%
  dplyr::distinct(mask_region, STRAIN, .keep_all=T) %>%
  tidyr::spread(STRAIN, GT, fill = "T")

test_ape <- comb_mask_matrix[,2:ncol(comb_mask_matrix)]

row.names(test_ape) <- comb_mask_matrix$mask_region

al_ape <- ape::as.alignment(test_ape)

str(al_ape)

db_ape <- as.DNAbin(al_ape)

D <- dist.dna(db_ape)
tre <- nj(D)
plot.phylo(tre,align.tip.label = F, cex = .6)
rootedtr <- root(tre, "XZ1516")
plot.phylo(rootedtr,align.tip.label = F, cex = .6)


h <- pegas::haplotype(d)
h <- sort(h, what = "label")
(net <- pegas::haploNet(h))


h <- pegas::haplotype(db_ape)

hn <- haploNet(h)

plot(hn, size=attr(hn, "freq"), fast = F)

ind.hap<-with(stack(setNames(attr(h, "index"), rownames(h))), table(hap=ind, individuals=samples))


samples <- colnames(comb_mask_matrix[,2:ncol(comb_mask_matrix)])

plot(hn, size=attr(hn, "freq"), scale.ratio=0.2, pie=ind.hap)



checkAlignment(db_ape, check.gaps = TRUE, plot = TRUE, what = 1:4)



x <- woodmouse[sample(15, size = 110, replace = TRUE), ]
(h <- haplotype(x))

x <- as.character(woodmouse)
x[, 1:20]

