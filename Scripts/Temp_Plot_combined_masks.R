library(maps)
library(ggthemes)
library(ggrepel)

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
  dplyr::filter(common_cluster_size>5000) %>%
  dplyr::distinct(CHROM, common_cluster_start, common_cluster_end, STRAIN, .keep_all = T) %>%
  dplyr::group_by(CHROM, common_cluster_start, common_cluster_end) %>%
  dplyr::mutate(ct = n())%>%
  ggplot()+
  geom_segment(aes(x = common_cluster_start/1e6, y = ct, xend = common_cluster_end/1e6, yend = ct, colour = common_cluster_size))+
  facet_grid(CHROM~., space="free")+
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
