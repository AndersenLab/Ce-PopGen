library(maps)
library(ggthemes)
library(ggrepel)

k = 14 
admix_pops <- data.table::fread("Processed_Data/Admixture/GATK-STRELKA_Intersection_Complete_Ancestral_Pops.tsv") %>%
  dplyr::filter(K== k)%>%
  dplyr::rename(strain=samples) %>%
  dplyr::filter(frac_cluster == max_frac) %>%
  dplyr::distinct()


isolation_info <- googlesheets::gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
  googlesheets::gs_read("WI C. elegans") %>%
  dplyr::filter(reference_strain == 1)%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude, state, country)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(admix_pops,.,by="strain")%>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::distinct(strain, long, lat, .keep_all = TRUE)

isolation_info$lat <- as.numeric(isolation_info$lat)
isolation_info$long <- as.numeric(isolation_info$long)


world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

ggplot()+ geom_map(data=world, map=world,
                                 aes(x=long, y=lat, map_id=region),
                                 color="black", fill="#7f7f7f", size=0.5)+
  scale_fill_manual(values = ancestry.colours,name = "Population")+
  theme_map()+
  geom_label_repel(aes(long, lat, label = strain, fill = cluster),
                   data=dplyr::filter(isolation_info, state == "Hawaii"),
                   fontface = 'bold', color = 'white',
                   segment.color = 'cyan')+
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-160,-155), ylim = c(19,22.3))


ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="black", fill="#7f7f7f", size=0.5)+
  scale_fill_manual(values = ancestry.colours,name = "Population")+
  theme_map()+
  geom_label_repel(aes(long, lat, label = strain, fill = cluster),
                   data=dplyr::filter(isolation_info, max_frac > 0.9),
                   fontface = 'bold', color = 'white',
                   segment.color = 'cyan')+
  theme(legend.position = "none")

ggsave(world_plot, filename = glue::glue("plots/K{K}_HW_MAP_Admix.pdf"),
       height = 10,
       width = 10)