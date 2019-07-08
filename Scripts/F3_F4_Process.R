#!/usr/bin/env Rscript

library(pophelper)
library(tidyverse)
library(ggthemes)
library(admixturegraph)

# setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# run script in base Ce-PopGen directory

# args
# 1 - analysis type - name of folder corresponding to what admixure analysis

# example - 
# args <- c("GATK-STRELKA_Intersection_Complete", "14")
args <- commandArgs(TRUE)

analysis_type <- args[1]
analysis_k <- args[2]

data_dir <- glue::glue("Data/ADMIXTURE/{analysis_type}/Treemix")

# define a color pallette with maximal contrast
# Colors from - https://gist.github.com/ollieglass/f6ddd781eeae1d24e391265432297538
ancestry.colours <- c("A"="gold2", "B"="plum4","C"= "darkorange1", 
                      "D"="lightskyblue2", "E"="firebrick","F"= "burlywood3", "G"="gray51", 
                      "H"="springgreen4", "I"="lightpink2", "J"="deepskyblue4", "K"="black", 
                      "L"="mediumpurple4","M"= "orange","N"= "maroon","O"= "yellow3","P"= "brown4", 
                      "Q"="yellow4", "R"="sienna4", "S"="chocolate", "T"="gray19")


f3_file <- grep(glue::glue("K-{analysis_k}_"), grep("threepop", list.files(glue::glue("{data_dir}")), value=T),value=T)
f4_file <- grep(glue::glue("K-{analysis_k}_"), grep("fourpop", list.files(glue::glue("{data_dir}")), value=T),value=T)

fthree <- data.table::fread(glue::glue("{data_dir}/{f3_file}") , sep = " ")%>%
  tidyr::separate(V1, into = c("Outgroup","otherpops"), sep = ";")%>%
  tidyr::separate(otherpops, into = c("PopA","PopB"), sep = ",")%>%
  dplyr::rename(F3 = V2, SE = V3, Z.value = V4)%>%
  dplyr::group_by(Outgroup)%>%
  dplyr::mutate(mf3 = median(F3))%>%
  dplyr::ungroup()%>%
  dplyr::arrange(desc(mf3))%>%
  dplyr::mutate(Outgroup2 = factor(Outgroup,levels = unique(Outgroup), labels = unique(Outgroup)))

outpop <- "B"
midpop1 <- "J"
midpop2 <- "N"

fthree %>%
  dplyr::filter(Outgroup == outpop) %>%
  dplyr::filter(PopA %in%c(midpop1,midpop2) | PopB %in%c(midpop1,midpop2))

hist(fthree$F3)


ggplot(fthree)+
  aes(x = Outgroup2, y = F3, fill = Outgroup)+
  geom_boxplot(alpha = 0.85)+
  theme_bw()+
  scale_fill_manual(values=ancestry.colours)+
  labs(x = "Outgroup", y = "F3 Statistic")+
  theme(axis.text.x = ggplot2::element_text(size = 16),
        axis.text.y = ggplot2::element_text(size = 16),
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3))


ffour <- data.table::fread(glue::glue("{data_dir}/{f4_file}") , sep = " ")%>%
  tidyr::separate(V1, into = c("A","B"), sep = ";")%>%
  tidyr::separate(A, into = c("W","X"), sep = ",")%>%
  tidyr::separate(B, into = c("Y","Z"), sep = ",")%>%
  dplyr::rename(D = V2, SE = V3, Z.value = V4)

plot(f4stats(ffour))

fthree <- data.table::fread(glue::glue("{data_dir}/{f3_file}") , sep = " ")%>%
  tidyr::separate(V1, into = c("W","otherpops"), sep = ";")%>%
  tidyr::separate(otherpops, into = c("X","Z"), sep = ",")%>%
  dplyr::rename(D = V2, SE = V3, Z.value = V4)
plot(f4stats(fthree))





leaves <- c(fthree$Outgroup[1],fthree$PopA[1],fthree$PopB[1])
inner_nodes <- c(glue::glue("{leaves[2]}{leaves[3]}"), glue::glue("{leaves[2]}{leaves[3]}{leaves[1]}"))
edges <- parent_edges(c(admixturegraph::edge(glue::glue("{leaves[2]}"), glue::glue("{leaves[2]}{leaves[3]}")),
                        admixturegraph::edge(glue::glue("{leaves[3]}"), glue::glue("{leaves[2]}{leaves[3]}")),
                        admixturegraph::edge(glue::glue("{leaves[2]}{leaves[3]}"), glue::glue("{leaves[2]}{leaves[3]}{leaves[1]}")),
                        admixturegraph::edge(glue::glue("{leaves[1]}"), glue::glue("{leaves[2]}{leaves[3]}{leaves[1]}"))))
graph <- agraph(leaves, inner_nodes, edges)
plot(graph)
