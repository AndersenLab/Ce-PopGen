library(tidyverse)
library(extrafont)
library(ggbeeswarm)
library(linkagemapping)
library(data.table)
library(genetics)
library(cegwas)
library(maps)
library(ggthemes)
library(ggtree)
# format settings
# colors
axis_color <- "#000F08"
highlight_color <- "#D7263D"
background_color <- "white"
# strain colors
strain_colors <- c("N2" = "#F9A227", "CB4856" = "#2790F9", "Other" = "gray50", "JU775" = "#999999", "DL238" = "hotpink3",
                   "REF"= "gray60", "ALT" = "#D7263D")
col_blind_colors <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

# font
number_font <- "Itim"
axes_text_size <- 20
axes_title_font <- "Montserrat ExtraBold"
axes_title_size <- 20
title_size <- 20

# plot features
boxplot_alpha <- 0.7
point_alpha <- 0.4
point_size <- 2
point_highlight_size <- 4

# import fonts
#font_import()


base_theme <- theme(
  line = element_line(colour = axis_color, size = 0.5, linetype = 1, lineend = "butt"), 
  rect = element_rect(fill = background_color, colour = axis_color, size = 0.5, linetype = 1), 
  text = element_text(family = axes_title_font, size = axes_text_size), 
  
  axis.text = element_text(family = number_font,size = rel(0.8), colour = "grey30", margin = unit(0.1, "cm")),
  strip.text = element_text(size = rel(0.8)), 
  strip.background = element_blank(),
  axis.text.x = element_text(vjust = 1), 
  axis.text.y = element_text(hjust = 1), 
  axis.ticks = element_line(colour = "gray90"), 
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), angle = 90), 
  axis.ticks.length = unit(0.15, "cm"),
  
  legend.background = element_rect(fill = background_color, colour = NA), 
  legend.spacing = unit(0.2, "cm"), 
  legend.key = element_rect(fill = NA, colour = NA), 
  legend.key.size = unit(1.2, "lines"), 
  legend.key.height = NULL, 
  legend.key.width = NULL, 
  legend.text = element_text(size = rel(0.8)), 
  legend.text.align = NULL, 
  legend.title = element_text(size = rel(0.8), hjust = 0), 
  legend.title.align = NULL, 
  legend.position = "right", 
  legend.direction = NULL, 
  legend.justification = "center", 
  legend.box = NULL, 
  
  plot.background = element_rect(fill = background_color),
  
  panel.background = element_rect(fill = background_color, colour = NA), 
  panel.border = element_blank(), 
  panel.grid.major = element_line(colour = "gray90"), 
  panel.grid.minor = element_blank(), 
  panel.spacing = unit(1, "lines"), 
  panel.margin.x = NULL, 
  panel.margin.y = NULL)

