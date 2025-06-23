library(tidyverse)
library(ape)
library(ggtree)
library(svglite)
library(extrafont)
# import arial (only needed once)
# font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")

setwd("/Users/jule/Desktop/genomic-novelty/")

# get data on used species
meta_turtles <- read_tsv("Data/2_comparison_outliers/metadata_habitat_reptraits.tsv")
meta_turtles$Habitat_factor <- factor(meta_turtles$Microhabitat, 
                                      levels=c("Marine", "Aquatic", 
                                               "Aquatic_Terrestrial", "Terrestrial", 
                                               "Outgroup"))
meta_turtles <- meta_turtles %>% filter(Microhabitat != "Outgroup")

# load species tree with branch lengths from Thomson et al. (2021)
species_tree_plot <- read.nexus("Data/2_comparison_outliers/bd.mcc.median_heights.tre")
species_tree_plot <- drop.tip(species_tree_plot, 
                              setdiff(species_tree_plot$tip.label, meta_turtles$ID))
species_tree_plot$tip.label <- meta_turtles$Species[
  match(species_tree_plot$tip.label, meta_turtles$ID)]

### BUSCO + N50 quality check

stats <- read_tsv("Data/4_additional_plots/Assembly_stats.tsv")
stats$`BUSCO S` <- stats$`BUSCO S` * 100

# plot for our tree
plot_tree <- ggtree::ggtree(species_tree_plot) + ggtree::xlim_tree(950)
plot_tree <- plot_tree +
  ggtree::geom_tiplab(size=4, offset=4, fontface = "italic") + 
  theme_tree2() +
  vexpand(0.01, direction = -1)
plot_tree

facet_stats <- facet_plot(plot_tree,
                          data=stats,
                          geom=geom_bar,
                          mapping = aes(x=`Assembly length (Gb)`, fill=`Assembly length (Gb)`),
                          stat="identity",
                          orientation='y',
                          panel="Assembly length (Gb)",
                          position = position_stack(reverse = TRUE),
                          show.legend = F)

facet_stats <- facet_plot(facet_stats,
                          data=stats,
                          geom=geom_bar,
                          mapping = aes(x=`Scaffold N50 (Mb)`, fill=`Scaffold N50 (Mb)`),
                          stat="identity",
                          orientation='y',
                          panel="Scaffold N50 (Mb)",
                          position = position_stack(reverse = TRUE),
                          show.legend = F)

facet_stats <- facet_plot(facet_stats,
                          data=stats,
                          geom=geom_bar,
                          mapping = aes(x=`BUSCO S`, fill=`BUSCO S`),
                          stat="identity",
                          orientation='y',
                          panel="BUSCO S (%)",
                          position = position_stack(reverse = TRUE),
                          show.legend = F)

facet_stats + theme(text = element_text(family = "Arial"))

facet_stats

svglite("Plots/4_additional_plots/assembly_stats.svg", width = 11, height = 5)
print(facet_stats)
dev.off()

ggsave("Plots/4_additional_plots/assembly_stats.pdf", width = 11, height = 5)
ggsave("Plots/4_additional_plots/assembly_stats.png", width = 11, height = 5)
