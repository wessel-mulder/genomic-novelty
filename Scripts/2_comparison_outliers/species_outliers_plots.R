library(tidyverse)
library(fishualize)
library(reshape2) # for melt
library(ape)
library(ggtree)
library(deeptime)
library(phangorn)
library(svglite)
library(extrafont)
library(Cairo)
# import arial (only needed once)
# font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")

setwd("/Users/jule/Desktop/genomic-novelty/")

colours_classes <- fish(n=1,option="Balistoides_conspicillum", end=0.9, 
                        begin=0.9)
colours_classes[2] <- fish(n=1,option="Balistoides_conspicillum", end=0.8, 
                           begin=0.8)
colours_classes[3] <- fish(n=1,option="Balistoides_conspicillum", end=0.4, 
                           begin=0.4)
colours_classes[4] <- fish(n=1,option="Balistoides_conspicillum", end=0.2, 
                           begin=0.2)

gene_names <- readLines("Results/0_preprocessing/list_final.txt")

outliers_dn_species <- read.csv("Results/1_mahalanobis_outliers/outliers_species/dn_chi95.csv")
colnames(outliers_dn_species) <- c("n", "species", "gene")
outliers_ds_species <- read.csv("Results/1_mahalanobis_outliers/outliers_species/ds_chi95.csv")
colnames(outliers_ds_species) <- c("n", "species", "gene")
outliers_dnds_species <- read.csv("Results/1_mahalanobis_outliers/outliers_species/dnds_chi95.csv")
colnames(outliers_dnds_species) <- c("n", "species", "gene")
outliers_raw_species <- read.csv("Results/1_mahalanobis_outliers/outliers_species/raw_chi95.csv")
colnames(outliers_raw_species) <- c("n", "species", "gene")

meta_turtles <- read_tsv("Data/2_comparison_outliers/metadata_habitat_reptraits.tsv")
meta_turtles$Habitat_factor <- factor(meta_turtles$Microhabitat, 
                                      levels=c("Marine", "Aquatic", 
                                               "Aquatic_Terrestrial", "Terrestrial", 
                                               "Outgroup"))
meta_turtles <- meta_turtles %>% filter(Microhabitat != "Outgroup")


#########################
### SPECIES TREE PLOT ###
#########################

# load species tree with branch lengths from Thomson et al. (2021)
species_tree_plot <- read.nexus("Data/2_comparison_outliers/bd.mcc.median_heights.tre")
species_tree_plot <- drop.tip(species_tree_plot, 
                              setdiff(species_tree_plot$tip.label, meta_turtles$ID))
species_tree_plot$tip.label <- meta_turtles$Species[
  match(species_tree_plot$tip.label, meta_turtles$ID)]

# plot for our tree
plot_tree <- ggtree::ggtree(species_tree_plot)
plot_tree <- plot_tree +
  ggtree::geom_tiplab(size=3, offset=0.5, fontface = "italic") + 
  theme_tree2() +
  coord_geo(xlim = c(-250, 100), ylim = c(-0.5, Ntip(species_tree_plot)+2),
            neg = TRUE, abbrv = list(FALSE), dat=list("periods"),
            fill = c("grey90", "grey85", "grey80", "grey75", "grey70", "grey65"),
            color= "black",
            lab_color = "black",
            pos = list("bottom"), size = "auto",
            height = list(unit(1, "lines"))) +
  scale_x_continuous(breaks = seq(-240, 0, 20), labels = abs(seq(-240, 0, 20))) +
  theme(panel.grid.major   = element_line(color="grey80", size=.2),
        panel.grid.major.y = element_blank(),
        text = element_text(family = "Arial2"))
revts(plot_tree)

svglite('Plots/2_comparison_outliers/species_tree_branch_lengths.svg', width = 8, height = 5)
revts(plot_tree)
dev.off()

CairoPDF("Plots/2_comparison_outliers/species_tree_branch_lengths.pdf", width = 8, height = 5)
revts(plot_tree)
dev.off()
ggsave("Plots/2_comparison_outliers/species_tree_branch_lengths.pdf", width = 8, height = 5)


###############################
### OVERLAP BETWEEN SPECIES ###
###############################

species_list <- meta_turtles$ID

df_species_id <- meta_turtles %>% dplyr::select(c("Species", "ID"))

species_ordered_list <- c("Malaclemys terrapin",
                          "Chrysemys picta",
                          "Terrapene mexicana",
                          "Platysternon megacephalum",
                          "Cuora mccordi",
                          "Cuora amboinensis",
                          "Gopherus agassizii",
                          "Chelonoidis abingdonii",
                          "Dermochelys coriacea",
                          "Chelonia mydas",
                          "Dermatemys mawii",
                          "Chelydra serpentina",
                          "Carettochelys insculpta",
                          "Pelodiscus sinensis",
                          "Pelusios castaneus",
                          "Podocnemis expansa",
                          "Emydura subglobosa",
                          "Mesoclemmys tuberculata")

# DN
data_heatmap_dn <- outliers_dn_species %>% dplyr::select(c("species", "gene"))

matrix_heatmap_dn <- matrix(0, nrow = length(species_list), ncol = length(species_list))
rownames(matrix_heatmap_dn) <- species_list
colnames(matrix_heatmap_dn) <- species_list

for (i in 1:length(species_list)) {
  for (j in 1:length(species_list)) {
    genes_i <- data_heatmap_dn$gene[data_heatmap_dn$species == species_list[i]]
    genes_j <- data_heatmap_dn$gene[data_heatmap_dn$species == species_list[j]]
    matrix_heatmap_dn[i, j] <- length(intersect(genes_i, genes_j))
    
    # if same species
    # if (i == j) {
    #   other_genes <- data_heatmap_dn$gene[data_heatmap_dn$species != species_list[i]]
    #   unique_genes <- setdiff(genes_i, other_genes)
    #   matrix_heatmap_dn[i, j] <- length(unique_genes)
    # }
  }
}

df_heatmap_dn <- melt(matrix_heatmap_dn)
df_heatmap_dn <- merge(df_heatmap_dn, df_species_id, by.x="Var1", by.y="ID", all.x=T)
colnames(df_heatmap_dn) <- c("Var1", "Var2", "value", "Species1")
df_heatmap_dn <- merge(df_heatmap_dn, df_species_id, by.x="Var2", by.y="ID", all.x=T)
colnames(df_heatmap_dn) <- c("Var2", "Var1", "value", "Species1", "Species2")
df_heatmap_dn$Species1 <- factor(df_heatmap_dn$Species1, 
                                 levels=species_ordered_list)
df_heatmap_dn$Species2 <- factor(df_heatmap_dn$Species2, 
                                 levels=species_ordered_list[rev(1:length(species_ordered_list))])

heatmap_dn <- ggplot(df_heatmap_dn, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value), color = "white") +
  labs(title = "dN",
       x = "Species",
       y = "Species") +
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#72315C", "#A6A867")) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = 'italic')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = 'italic')) +
  theme(text = element_text(family = "Arial"))

svglite('Plots/2_comparison_outliers/heatmap_dn.svg', width = 8, height = 5)
print(heatmap_dn)
dev.off()

heatmap_dn
ggsave("Plots/2_comparison_outliers/heatmap_dn.pdf", width = 8, height = 5)


# DS
data_heatmap_ds <- outliers_ds_species %>% dplyr::select(c("species", "gene"))

matrix_heatmap_ds <- matrix(0, nrow = length(species_list), ncol = length(species_list))
rownames(matrix_heatmap_ds) <- species_list
colnames(matrix_heatmap_ds) <- species_list

for (i in 1:length(species_list)) {
  for (j in 1:length(species_list)) {
    genes_i <- data_heatmap_ds$gene[data_heatmap_ds$species == species_list[i]]
    genes_j <- data_heatmap_ds$gene[data_heatmap_ds$species == species_list[j]]
    matrix_heatmap_ds[i, j] <- length(intersect(genes_i, genes_j))
    
    # if same species
    # if (i == j) {
    #   other_genes <- data_heatmap_ds$gene[data_heatmap_ds$species != species_list[i]]
    #   unique_genes <- setdiff(genes_i, other_genes)
    #   matrix_heatmap_ds[i, j] <- length(unique_genes)
    # }
  }
}

df_heatmap_ds <- melt(matrix_heatmap_ds)
df_heatmap_ds <- merge(df_heatmap_ds, df_species_id, by.x="Var1", by.y="ID", all.x=T)
colnames(df_heatmap_ds) <- c("Var1", "Var2", "value", "Species1")
df_heatmap_ds <- merge(df_heatmap_ds, df_species_id, by.x="Var2", by.y="ID", all.x=T)
colnames(df_heatmap_ds) <- c("Var2", "Var1", "value", "Species1", "Species2")
df_heatmap_ds$Species1 <- factor(df_heatmap_ds$Species1, 
                                 levels=species_ordered_list)
df_heatmap_ds$Species2 <- factor(df_heatmap_ds$Species2, 
                                 levels=species_ordered_list[rev(1:length(species_ordered_list))])

heatmap_ds <- ggplot(df_heatmap_ds, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value), color = "white") +
  labs(title = "dS",
       x = "Species",
       y = "Species") +
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#72315C", "#A6A867")) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = 'italic')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = 'italic')) +
  theme(text = element_text(family = "Arial"))

svglite('Plots/2_comparison_outliers/heatmap_ds.svg', width = 8, height = 5)
print(heatmap_ds)
dev.off()

heatmap_ds
ggsave("Plots/2_comparison_outliers/heatmap_ds.pdf", width = 8, height = 5)


# DNDS RELATIVE
data_heatmap_dnds <- outliers_dnds_species %>% dplyr::select(c("species", "gene"))

matrix_heatmap_dnds <- matrix(0, nrow = length(species_list), ncol = length(species_list))
rownames(matrix_heatmap_dnds) <- species_list
colnames(matrix_heatmap_dnds) <- species_list

for (i in 1:length(species_list)) {
  for (j in 1:length(species_list)) {
    genes_i <- data_heatmap_dnds$gene[data_heatmap_dnds$species == species_list[i]]
    genes_j <- data_heatmap_dnds$gene[data_heatmap_dnds$species == species_list[j]]
    matrix_heatmap_dnds[i, j] <- round(length(intersect(genes_i, genes_j)) / length(genes_i), digits=2)

    # if same species
    # if (i == j) {
    #   other_genes <- data_heatmap_dnds$gene[data_heatmap_dnds$species != species_list[i]]
    #   unique_genes <- setdiff(genes_i, other_genes)
    #   matrix_heatmap_dnds[i, j] <- length(unique_genes)
    # }
  }
}

df_heatmap_dnds <- melt(matrix_heatmap_dnds)
df_heatmap_dnds <- merge(df_heatmap_dnds, df_species_id, by.x="Var1", by.y="ID", all.x=T)
colnames(df_heatmap_dnds) <- c("Var1", "Var2", "value", "Species1")
df_heatmap_dnds <- merge(df_heatmap_dnds, df_species_id, by.x="Var2", by.y="ID", all.x=T)
colnames(df_heatmap_dnds) <- c("Var2", "Var1", "value", "Species1", "Species2")

df_heatmap_dnds$Species1 <- factor(df_heatmap_dnds$Species1, 
                                   levels=species_ordered_list)
df_heatmap_dnds$Species2 <- factor(df_heatmap_dnds$Species2, 
                                   levels=species_ordered_list[rev(1:length(species_ordered_list))])

heatmap_dnds <- ggplot(df_heatmap_dnds, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value), color = "white") +
  labs(title = "dN/dS",
       x = "Species",
       y = "Species") +
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#72315C", "#A6A867")) +
  # scale_fill_gradientn(name = "Overlapping genes", colours=c("#632A50", "#BDBF09")) +
  # scale_fill_gradientn(name = "Overlapping genes", colours=c("#632A50", "#60935D")) +
  # scale_fill_gradientn(name = "Overlapping genes", colours=c("#632A50", "#188FA7")) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = 'italic')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = 'italic')) +
  theme(text = element_text(family = "Arial"))

svglite('Plots/2_comparison_outliers/heatmap_dnds_relative.svg', width = 12, height = 7.5)
print(heatmap_dnds)
dev.off()

heatmap_dnds
ggsave("Plots/2_comparison_outliers/heatmap_dnds_relative.pdf", width = 16, height = 10)

# DNDS
data_heatmap_dnds <- outliers_dnds_species %>% dplyr::select(c("species", "gene"))

matrix_heatmap_dnds <- matrix(0, nrow = length(species_list), ncol = length(species_list))
rownames(matrix_heatmap_dnds) <- species_list
colnames(matrix_heatmap_dnds) <- species_list

for (i in 1:length(species_list)) {
  for (j in 1:length(species_list)) {
    genes_i <- data_heatmap_dnds$gene[data_heatmap_dnds$species == species_list[i]]
    genes_j <- data_heatmap_dnds$gene[data_heatmap_dnds$species == species_list[j]]
    matrix_heatmap_dnds[i, j] <- length(intersect(genes_i, genes_j))
    
    # if same species
    # if (i == j) {
    #   other_genes <- data_heatmap_dnds$gene[data_heatmap_dnds$species != species_list[i]]
    #   unique_genes <- setdiff(genes_i, other_genes)
    #   matrix_heatmap_dnds[i, j] <- length(unique_genes)
    # }
  }
}

df_heatmap_dnds <- melt(matrix_heatmap_dnds)
df_heatmap_dnds <- merge(df_heatmap_dnds, df_species_id, by.x="Var1", by.y="ID", all.x=T)
colnames(df_heatmap_dnds) <- c("Var1", "Var2", "value", "Species1")
df_heatmap_dnds <- merge(df_heatmap_dnds, df_species_id, by.x="Var2", by.y="ID", all.x=T)
colnames(df_heatmap_dnds) <- c("Var2", "Var1", "value", "Species1", "Species2")

df_heatmap_dnds$Species1 <- factor(df_heatmap_dnds$Species1, 
                                   levels=species_ordered_list)
df_heatmap_dnds$Species2 <- factor(df_heatmap_dnds$Species2, 
                                   levels=species_ordered_list[rev(1:length(species_ordered_list))])

heatmap_dnds <- ggplot(df_heatmap_dnds, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value), color = "white") +
  labs(title = "dN/dS",
       x = "Species",
       y = "Species") +
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#72315C", "#A6A867")) +
  # scale_fill_gradientn(name = "Overlapping genes", colours=c("#632A50", "#BDBF09")) +
  # scale_fill_gradientn(name = "Overlapping genes", colours=c("#632A50", "#60935D")) +
  # scale_fill_gradientn(name = "Overlapping genes", colours=c("#632A50", "#188FA7")) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = 'italic')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = 'italic')) +
  theme(text = element_text(family = "Arial"))

svglite('Plots/2_comparison_outliers/heatmap_dnds.svg', width = 8, height = 5)
print(heatmap_dnds)
dev.off()

heatmap_dnds
ggsave("Plots/2_comparison_outliers/heatmap_dnds.pdf", width = 8, height = 5)

# RAW
data_heatmap_raw <- outliers_raw_species %>% dplyr::select(c("species", "gene"))

matrix_heatmap_raw <- matrix(0, nrow = length(species_list), ncol = length(species_list))
rownames(matrix_heatmap_raw) <- species_list
colnames(matrix_heatmap_raw) <- species_list

for (i in 1:length(species_list)) {
  for (j in 1:length(species_list)) {
    genes_i <- data_heatmap_raw$gene[data_heatmap_raw$species == species_list[i]]
    genes_j <- data_heatmap_raw$gene[data_heatmap_raw$species == species_list[j]]
    matrix_heatmap_raw[i, j] <- length(intersect(genes_i, genes_j))
    
    # if same species
    # if (i == j) {
    #   other_genes <- data_heatmap_raw$gene[data_heatmap_raw$species != species_list[i]]
    #   unique_genes <- setdiff(genes_i, other_genes)
    #   matrix_heatmap_raw[i, j] <- length(unique_genes)
    # }
  }
}

df_heatmap_raw <- melt(matrix_heatmap_raw)
df_heatmap_raw <- merge(df_heatmap_raw, df_species_id, by.x="Var1", by.y="ID", all.x=T)
colnames(df_heatmap_raw) <- c("Var1", "Var2", "value", "Species1")
df_heatmap_raw <- merge(df_heatmap_raw, df_species_id, by.x="Var2", by.y="ID", all.x=T)
colnames(df_heatmap_raw) <- c("Var2", "Var1", "value", "Species1", "Species2")

df_heatmap_raw$Species1 <- factor(df_heatmap_raw$Species1, 
                                   levels=species_ordered_list)
df_heatmap_raw$Species2 <- factor(df_heatmap_raw$Species2, 
                                   levels=species_ordered_list[rev(1:length(species_ordered_list))])

heatmap_raw <- ggplot(df_heatmap_raw, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile() +
  geom_text(aes(label=value), color = "white") +
  labs(title = "RAW",
       x = "Species",
       y = "Species") +
  scale_fill_gradientn(name = "Overlapping genes", colours=c("#72315C", "#A6A867")) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = 'italic')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = 'italic')) +
  theme(text = element_text(family = "Arial"))

svglite('Plots/2_comparison_outliers/heatmap_raw.svg', width = 8, height = 5)
print(heatmap_raw)
dev.off()

heatmap_raw
ggsave("Plots/2_comparison_outliers/heatmap_raw.pdf", width = 8, height = 5)


################
### DOTPLOTS ###
################

### PAIRWISE DIVERGENCE TIMES

# just tips
div_times <- cophenetic.phylo(species_tree_plot)
div_times[lower.tri(div_times)] <- NA

df_div_times <- melt(div_times)
colnames(df_div_times) <- c("Species1", "Species2", "divergence_time")
# change from distance to divergence time
df_div_times <- df_div_times %>% filter(!is.na(divergence_time))
df_div_times$divergence_time <- df_div_times$divergence_time/2
df_div_times$Species_pair <- with(df_div_times, paste0(Species1,", ", Species2))

df_heatmap_dnds$Species_pair <- with(df_heatmap_dnds, paste0(Species1,", ", Species2))
df_heatmap_dnds_reduced <- df_heatmap_dnds%>% 
  dplyr::select(Species_pair, value)
df_div_times <- merge(df_div_times, df_heatmap_dnds_reduced, by="Species_pair")
df_div_times <- df_div_times %>% filter(Species1 != Species2)

dotplot_tips_divtimes_overlaps <- ggplot(df_div_times, 
                                         aes(x=divergence_time, y=value, label=Species_pair)) +
  geom_point(colour="#72315C") +
  theme_minimal() +
  labs(title = "Species pairs",
       x = "Divergence time",
       y = "Number of overlapping genes") +
  theme(text = element_text(family = "Arial"))

svglite('Plots/2_comparison_outliers/dotplot_tips_divtimes_overlaps.svg', width = 8, height = 5)
print(dotplot_tips_divtimes_overlaps)
dev.off()

dotplot_tips_divtimes_overlaps
ggsave("Plots/2_comparison_outliers/dotplot_tips_divtimes_overlaps.pdf", width = 8, height = 5)


### INTERNAL NODES

outliers_dnds_species <- merge(outliers_dnds_species, df_species_id, by.x="species", by.y="ID", all.x=T)
# divergence times all nodes
div_times_nodes <- dist.nodes(species_tree_plot)[19,]
div_times_nodes <- round(max(div_times_nodes) - div_times_nodes, digits = 4)

# get all internal nodes
internal_nodes <- (length(species_tree_plot$tip.label) + 1):max(species_tree_plot$edge)

get_tips <- function(tree, node) {
  tips <- Descendants(tree, node, type = "tips")[[1]]
  # return species names for all tips
  return(tree$tip.label[tips])
}

get_overlapping_genes <- function(tips, genes_species) {
  tip_genes <- genes_species %>% filter(Species %in% tips)
  split_gene_lists <- split(tip_genes$gene, tip_genes$Species)
  overlapping_genes <- Reduce(intersect, split_gene_lists)
  return(length(overlapping_genes))
}

results <- data.frame(node = integer(), divergence_time = numeric(), 
                      num_overlapping_genes = integer(), num_tips = integer())

for (node in internal_nodes) {
  tips <- get_tips(species_tree_plot, node)
  overlapping_genes <- get_overlapping_genes(tips, outliers_dnds_species)
  results <- rbind(results, data.frame(
    node = node,
    divergence_time = div_times_nodes[node],
    num_overlapping_genes = overlapping_genes,
    num_tips = length(tips)
  ))
}

species_richness <- data.frame(ltt.plot.coords(species_tree_plot))
species_richness$time <- round(species_richness$time * -1, digits = 4)
data_plot <- merge(results, species_richness, 
                   by.x = "divergence_time", by.y = "time",
                   all.x = TRUE, all.y = TRUE)


dotplot_internal_divtimes_overlaps <- ggplot(data_plot, 
                                             aes(x=divergence_time, y=num_overlapping_genes)) +
  geom_step(aes(x=divergence_time, y=N), color="grey") +
  geom_point(color="#72315C") +
  scale_y_continuous(
    name = "Number of outlier genes present in all associated tips",
    sec.axis = sec_axis(~ . , name = "Species Richness")
  ) +
  theme_minimal() +
  scale_x_reverse(breaks = seq(0, 240, 20), labels = abs(seq(0, 240, 20)),
                  limits = c(240, 0)) +
  theme(panel.grid.major   = element_line(color="grey80", size=.2),
        panel.grid.minor.x = element_blank()) +
  labs(title = "Internal nodes",
       x = "Divergence time") +
  theme(text = element_text(family = "Arial"))

svglite('Plots/2_comparison_outliers/dotplot_internal_divtimes_overlaps.svg', width = 8, height = 5)
print(dotplot_internal_divtimes_overlaps)
dev.off()

dotplot_internal_divtimes_overlaps
ggsave("Plots/2_comparison_outliers/dotplot_internal_divtimes_overlaps.pdf", width = 8, height = 5)

