library(ClockstaRX)
library(tidyverse)
library(phytools)
library(svglite)
library(ggvenn)

setwd('/Users/jule/Desktop/genomic-novelty/')

# load data
data <- read_tsv('Data/2_comparison_outliers/metadata_habitat.tsv')

# load species tree
species.tree <- read.nexus('Data/2_comparison_outliers/bd.mcc.median_heights.tre')
species.tree <- drop.tip(species.tree, 
                         setdiff(species.tree$tip.label, data$ID))

# load dn_ds matrices
dNdS_impute <- readRDS('Results/1_mahalanobis_outliers/dNdS_impute.rds')

# create locus trees with species tree topology and branch length defined by dNdS
locus.trees <- list()
for (i in seq_along(dNdS_impute)) {
  locus.tree <- nnls.tree(dm = dNdS_impute[[i]], tree = species.tree)
  locus.trees[[i]] <- locus.tree
}

class(locus.trees) <- "multiPhylo"

# filter that trees have at least 3 edges, meaning at least 2 tips
locus.trees <- locus.trees[Ntip(locus.trees) > 2]


################
### SKIP !!! ###
################

# run clockstaRX
analysis <- diagnose.clocks(loctrs = locus.trees, 
                            sptr = species.tree, 
                            pdf.file = 'Results/4_additional_plots/clockstarx')

# save results to be able to load later
saveRDS(analysis, file = 'Results/4_additional_plots/clockstar_analysis.rds')


################
### HERE !!! ###
################

# load results
analysis <- readRDS("Results/4_additional_plots/clockstar_analysis.rds")

# Extract residual PC1, PC2, PC3
pc1.resid <- analysis$weighted.pca.clock.space$empPCA$x[,1]
pc2.resid <- analysis$weighted.pca.clock.space$empPCA$x[,2]
pc3.resid <- analysis$weighted.pca.clock.space$empPCA$x[,3]
names <- names(dNdS_impute)

# Quantiles

# Calculate the quantile threshold (absolute for two-sided test)
threshold <- quantile(abs(pc1.resid), 0.95)
# Identify outlier indices
outlier_indices <- which(abs(pc1.resid) > threshold)
# Extract outlier matrices
outlier_genes_pc1_resid <- names[outlier_indices]
write.csv(outlier_genes_pc1_resid, 'Results/4_additional_plots/outliers_genes/pc1_resid_q95.csv', row.names = FALSE)

# Calculate the quantile threshold (absolute for two-sided test)
threshold <- quantile(abs(pc2.resid), 0.95)
# Identify outlier indices
outlier_indices <- which(abs(pc2.resid) > threshold)
# Extract outlier matrices
outlier_genes_pc2_resid <- names[outlier_indices]
write.csv(outlier_genes_pc2_resid, 'Results/4_additional_plots/outliers_genes/pc2_resid_q95.csv', row.names = FALSE)

# Calculate the quantile threshold (absolute for two-sided test)
threshold <- quantile(abs(pc3.resid), 0.95)
# Identify outlier indices
outlier_indices <- which(abs(pc3.resid) > threshold)
# Extract outlier matrices
outlier_genes_pc3_resid <- names[outlier_indices]
write.csv(outlier_genes_pc3_resid, 'Results/4_additional_plots/outliers_genes/pc3_resid_q95.csv', row.names = FALSE)

gene_dnds <- read.csv('Results/1_mahalanobis_outliers/outliers_genes/dnds_q95.csv',
                      col.names = c('gene', 'p'))

# Venn diagrams

venn_data_gene <- setNames(
  list(
    unique(gene_dnds$gene),
    unique(outlier_genes_pc1_resid)
  ),
  list("Mahalanobis", "Residual PC1")
)
venn_gene <- ggvenn(venn_data_gene, 
                    stroke_size = 0, set_name_size = 7, text_size = 7, 
                    show_percentage = F, stroke_color = "grey", text_color = "black",
                    fill_color = c("#72315C", "#A6A867")
)
venn_gene
svglite('Plots/4_additional_plots/venn_resid_pc1.svg', width = 7, height = 5)
print(venn_gene)
dev.off()

venn_data_gene <- setNames(
  list(
    unique(gene_dnds$gene),
    unique(outlier_genes_pc2_resid)
  ),
  list("Mahalanobis", "Residual PC2")
)
venn_gene <- ggvenn(venn_data_gene, 
                    stroke_size = 0, set_name_size = 7, text_size = 7, 
                    show_percentage = F, stroke_color = "grey", text_color = "black",
                    fill_color = c("#72315C", "#A6A867")
)
venn_gene
svglite('Plots/4_additional_plots/venn_resid_pc2.svg', width = 7, height = 5)
print(venn_gene)
dev.off()

venn_data_gene <- setNames(
  list(
    unique(gene_dnds$gene),
    unique(outlier_genes_pc3_resid)
  ),
  list("Mahalanobis", "Residual PC3")
)
venn_gene <- ggvenn(venn_data_gene, 
                    stroke_size = 0, set_name_size = 7, text_size = 7, 
                    show_percentage = F, stroke_color = "grey", text_color = "black",
                    fill_color = c("#72315C", "#A6A867")
)
venn_gene
svglite('Plots/4_additional_plots/venn_resid_pc3.svg', width = 7, height = 5)
print(venn_gene)
dev.off()

## Statistical tests

# Fisher's exact test
mahalanobis_length <- length(unique(gene_dnds$gene))

pc1_length <- length(unique(outlier_genes_pc1_resid))
overlap_pc1_length <- length(intersect(unique(gene_dnds$gene),
                                       unique(outlier_genes_pc1_resid)))
input_pc1 <- matrix(c(overlap_pc1_length, 
                      pc1_length - overlap_pc1_length, 
                      mahalanobis_length - overlap_pc1_length, 
                      length(locus.trees) - pc1_length - mahalanobis_length + overlap_pc1_length), 
                    nrow = 2)
fisher.test(input_pc1, alternative = "greater")

pc2_length <- length(unique(outlier_genes_pc2_resid))
overlap_pc2_length <- length(intersect(unique(gene_dnds$gene),
                                       unique(outlier_genes_pc2_resid)))
input_pc2 <- matrix(c(overlap_pc2_length, 
                      pc1_length - overlap_pc2_length, 
                      mahalanobis_length - overlap_pc2_length, 
                      length(locus.trees) - pc2_length - mahalanobis_length + overlap_pc2_length), 
                    nrow = 2)
fisher.test(input_pc2, alternative = "greater")

pc3_length <- length(unique(outlier_genes_pc3_resid))
overlap_pc3_length <- length(intersect(unique(gene_dnds$gene),
                                       unique(outlier_genes_pc3_resid)))
input_pc3 <- matrix(c(overlap_pc3_length, 
                      pc1_length - overlap_pc3_length, 
                      mahalanobis_length - overlap_pc3_length, 
                      length(locus.trees) - pc3_length - mahalanobis_length + overlap_pc3_length), 
                    nrow = 2)
fisher.test(input_pc3, alternative = "greater")

