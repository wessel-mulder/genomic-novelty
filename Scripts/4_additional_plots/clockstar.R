library(ClockstaRX)
library(tidyverse)
library(phytools)
library(svglite)
library(ggvenn)

setwd('/Users/jule/Desktop/genomic-novelty/')

# load data
data <- read_tsv('Data/2_comparison_outliers/metadata_habitat_reptraits.tsv')

# remove old evolutionary rates
data <- data %>% dplyr::select(!contains(c("evolutionary")))

# load species tree
species.tree <- read.nexus('Data/2_comparison_outliers/bd.mcc.median_heights.tre')
species.tree <- drop.tip(species.tree, 
                         setdiff(species.tree$tip.label, data$ID))

# load dn_ds matrices
dNdS_impute <- readRDS('Results/1_mahalanobis_outliers/dNdS_impute.rds')

# load UCE/gene trees
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

analysis <- diagnose.clocks(loctrs = locus.trees, 
                            sptr = species.tree, 
                            pdf.file = 'Results/4_additional_plots/clockstarx')

saveRDS(analysis, file = 'Results/4_additional_plots/clockstar_analysis.rds')


################
### HERE !!! ###
################

analysis <- readRDS("Results/4_additional_plots/clockstar_analysis.rds")
analysis$pca.best.k
# 1


# Extract PC1 for two Euclidean spaces
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

# Calculate the quantile threshold (absolute for two-sided test)
threshold <- quantile(abs(pc2.resid), 0.95)
# Identify outlier indices
outlier_indices <- which(abs(pc2.resid) > threshold)
# Extract outlier matrices
outlier_genes_pc2_resid <- names[outlier_indices]

# Calculate the quantile threshold (absolute for two-sided test)
threshold <- quantile(abs(pc3.resid), 0.95)
# Identify outlier indices
outlier_indices <- which(abs(pc3.resid) > threshold)
# Extract outlier matrices
outlier_genes_pc3_resid <- names[outlier_indices]

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


binom.test(length(intersect(unique(gene_dnds$gene),
                            unique(outlier_genes_pc1_resid))), 
           length(gene_dnds$gene), p = 0.05, alternative='greater')

binom.test(length(intersect(unique(gene_dnds$gene),
                            unique(outlier_genes_pc2_resid))), 
           length(gene_dnds$gene), p = 0.05, alternative='greater')

binom.test(length(intersect(unique(gene_dnds$gene),
                            unique(outlier_genes_pc3_resid))), 
           length(gene_dnds$gene), p = 0.05, alternative='greater')

pc1_length <- length(unique(outlier_genes_pc1_resid))
mahalanobis_length <- length(unique(gene_dnds$gene))
overlap_pc1_length <- length(intersect(unique(gene_dnds$gene),
                                       unique(outlier_genes_pc1_resid)))
input_pc1 <- matrix(c(overlap_pc1_length, 
                      pc1_length - overlap_pc1_length, 
                      mahalanobis_length - overlap_pc1_length, 
                      length(locus.trees) - pc1_length - mahalanobis_length + overlap_pc1_length), 
                    nrow = 2)
input_pc1
fisher.test(input_pc1, alternative = "greater")

