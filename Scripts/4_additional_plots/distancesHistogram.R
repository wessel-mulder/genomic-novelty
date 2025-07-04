input <- 'Results/1_mahalanobis_outliers'
output <- 'Plots/4_additional_plots/'

df_allgenes <- read.csv(file.path(input,'distances_allgenes/dnds.csv'))
df_outliers <- read.csv(file.path(input,'outliers_genes/dnds_q95.csv'))
breaks <- seq(0,1,by=0.005)
hist(df_allgenes$pval,
     breaks=breaks,
     xlim=c(0,0.05))
abline(v=0.05,col='red')

threshold <- quantile(df_allgenes$distance,0.95)

pdf(file.path(output,'dnds_distances.pdf'),
    width = 8,
    height = 4)
hist(df_allgenes$distance,
     main = NULL,
     xlab = 'Mahalanobis distance from multivariate distribution of pairwise comparisons across taxa')
abline(v=threshold,col='red')
dev.off()


