input <- 'Results/1_mahalanobis_outliers'

df_allgenes <- read.csv(file.path(input,'distances_allgenes/distances_dnds.csv'))
df_outliers <- read.csv(file.path(input,'outliers_genes/outliers_dnds_q95.csv'))
breaks <- seq(0,1,by=0.0005)
hist(df_allgenes$pval,
     breaks=breaks,
     xlim=c(0,0.05))

