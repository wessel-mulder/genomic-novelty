library(ggtree)
library(dplyr)
library(readr)
library(ape)
library(phytools)


# COPY JULES SCRIPT -------------------------------------------------------
library(tidyverse)
library(fishualize)
library(reshape2) # for melt
library(ape)
library(ggtree)
library(deeptime)
library(phangorn)
library(svglite)
library(extrafont)
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

# read outlier genes 
### GET NUMBER OF DNDS OUTLIERS
outliers_dnds_species <- read.csv("Results/1_mahalanobis_outliers/outliers_species/dnds_chi95.csv")[,-1]
colnames(outliers_dnds_species) <- c("ID", "gene")
# n_genes
n_genes <- outliers_dnds_species %>% 
  dplyr::select(ID, gene) %>%
  group_by(ID) %>%
  summarise(unique_outliers = n_distinct(gene), .groups = "drop")
rm(outliers_dnds_species)
# merge for species names into metaturtles
meta_turtles <- merge(meta_turtles,n_genes,by='ID')






#########################
### SPECIES TREE PLOT ###
#########################

# load species tree with branch lengths from Thomson et al. (2021)
species_tree_plot <- read.nexus("Data/2_comparison_outliers/bd.mcc.median_heights.tre")
species_tree_plot <- drop.tip(species_tree_plot, 
                              setdiff(species_tree_plot$tip.label, meta_turtles$ID))
species_tree_plot$tip.label <- meta_turtles$Species[
  match(species_tree_plot$tip.label, meta_turtles$ID)]

metadata <- data.frame(taxa=meta_turtles$Species, Outliers=meta_turtles$un)

# plot for our tree
library('ggtree')
# Ensure the Outliers are assigned correctly by species
species_tree_plot$Outliers <- metadata$Outliers[match(species_tree_plot$tip.label, metadata$taxa)]

metadata$tip.label <- metadata$taxa

# Start plotting
plot_tree <- ggtree(species_tree_plot) %<+% metadata +
  geom_tiplab(aes(label = tip.label), size = 3, offset = 2, fontface = "italic") +
  geom_tippoint(aes(color = Outliers), size = 3) +
  scale_color_continuous(name = 'dN/dS', limits = c(0, 1.5),
                         oob = scales::squish,
                         low = 'darkgreen', high = 'red') +
  theme_tree2() +
  theme(text = element_text(family = "Arial"))

plot_tree 
revts(plot_tree)

  # make points a bit bigger
  scale_color_viridis_c(option = "C") +             # continuous color scale, nice for numeric
  theme_tree2() +
  coord_geo(xlim = c(-250, 100), ylim = c(-0.5, Ntip(species_tree_plot)+2),
            neg = TRUE, abbrv = list(FALSE), dat = list("periods"),
            fill = c("grey90", "grey85", "grey80", "grey75", "grey70", "grey65"),
            color = "black",
            lab_color = "black",
            pos = list("bottom"), size = "auto",
            height = list(unit(1, "lines"))) +
  scale_x_continuous(breaks = seq(-240, 0, 20), labels = abs(seq(-240, 0, 20))) +
  theme(panel.grid.major = element_line(color = "grey80", size = 0.2),
        panel.grid.major.y = element_blank(),
        text = element_text(family = "Arial"))


plot_tree <- plot %<+% metadata +
  geom_tiplab(size = 3, offset = 2, fontface = "italic") + 
  geom_tippoint(aes(color = Outliers), size = 3) +  # points at tips
  geom_text2(aes(label = Outliers), hjust = -0.3, size = 3, color = "black") +  # values on branches
  scale_color_viridis_c(option = "C") +
  theme_tree2() +
  coord_geo(xlim = c(-250, 100), ylim = c(-0.5, Ntip(species_tree_plot)+2),
            neg = TRUE, abbrv = list(FALSE), dat = list("periods"),
            fill = c("grey90", "grey85", "grey80", "grey75", "grey70", "grey65"),
            color = "black",
            lab_color = "black",
            pos = list("bottom"), size = "auto",
            height = list(unit(1, "lines"))) +
  scale_x_continuous(breaks = seq(-240, 0, 20), labels = abs(seq(-240, 0, 20))) +
  theme(panel.grid.major = element_line(color = "grey80", size = 0.2),
        panel.grid.major.y = element_blank(),
        text = element_text(family = "Arial"))

revts(plot_tree)


# LOADING DATA  -----------------------------------------------------------
outliers_vec <- setNames(meta_turtles$unique_outliers, meta_turtles$Species)

#test
phylosig(species_tree_plot,outliers_vec,test=T)
phylosig(species_tree_plot,outliers_vec,method='lambda',test=T)

contMap(species_tree_plot, outliers_vec)   # plots trait along the tree branches



# TRYING TO PLOT ----------------------------------------------------------


# function to plot probability or trait value by branch
# written by Liam J. Revell 2013
plotBranchbyTrait(species_tree_plot,outliers_vec,method='tips')

plotBranchbyTrait<-function(tree,x,mode=c("edges","tips","nodes"),palette="rainbow",legend=TRUE,xlims=NULL,...){
  mode<-mode[1]
  if(mode=="tips"){
    x<-c(x[tree$tip.label],fastAnc(tree,x))
    names(x)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  } else if(mode=="nodes"){
    XX<-matrix(x[tree$edge],nrow(tree$edge),2)
    x<-rowMeans(XX)
  }
  if(hasArg(tol)) tol<-list(...)$tol
  else tol<-1e-6
  if(palette=="heat.colors") cols<-heat.colors(n=1000)
  if(palette=="gray") cols<-gray(1000:1/1000)
  if(palette=="rainbow")	cols<-rainbow(1000,start=0.7,end=0) # blue->red
  if(is.null(xlims)) xlims<-range(x)+c(-tol,tol)
  breaks<-0:1000/1000*(xlims[2]-xlims[1])+xlims[1]
  whichColor<-function(p,cols,breaks){
    i<-1
    while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
    cols[i]
  }
  colors<-sapply(x,whichColor,cols=cols,breaks=breaks)
  par(lend=2)
  plot.phylo(tree,no.margin=TRUE,edge.width=4,edge.color=colors,label.offset=0.01*max(nodeHeights(tree)),lend=2,new=FALSE,...)
  if(legend==TRUE&&is.logical(legend)) legend<-round(0.3*max(nodeHeights(tree)),2)
  if(legend){
    if(hasArg(title)) title<-list(...)$title
    else title<-NULL
    if(hasArg(digits)) digits<-list(...)$digits
    else digits<-1
    add.color.bar(legend,cols,title,xlims,digits,prompt=FALSE)
  }
}

# function to add color bar
# written by Liam J. Revell 2013

add.color.bar<-function(leg,cols,title=NULL,lims=c(0,1),digits=1,prompt=TRUE){
  cat("Click where you want to draw the bar\n")
  x<-unlist(locator(1))
  y<-x[2]
  x<-x[1]
  X<-x+cbind(0:(length(cols)-1)/length(cols),1:length(cols)/length(cols))*(leg)
  Y<-cbind(rep(y,length(cols)),rep(y,length(cols))) 		
  lines(c(X[1,1],X[nrow(X),2]),c(Y[1,1],Y[nrow(Y),2]),lwd=4+2,lend=2) 
  for(i in 1:length(cols)) lines(X[i,],Y[i,],col=cols[i],lwd=4,lend=2)
  text(x=x,y=y,round(lims[1],digits),pos=3)
  text(x=x+leg,y=y,round(lims[2],digits),pos=3)
  if(is.null(title)) title<-"P(state=1)"
  text(x=(2*x+leg)/2,y=y,title,pos=3)
  text(x=(2*x+leg)/2,y=y,paste("length=",round(leg,3),sep=""),pos=1)
}
