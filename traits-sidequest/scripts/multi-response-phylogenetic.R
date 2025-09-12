# GETTING STARTED ---------------------------------------------------------
# packages
pacman::p_load(tidyverse,brms,MCMCglmm,parallel,
               knitr,geiger,ape,phytools,caper,bayesplot,
               ggtree,gridExtra,MCMCpack,dplyr)

# tutorial functions file
source("MR-PMM_tutorial_functions.R")

### 

# LOADING DATA  -----------------------------------------------------------
# TURTLE METADATA
meta_turtles <- read_tsv("Data/2_comparison_outliers/metadata_habitat_reptraits.tsv")
meta_turtles$Habitat_factor <- factor(meta_turtles$Microhabitat, 
                                      levels=c("Marine", "Aquatic", 
                                               "Aquatic_Terrestrial", "Terrestrial", 
                                               "Outgroup"))
meta_turtles <- meta_turtles %>% filter(Microhabitat != "Outgroup")

### GET NUMBER OF DNDS OUTLIERS
outliers_dnds_species <- read.csv("Results/1_mahalanobis_outliers/outliers_species/dnds_chi95.csv")[,-1]
colnames(outliers_dnds_species) <- c("ID", "gene")
# n_genes
n_genes <- outliers_dnds_species %>% 
  dplyr::select(ID, gene) %>%
  group_by(ID) %>%
  summarise(unique_outliers = n_distinct(gene), .groups = "drop")
rm(outliers_dnds_species)
# merge for species names 
n_genes <- left_join(n_genes,
                     meta_turtles[,colnames(meta_turtles) %in% c('ID','Species')],
                     by='ID')


### PREPARE PHYLOGONY
species_tree_plot <- read.nexus("Data/2_comparison_outliers/bd.mcc.median_heights.tre")
species_tree_plot <- drop.tip(species_tree_plot, 
                              setdiff(species_tree_plot$tip.label, meta_turtles$ID))
species_tree_plot$tip.label <- meta_turtles$Species[
  match(species_tree_plot$tip.label, meta_turtles$ID)]
# plot for sanity
plot(species_tree_plot)
# get cov matrix 
C <- ape::vcv.phylo(species_tree_plot,corr=T)

### GET TRAITS 
traits <- readxl::read_xlsx('traits-sidequest/data/CheloniansTraits dataset.xlsx',
                            sheet ='CheloniansTraits Data',
                            skip=1)
#traits <- traits[traits$`In our set`==1,]
setdiff(species_tree_plot$tip.label,traits$Species)
# 1 mismatch, subspecies vs species name 
traits$Species[traits$Species=='Chelonoidis niger'] <- 'Chelonoidis abingdonii'
setdiff(species_tree_plot$tip.label,traits$Species)
traits_subset <- traits[traits$Species%in%species_tree_plot$tip.label,]


#rm(traits,traits_ours)
# no mismatches 
# tidy & subset
colnames(traits_subset) <- gsub(" ", "_",
                                tolower(gsub("(.)([A-Z])", "\\1 \\2", 
                                             colnames(traits_subset))))
colnames(traits) <- gsub(" ", "_",
                                tolower(gsub("(.)([A-Z])", "\\1 \\2", 
                                             colnames(traits))))

# IMPUTE TRAITS  ------------------------------------------------------------------
library(dplyr)
col <- 'maximum_mass_g'
df <- traits
# numeric
impute_numeric_by_genus <- function(df, col, by = "genus") {
  col_sym <- sym(col)
  by_sym <- sym(by)
  
  df %>%
    group_by(!!by_sym) %>%
    mutate(!!col_sym := if_else(is.na(!!col_sym),
                                mean(!!col_sym, na.rm = TRUE),
                                !!col_sym)) %>%
    ungroup()
}

# # Function to impute missing values in a categorical column by the mode within genus (or any other group)
# impute_categorical_by_genus <- function(df, col, by = "genus") {
#   col_sym <- sym(col)
#   by_sym <- sym(by)
#   
#   mode_fn <- function(x) {
#     ux <- na.omit(x)
#     if (length(ux) == 0) return(NA)
#     ux[which.max(tabulate(match(ux, ux)))]
#   }
#   
#   df %>%
#     group_by(!!by_sym) %>%
#     mutate(!!col_sym := if_else(is.na(!!col_sym), mode_fn(!!col_sym), !!col_sym)) %>%
#     ungroup()
# }

### GET TRAITS  
# get appropriate data 
df <- traits %>% 
  dplyr::select(species,genus,family,order,
                `maximum_mass_(g)`,
         `mean__clutch_size`,
         `mean_number_of_clutches_per_year`,
         microhabitat,
         `maximum_lifespan_(year)`,
         neck_retraction_mechanism)

numeric_cols <- c("maximum_mass_(g)", "mean__clutch_size", 
                     "mean_number_of_clutches_per_year","maximum_lifespan_(year)")
df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)

all_cols <- c("maximum_mass_(g)", "mean__clutch_size", 
                  "mean_number_of_clutches_per_year","maximum_lifespan_(year)",
              "microhabitat",'neck_retraction_mechanism')
# check for NAs
check_names_na_cols <- function(x){
  names(which(colSums(is.na(x[x$species %in% species_tree_plot$tip.label,all_cols]))>0))
}
NA_cols_gen<-check_names_na_cols(df)
NA_cols_gen
# 3 cols with NAs for our species and columns of interest 
# impute by genus 
df_2_iter <- df
for (col in NA_cols_gen) {
  df_2_iter <- impute_numeric_by_genus(df_2_iter, col,by = 'genus')
}
head(df_2_iter,n=10)
head(df,n=10)

# repeat 
NA_cols_fam <- check_names_na_cols(df_2_iter)
NA_cols_fam
# 2 cols with NAs for our species and columns of interest 
# impute by family 
for (col in NA_cols_fam) {
  df_2_iter <- impute_numeric_by_genus(df_2_iter, col,by = 'family')
}
NA_cols_ord<-check_names_na_cols(df_2_iter)
NA_cols_ord
# no NAs left 

### NOW SUBSET FOR OUR SPECIES 
df_ours <- df_2_iter[df_2_iter$species%in%species_tree_plot$tip.label,]
check_names_na_cols(df_ours)
colnames(df_ours) <- c('Species','Genus','Family','Order',
                       'Mass',
                       'Clutchsize','Clutch_frequency',
                       'Microhabitat',
                       'Lifespan','Neck_Retraction')

df <- left_join(df_ours,n_genes[c('Species','unique_outliers')],by='Species')

# construct formula
# predictors <- c("Clutchsize")
# fmla <- as.formula(paste("unique_outliers ~", paste(predictors, collapse = " + ")))

df$Microhabitat <- as.factor(df$Microhabitat)
df$Neck_Retraction <- as.factor(df$Neck_Retraction)

# test brm
b.1 <- brm(unique_outliers ~ Mass + Clutchsize + Clutch_frequency + Lifespan + Microhabitat + Neck_Retraction + (1|gr(Species, cov = C)), 
           data   = df,
           family = gaussian(), 
           data2  = list(C = C),
           cores  = 4, 
           chains = 4,
           iter = 100000,
           thin = 10
           
)

saveRDS(b.1,
        file='traits-sidequest/results/mod1.rds')
plot(b.1)   
b.1

# Get posterior summary of all parameters
library(dplyr)
library(tibble)
library(ggplot2)

# Extract posterior summary
param_df <- posterior_summary(b.1) %>%
  as_tibble(rownames = "parameter") %>%
  dplyr::select(parameter, Estimate, lower = Q2.5, upper = Q97.5)

desired_order <- c(
#  "b_Intercept",                 # Intercept first
  "b_Lifespan",                  # Life history
  "b_Clutchsize",
  "b_Clutch_frequency",
  "b_Mass",                      # Mass
  "b_Neck_RetractionSideMnecked", # Neck retraction
  "b_MicrohabitatAquaticPMarine", # Microhabitat
  "b_MicrohabitatMarine",
  "b_MicrohabitatTerrestrial",
  "sd_Species__Intercept",       # Phylogenetic SD
  "sigma"                        # Residual SD
)

param_df <- param_df %>%
  filter(parameter %in% desired_order)

param_df <- param_df %>%
  mutate(parameter = factor(parameter, levels = rev(desired_order))) 

param_df <- param_df %>%
  mutate(group = case_when(
    grepl("^b_Intercept", parameter) ~ "Intercept",
    grepl("^b_Micro", parameter) ~ "Microhabitat",
    grepl("^b_Neck", parameter) ~ "Neck-retraction mechanism",
    grepl("^b_Mass", parameter) ~ "Morphology",
    grepl("^b_Clu|Lif", parameter) ~ "Life-history",
    grepl("^(sd_)", parameter) ~ "Random effects",
    grepl("^(sigma)", parameter) ~ "Random effects",
    TRUE ~ "Other"
  ))
param_df <- param_df %>%
  mutate(group = factor(group, levels = c('Life-history',
                                              'Morphology',
                                          'Neck-retraction mechanism',
                                              'Microhabitat',
                                              'Random effects')))

labels <- c('Lifespan','Clutchsize','Clutch frequency',
            'Mass','Side-necked',
            'Aquatic + Marine','Marine','Terrestrial',
            'Phylogenetic','Residual')

ggplot(param_df, aes(x = parameter, y = Estimate,
                     ymin = lower, ymax = upper,
                     color = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +  # line at 0
  geom_pointrange(size = 0.8) +
  coord_flip() +  # horizontal layout
  theme_minimal(base_size = 14) +
  scale_x_discrete(labels=rev(labels)) +
  labs(x = "Parameter",
       y = "Posterior median ± 95% CI",
       color = NULL,
       title = "Posterior Estimates of Model Parameters")+
  theme(legend.position = c(0.8,0.8),
        legend.box.background = element_rect(fill='white',
                                             linewidth = 0.5))

