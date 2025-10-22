# GETTING STARTED ---------------------------------------------------------
# packages

library(devtools)
devtools::install_github("stan-dev/cmdstanr")
devtools::install_github("rmcelreath/rethinking")

pacman::p_load(tidyverse,brms,MCMCglmm,parallel,
               knitr,geiger,ape,phytools,caper,bayesplot,
               ggtree,gridExtra,MCMCpack,dplyr,tidyr,vegan,
               rethinking,cmdstanr)



theme_set(theme_tidybayes() + panel_border())

# tutorial functions file
source("MR-PMM_tutorial_functions.R")

### 

# LOADING DATA  -----------------------------------------------------------
# TURTLE METADATA
meta_turtles <- read_tsv("Data/2_comparison_outliers/metadata_habitat.tsv")
meta_turtles <- meta_turtles %>% filter(Microhabitat != "Outgroup")
meta_turtles$Habitat_factor <- factor(meta_turtles$Microhabitat, 
                                      levels=c("Marine", "Aquatic_Marine", 
                                               "Aquatic", "Terrestrial"))

### GET NUMBER OF DNDS OUTLIERS
outliers_dnds_species <- read.csv("Results/1_mahalanobis_outliers/outliers_species/dnds_chi95.csv")[,-1]
colnames(outliers_dnds_species) <- c("ID", "gene")

# make matrix of genes
mat <- outliers_dnds_species %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = gene,
    values_from = present,
    values_fill = 0
  ) %>%
  column_to_rownames("ID") %>%
  as.matrix()
head(mat,n=16)
head(outliers_dnds_species,n=16)

# jaccard distance
dist_jaccard <- as.matrix(vegdist(mat, method = "jaccard"))
dist_jaccard[dist_jaccard == 0] <- NA
uniqueness <- rowMeans(dist_jaccard,na.rm=T)

# merge with meta file 
meta_turtles <- merge(meta_turtles,data.frame('ID' = names(uniqueness),
                                        'Jaccard' = uniqueness),
                by = 'ID')

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

df <- left_join(df_ours,meta_turtles[c('Species','Jaccard')],by='Species')

# calculate jaccard index 

# MULTI-RESPONSE  ---------------------------------------------------------
df2 <- df
df2 <- df2 %>%
  mutate(
    Microhabitat   = gsub("\\+","_", Microhabitat),   # replace + with _
    Microhabitat   = gsub(" ", "", Microhabitat),     # remove spaces
    Neck_Retraction = gsub("-", "_", Neck_Retraction) # replace dash with underscore
  )

# Microhabitat into binary dummies
microhab_dummies <- model.matrix(~ Microhabitat - 1, data = df2)

# Neck retraction into binary dummies
neck_dummies <- model.matrix(~ Neck_Retraction - 1, data = df2)

# Bind them back
df2 <- cbind(df2, microhab_dummies, neck_dummies)

#standardize 
df2$Jaccard <- standardize(df2$Jaccard)
df2$Mass <- standardize(df2$Mass)
df2$Clutchsize <- standardize(df2$Clutchsize)
df2$Clutch_frequency <- standardize(df2$Clutch_frequency)
df2$Lifespan <- standardize(df2$Lifespan)

head(df2)
# add observation 
df2$obs <- seq_len(nrow(df2))

# Gaussian responses
f1 <- brmsformula(Jaccard ~ (1|u|gr(Species, cov = C)) + (1|e|obs), sigma = 0.1) + gaussian()
f2 <- brmsformula(Mass            ~ (1|u|gr(Species, cov = C)) + (1|e|obs), sigma = 0.1) + gaussian()
f3 <- brmsformula(Clutchsize      ~ (1|u|gr(Species, cov = C)) + (1|e|obs), sigma = 0.1) + gaussian()
f4 <- brmsformula(Clutch_frequency~ (1|u|gr(Species, cov = C)) + (1|e|obs), sigma = 0.1) + gaussian()
f5 <- brmsformula(Lifespan        ~ (1|u|gr(Species, cov = C)) + (1|e|obs), sigma = 0.1) + gaussian()
# Bernoulli responses
# Microhabitat dummies
f6  <- brmsformula(MicrohabitatAquatic        ~ (1|u|gr(Species, cov = C)) + (1|e|obs)) + bernoulli()
f7 <- brmsformula(MicrohabitatAquatic_Marine ~ (1|u|gr(Species, cov = C)) + (1|e|obs)) + bernoulli()
f8 <- brmsformula(MicrohabitatMarine         ~ (1|u|gr(Species, cov = C)) + (1|e|obs)) + bernoulli()
f9  <- brmsformula(MicrohabitatTerrestrial    ~ (1|u|gr(Species, cov = C)) + (1|e|obs)) + bernoulli()
# Neck retraction dummies
f10 <- brmsformula(Neck_RetractionHidden_necked ~ (1|u|gr(Species, cov = C)) + (1|e|obs)) + bernoulli()
f11   <- brmsformula(Neck_RetractionSide_necked   ~ (1|u|gr(Species, cov = C)) + (1|e|obs)) + bernoulli()
# Combine all into one mvbrmsformula
b.2.f <- f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11

# test brm
b.2 <- brm(b.2.f + set_rescor(F),
           data   = df2,
           data2  = list(C = C),
           cores  = 4, 
           chains = 4,
           iter = 10000,
           thin = 5
)

summary(b.2)
colnames(as_tibble(b.2))

param_df <- posterior_summary(b.2) %>%
  as_tibble(rownames = "parameter") %>%
  dplyr::select(parameter, Estimate, lower = Q2.5, upper = Q97.5)

# keep correlation parameters with "uniqueoutliers"
sel <- grepl("^(cor|rescor).*Jaccard", param_df$parameter)
outliers_corr <- param_df[sel, ]

outliers_corrs <- outliers_corr %>%
  # make nice groups 
  mutate(Trait = case_when( 
    grepl("Mass|NeckRetraction", parameter) ~ "Morphology",
    grepl("Clutchsize|Clutchfrequency|Lifespan", parameter) ~ "Life-history",
    grepl("Microhabitat", parameter) ~ "Habitat",
    TRUE ~ "Other"  # fallback if needed
  )) %>%
  # make phylo / resid correlations
  mutate(
    Correlation = case_when(
      str_starts(parameter, "cor_obs__") ~ "Residual",
      str_starts(parameter, "cor_Species__") ~ "Phylogenetic",
      TRUE ~ NA_character_  # just in case there are other parameters
    )
  ) %>%
  #  make nice names 
  mutate(unique_param = paste0(
    "Jaccard + ",
    str_extract(parameter, "(?<=Intercept__)[A-Za-z0-9]+")
  )) %>%
  mutate(unique_param = str_remove_all(unique_param, "Microhabitat|NeckRetraction")) %>%
  mutate(unique_param = unique_param %>%
           # Add dashes for known multi-word cases
           str_replace_all("Clutchfrequency", "Clutch-frequency") %>%
           str_replace_all("AquaticMarine", "Aquatic-Marine") %>%
           str_replace_all("Hiddennecked", "Hidden-necked") %>%
           str_replace_all("Sidenecked", "Side-necked")
  ) 


# Reorder parameter by col

outliers_corrs <- outliers_corrs %>%
  mutate(Trait = factor(Trait, 
                      levels = c('Morphology','Life-history','Habitat'))) %>%
  mutate(Correlation = factor(Correlation, 
                        levels = c('Phylogenetic','Residual'))) %>%
  mutate(unique_param = factor(unique_param,
                               levels = unique(unique_param[order(Trait)])))

# and plot
p<-ggplot(outliers_corrs, aes(x = fct_rev(unique_param), y = Estimate,
                     ymin = lower, ymax = upper,
                     color = Trait,linetype = Correlation,
                     shape = Correlation)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +  # line at 0
  geom_pointrange(size = 0.5, 
                  position = position_dodge(width = -0.6)) +  # <-- dodge corr_type
  coord_flip() +  # horizontal layout
  theme_minimal(base_size = 16) +
  labs(x = NULL,
       y = "Posterior median ± 95% CI",
       title = "Posterior estimates of trait correlations")+
  scale_color_manual(values = c(
    "Morphology" = "#d95f02",       # greenish
    "Life-history" = "#7570b3",     # orange
    "Habitat" = "#1b9e77"           # purple
  )) +
  theme(#legend.position = c(1,1),
        legend.box.background = element_rect(fill='white',
                                             linewidth = 0.1))
print(p)

pdf(file = 'traits-sidequest/figs/posterior-estimates-traits.pdf',
    width = 10,
    height = 10)
print(p)
dev.off()

write.csv(df2,
          file = 'traits-sidequest/data/subset-traits.csv')

saveRDS(b.2,
        file = 'traits-sidequest/results/jacard-mod.rds')


