
input <- 'Results/0_preprocessing'
output <- 'Results/1_mahalanobis_outliers/'
scripts <- 'Scripts/1_mahalanobis_outliers/'

# LOAD DATA ----------------------------------------------------------------
### READ PWD
list <- list.files(file.path(input,'dist_matrices'),
                   pattern = '*.dist',
                   full.names = T)

matrices <- list()
for(i in list){
  base <- basename(i)
  name <- strsplit(base,'.dist')[[1]]
  dist <- read.table(i)
  matrices[[name]] <- as.matrix(dist)
}

matrices_nonzero <- lapply(matrices, function(mat) {
  mat[mat == 0] <- NA
  return(mat)
})

### READ dNdS
# Define the folder containing the data
data_folder <- file.path(input,"dN_dS_final_data")

# List all _dN.dist and _dS.dist files
dN_files <- list.files(data_folder, pattern = "_dN\\.dist$", full.names = TRUE)
dS_files <- list.files(data_folder, pattern = "_dS\\.dist$", full.names = TRUE)


# Ensure the files are paired and match
if (!all(gsub("_dN\\.dist$", "", basename(dN_files)) == gsub("_dS\\.dist$", "", basename(dS_files)))) {
  stop("Mismatch between dN and dS files!")
}

# Initialize empty lists to store matrices
dN_matrices <- list()
dS_matrices <- list()

# Loop through all files and read data
for (i in seq_along(dN_files)) {
  # Check if the file is empty
  if (file.info(dN_files[i])$size == 0 || file.info(dS_files[i])$size == 0) {
    warning(paste("Skipping empty file pair:", dN_files[i], "and", dS_files[i]))
    next
  }
  
  matrix_name <- gsub("_dN\\.dist$", "", basename(dN_files[i]))
  
  
  
  # Read the number of species from the first line
  n_species <- as.integer(readLines(dN_files[i], n = 1))
  
  # Load dN and dS matrices
  dN <- read.table(dN_files[i], row.names = 1, skip = 1, fill = TRUE, 
                   col.names = paste0("V", 1:(n_species + 1)))
  dS <- read.table(dS_files[i], row.names = 1, skip = 1, fill = TRUE, 
                   col.names = paste0("V", 1:(n_species + 1)))
  
  # Adjust column names to match species names
  colnames(dN) <- rownames(dN)
  colnames(dS) <- rownames(dS)
  
  # Store matrices in lists
  dN_matrices[[matrix_name]] <- as.matrix(dN)
  dS_matrices[[matrix_name]] <- as.matrix(dS)
}

# Check the data
summary(dN_matrices)
summary(dS_matrices)

# Initialize an empty list to store the dN/dS results
dN_matrices_nonzero <- list()
dS_matrices_nonzero <- list()

dN_nonzero <- lapply(dN_matrices, function(mat) {
  mat[mat == 0] <- NA
  return(mat)
})

dS_nonzero <- lapply(dS_matrices, function(mat) {
  mat[mat == 0] <- NA
  return(mat)
})


dN_dS_matrices_nonzero <- list()
# Loop through the names of dN_matrices
for (name in names(dN_matrices)) {
  # Check if there is a corresponding dS matrix
  if (!name %in% names(dS_matrices)) {
    warning(paste("No matching dS matrix found for", name))
    next
  }
  print(name)
  
  
  
  # Perform element-wise division
  dN_dS <- dN_nonzero[[name]] / dS_nonzero[[name]]
  
  # Store the result in the list
  dN_dS_matrices_nonzero[[name]] <- dN_dS
}

### check if all goes well
dN_nonzero[['121468at32523']]
dS_nonzero[['121468at32523']]
dN_dS_matrices_nonzero[['121468at32523']]

# EXTRACT OUTLIERS -----------------------------------------------
# Required Libraries
library(MASS)  # For Mahalanobis distance
library(stats) # For Chi-squared distribution

### LOAD MATRIX OF INTEREST 

source(file.path(scripts,'helperFunctions.R'))

#dN_dS_matrices[[1]]
dNdS_checked <- matrixCheck(dN_dS_matrices_nonzero)
dNdS_checked[['121468at32523']]

dN_checked <- matrixCheck(dN_nonzero)
dN_checked[['121468at32523']]

dS_checked <- matrixCheck(dS_nonzero)
dS_checked[['121468at32523']]

pure_checked <- matrixCheck(matrices_nonzero)
pure_checked[['121468at32523']]

#impute from average across matrix
source(file.path(scripts,'impMeanMatrixAverage.R'))
dNdS_impute <- impMean_matrices(dNdS_checked)
dN_impute <- impMean_matrices(dN_checked)
dS_impute <- impMean_matrices(dS_checked)
pure_impute <- impMean_matrices(pure_checked)

dNdS_impute[['121468at32523']]
dN_impute[['121468at32523']]
dS_impute[['121468at32523']]
pure_impute[['121468at32523']]

# make combined dataframe
source(file.path(scripts,'stitchMatricesToDataFrame.R'))
dNdS_2d <- stitchMatricesToDataFrame(dNdS_impute)
dN_2d <- stitchMatricesToDataFrame(dN_impute)
dS_2d <- stitchMatricesToDataFrame(dS_impute)
pure_2d <- stitchMatricesToDataFrame(pure_impute)

dNdS_log <- log(dNdS_2d)
dN_log <- log(dN_2d)
dS_log <- log(dS_2d)
pure_log <- log(pure_2d)

#detect outliers
source(file.path(scripts,'detectOutliers.R'))

# Detect outliers and extract quantiles and chi-square for meangene
dNdS_outliers <- detect_outliers_and_extract_quantiles(dNdS_impute, dNdS_log, quantile_threshold = 0.95)
dN_outliers <- detect_outliers_and_extract_quantiles(dN_impute, dN_log, quantile_threshold = 0.95)
dS_outliers <- detect_outliers_and_extract_quantiles(dS_impute, dS_log, quantile_threshold = 0.95)
pure_outliers <- detect_outliers_and_extract_quantiles(pure_impute, pure_log, quantile_threshold = 0.95)

# Your results
outlier_obs <- list(
  dfs = list(dNdS_outliers, dN_outliers, dS_outliers, pure_outliers),
  names = c('dnds', 'dn', 'ds', 'raw')  # no need for list()
)

df <- dNdS_outliers
# Use Map to loop through both dfs and names simultaneously
Map(function(df, name) {
  if (!is.null(df$indices)) {
    outlier_indices <- data.frame(
      busco_id = df$indices,
      pval = df$pchisq_all[df$indices]
    )
    write.csv(outlier_indices, file.path(output,'outliers_genes',paste0(name, '_q95.csv')), row.names = FALSE)
    all_indices <- data.frame(
      busco_id = names(df$distances_all),
      distance = df$distances_all,
      pval = df$pchisq_all[names(df$distances_all)]
    )
    write.csv(all_indices,file.path(output,'distances_allgenes',paste0(name,'.csv')), row.names = FALSE)
  }
}, outlier_obs$dfs, outlier_obs$names)


# outlier species - just outlier genes ------------------------------------------------
# keep entire thing


outlierSpeciesTranspose <- function(matrices_list) {
  # Extract the gene names and the list of matrices
  gene_names <- matrices_list[[1]]
  matrices <- matrices_list[[2]]
  
  # Initialize an empty list to store the renamed matrices
  matrices_renamed <- list()
  
  # Loop through the matrices and rename rows based on gene names
  for (i in 1:length(matrices)) {
    # Get the current matrix
    mat <- matrices[[i]]
    
    # Assign the new row names by appending the gene name
    rownames(mat) <- paste0(rownames(mat), '_', names(matrices)[i])
    
    # Add the renamed matrix to the list
    matrices_renamed[[i]] <- mat
  }
  
  # Combine the matrices by rows
  df <- do.call(rbind, matrices_renamed)
  
  return(df)
}

# transpose matrices 
dNdS_df <- outlierSpeciesTranspose(dNdS_outliers)
dN_df <- outlierSpeciesTranspose(dN_outliers)
dS_df <- outlierSpeciesTranspose(dS_outliers)
pure_df <- outlierSpeciesTranspose(pure_outliers)

# extract significant outliers
dNdS_outliers_species <- detect_outliers_and_extract_chisq(dNdS_impute, dNdS_df, conf = 0.95)
dN_outliers_species <- detect_outliers_and_extract_chisq(dN_impute, dN_df, conf = 0.95)
dS_outliers_species <- detect_outliers_and_extract_chisq(dS_impute, dS_df, conf = 0.95)
pure_outliers_species <- detect_outliers_and_extract_chisq(pure_impute, pure_df, conf = 0.95)

text <- dNdS_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(species = x[1], busco_id = x[2])))
write.csv(df,file.path(output,'outliers_species','dnds_chi95.csv'))

text <- dN_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(species = x[1], busco_id = x[2])))
write.csv(df,file.path(output,'outliers_species','dn_chi95.csv'))

text <- dS_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(species = x[1], busco_id = x[2])))
write.csv(df,file.path(output,'outliers_species','ds_chi95.csv'))

text <- pure_outliers_species$indices
split <- strsplit(text,'_')
df <- do.call(rbind, lapply(split, function(x) data.frame(species = x[1], busco_id = x[2])))
write.csv(df,file.path(output,'outliers_species','raw_chi95.csv'))

