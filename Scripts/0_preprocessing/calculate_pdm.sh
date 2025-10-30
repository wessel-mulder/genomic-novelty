# Run R script to calculate pairwise distance matrices on trimmed multiple sequence alignments

# Create output directory
mkdir ../../Data/0_preprocessing/dist_matrices

# For each gene in list of final genes (after all filtering steps)
cat ../../Data/0_preprocessing/list_final.txt | while read gene
do

	# copy prepared R script
	cp calculate_pdm.R calculate_pdm_gene.R
	# replace initial gene ID with current gene ID in duplicated script
	sed -i '' "s/1000at32523/${gene}/g" calculate_pdm_gene.R
	# run R script to calculate pairwise distance matrix
	Rscript calculate_pdm_gene.R

done