# Run R script to trim bad codons (meaning many gaps or high heterozygosity)

# Create output directory
mkdir ../../Data/0_preprocessing/BUSCO_alignments_trimmed

# For each gene
cat ../../Data/0_preprocessing/list_enough_taxa.txt | while read gene
do

	# copy prepared R script
	cp trim_bad_codons.R trim_bad_codons_gene.R
	# replace initial gene ID with current gene ID in duplicated script
	sed -i '' "s/1000at32523/${gene}/g" trim_bad_codons_gene.R
	# run R script to remove gappy or heterozygous codons
	Rscript trim_bad_codons_gene.R

done
