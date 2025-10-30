# Check that all genes have at least 9 taxa

# Create output directory
mkdir ../../Data/0_preprocessing/BUSCO_alignments_enough_taxa/

# For each gene
cat ../../Data/0_preprocessing/list_all.txt | while read gene
do
	# Get number of taxa in alignment as number of IDs (number of ">")
	n_taxa=($(grep ">" ../../Data/0_preprocessing/BUSCO_alignments_no_outgroups/${gene}.fa | wc -l))
	# if number of taxa > 8:
	if [ ${n_taxa} -gt 8 ]; then
		# keep gene
		cp ../../Data/0_preprocessing/BUSCO_alignments_no_outgroups/${gene}.fa ../../Data/0_preprocessing/BUSCO_alignments_enough_taxa/${gene}.fa
	fi
done
