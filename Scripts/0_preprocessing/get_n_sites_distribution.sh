# Get number of sites for all genes to look at distribution and decide cut-off

# For each gene
cat ../../Data/0_preprocessing/list_enough_taxa.txt | while read gene
do

	# get the number of sites as couting the characters of one row that is not a header but a sequence
	n_sites=($(seqkit seq -w 0 ../../Data/0_preprocessing/BUSCO_alignments_trimmed/${gene}.fa | head -2 | tail -1 |  wc -m))

	# write number of sites in file (can then be plotted in R)
    echo $n_sites >> ../../Data/0_preprocessing/all_n_sites.txt
done