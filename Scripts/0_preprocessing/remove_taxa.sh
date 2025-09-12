# Script to remove a list of taxa (here outgroups) from multiple sequence alignment

# Create output directory
mkdir ../../Data/0_preprocessing/BUSCO_alignments_no_outgroups/

# For each gene
cat ../../Data/0_preprocessing/list_all.txt | while read gene
do
	# Remove all outgroups listed in outgroups.txt
	# grep -v inverts: keep all IDs that do not match to outgroups.txt
	seqkit grep -v -f ../../Data/0_preprocessing/outgroups.txt ../../Data/0_preprocessing/Trimmed_BUSCO_alignments/${gene}.trim > ../../Data/0_preprocessing/BUSCO_alignments_no_outgroups/${gene}.fa
done
