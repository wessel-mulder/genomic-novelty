# Script to get rid of internal stop codons with macse (replace with NNN)

# Create output directories
mkdir ../../Data/0_preprocessing/BUSCO_alignments_no_stop
mkdir ../../Data/0_preprocessing/no_stop_AA

# For each gene
cat ../../Data/0_preprocessing/list_enough_taxa.txt | while read gene
do
	
	# run macse to replace internal stop codons in given alignment with NNN
	java -jar macse_v2.07.jar -prog exportAlignment -align ../../Data/0_preprocessing/BUSCO_alignments_correct_n_codons/${gene}.fa -codonForInternalStop NNN -codonForFinalStop --- -codonForInternalFS --- -charForRemainingFS - -out_NT ../../Data/0_preprocessing/BUSCO_alignments_no_stop/${gene}.fa -out_AA ../../Data/0_preprocessing/no_stop_AA/${gene}.fa
done