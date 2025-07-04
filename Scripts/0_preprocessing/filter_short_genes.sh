# Check that all genes have at least 500 sites

# For each gene
cat ../../Data/0_preprocessing/list_enough_taxa.txt | while read gene
do

    # get the number of sites as couting the characters of one row that is not a header but a sequence
	n_sites=($(seqkit seq -w 0 ../../Data/0_preprocessing/BUSCO_alignments_trimmed/${gene}.fa | head -2 | tail -1 |  wc -m))

    # if number of sites > 499:
    if [ ${n_sites} -gt 499 ]; then
        # keep gene
        cp ../../Data/0_preprocessing/BUSCO_alignments_trimmed/${gene}.fa ../../Data/0_preprocessing/BUSCO_alignments_trimmed_final/${gene}.fa
    fi
done
