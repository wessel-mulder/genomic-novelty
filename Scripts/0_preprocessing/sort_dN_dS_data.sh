# Script to copy all relevant codeml output into one folder, rename by gene

# Create output directory
mkdir ../../Data/0_preprocessing/dN_dS_data

# For each gene
cat ../../Data/0_preprocessing/list_final.txt | while read gene
do

	# copy 2ML.dN and 2ML.dS into new folder, rename according to gene ID
	cp ../../Data/0_preprocessing/codeml_out/${gene}/2ML.dN ../../Data/0_preprocessing/dN_dS_data/${gene}_dN.dist
	cp ../../Data/0_preprocessing/codeml_out/${gene}/2ML.dS ../../Data/0_preprocessing/dN_dS_data/${gene}_dS.dist
done