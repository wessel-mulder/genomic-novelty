# Create output directory
mkdir ../../Data/0_preprocessing/codeml_out/

# for each gene in list of final genes (after all filtering steps)
cat ../../Data/0_preprocessing/list_final.txt | while read gene
do
	# create output folder
	mkdir ../../Data/0_preprocessing/codeml_out/${gene}
	# copy control file
	cp dnds.ctl ../../Data/0_preprocessing/codeml_out/${gene}/dnds_${gene}.ctl
	# replace initial gene ID with ID of current gene
	sed -i '' -e "s/66at32523/${gene}/g" ../../Data/0_preprocessing/codeml_out/${gene}/dnds_${gene}.ctl
	cd ../../Data/0_preprocessing/codeml_out/${gene}/
	# run codeML as defined in control file, keep output in logfile
	../../../../../dNdS_analysis/paml4.8/bin/codeml dnds_${gene}.ctl | tee logfile.txt
	cd ../../../../Scripts/0_preprocessing/
done