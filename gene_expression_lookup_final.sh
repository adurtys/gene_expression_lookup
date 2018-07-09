# Date Created: 9 July 2018
# Date Last Modified: 9 July 2018
# TODO: Execution: python gene_expression_lookup_final.sh
# TODO: Description:

#!/usr/bein/env bash

if [ -f ./gene_annotations.txt ] && [ -f ./normalizedGTEx.tstat.txt ]
then
	echo "Both the GENCODE GTF file and the t-statistics file have been processed."
	echo "Now starting snp lookup."

	# check number of command line arguments
	if [ $# != 0 ]
	then
		echo "something is wrong..."
	else
		echo "command line contains $# arguments."
	fi
	
	# ./expression_lookup_final.py

elif [ -f ./gene_annotations.txt]
then
	echo "GENCODE GTF file has been processed."
	echo "Now processing t-statistics file."

	./tstat_normalization_final.py

	echo "Now starting snp lookup."
else [ -f ./normalizedGTEx.tstat.txt]
	echo "GTEx t-statistics file has been processed."
	echo "Now processing GENCODE GTF file"
	
	./GTF_processing_final.py

	echo "Now starting snp lookup."

	# ./expression_lookup_final.py

fi
