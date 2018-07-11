# Date Created: 9 July 2018
# Date Last Modified: 11 July 2018
# Execution: python gene_expression_lookup_final.sh snpFile numGenesToSearch distanceFromSnp expressionThreshold outFile
# snpFile is the name of the file containing the snps for which to search
# numGenesToSearch is the number of genes closest to the snp for which to look up tissue expression
# distanceFromSnp is the maximum distance upstream and downstream from the snp for which to search for genes, in kbp
# if distanceFromSnp > 1000 kbp, distanceFromSnp = 1000kbp (max distance from gene that will be searched is 1 mbp)
# expressionThreshold is the percentage (in decimal) of top rank-ordered t-statistics that should be considered "highly
# expressed" for each tissue
# Description: Shell script for tissue expression lookup of genes near specified snps. First, ensures that both the 
# GENCODE GTF file and GTEx t-statistics file have been processed. Then analyzes expression of 

#!/usr/bin/env bash

if [ -f ./gene_annotations.txt ] && [ -f ./normalizedGTEx.tstat.txt ]
then
	echo "Both the GENCODE GTF file and the t-statistics file have been processed."
	echo "Now starting snp lookup."

	# check number of command-line arguments
	if [ $# = 0 ]
	then
		echo "ERROR: incorrect number of command-line arguments."
		echo "This program requires 4 command-line arguments in addition to the name of the command."
	else
		# store command-line arguments as variables
		snpFile="$1"
		numGenesToSearch="$2"
		distanceFromSnp="$3"
		expressionThreshold="$4"
		
		outFile=geneExpressionLookupResults_$snpFile

		echo "snpFile: $1"
		echo "numGenesToSearch: $2"
		echo "distanceFromSnp: $3"
		echo "expressionThreshold: $4"
		echo $outFile

		python ./expression_lookup_final.py $snpFile $numGenesToSearch $distanceFromSnp $expressionThreshold $outFile
	fi

elif [ -f ./gene_annotations.txt]
then
	echo "GENCODE GTF file has been processed."
	echo "Now processing t-statistics file."

	python ./tstat_normalization_final.py

	echo "Now starting snp lookup."
	
	# check number of command-line arguments
	if [ $# = 0 ]
	then
		echo "ERROR: incorrect number of command-line arguments."
		echo "This program requires 4 command-line arguments in addition to the name of the command."
	else
		# store command-line arguments as variables
		snpFile="$1"
		numGenesToSearch="$2"
		distanceFromSnp="$3"
		expressionThreshold="$4"

		echo "snpFile: $1"
		echo "numGenesToSearch: $2"
		echo "distanceFromSnp: $3"
		echo "expressionThreshold: $4"

		python ./expression_lookup_final.py $snpFile $numGenesToSearch $distanceFromSnp $expressionThreshold
	fi

else [ -f ./normalizedGTEx.tstat.txt]
	echo "GTEx t-statistics file has been processed."
	echo "Now processing GENCODE GTF file"
	
	python ./GTF_processing_final.py

	echo "Now starting snp lookup."
	
	# check number of command-line arguments
	if [ $# = 0 ]
	then
		echo "ERROR: incorrect number of command-line arguments."
		echo "This program requires 4 command-line arguments in addition to the name of the command."
	else
		# store command-line arguments as variables
		snpFile="$1"
		numGenesToSearch="$2"
		distanceFromSnp="$3"
		expressionThreshold="$4"

		echo "snpFile: $1"
		echo "numGenesToSearch: $2"
		echo "distanceFromSnp: $3"
		echo "expressionThreshold: $4"

		python ./expression_lookup_final.py $snpFile $numGenesToSearch $distanceFromSnp $expressionThreshold
	fi
fi
