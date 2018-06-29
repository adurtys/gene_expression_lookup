# Date Created: 27 June 2018
# Date Last Modified: 29 June 2018
# Execution: python GTF_processing_tests.py
# This program proceses GENCODE Comprehensive Gene Annotation GTF File for the start and end locations of protein-coding genes. 
# The output is a file containing locations of protein-coding genes, and maintains both GeneID (ENSG#) and name. 

#!/usr/bin/env python

# create a new file with four columns corresponding to GeneID, Gene Name, start location and end location
# read in the GTF file
inFilename = "gencode.v19.annotation.gtf"
outFilename = "gene_annotations.txt"

inFile = open(inFilename, 'r')
outFile = open(outFilename, 'w')

tab = "\t"
newline = "\n"

# skip first five lines
for i in range(5):
        inFile.readline()

for line in inFile:
   	line = line.rstrip('\r\n')
	
	# split line on tab, return list of columns
	columns = line.split('\t')
    
    chromosome = columns[0]
    featureType = columns[2]    
	geneStart = columns[3]
	geneEnd = columns[4]
	
	# process otherInfo column by getting rid of label and quotes
	otherInfo = columns[8]
	otherInfo = otherInfo.split(';')

	geneId = otherInfo[0]
	geneId = geneId.strip('gene_id ')
	geneId = geneId.strip('"')

	geneName = otherInfo[4]
	geneName = geneName.strip('gene_name ')
	geneName = geneName.strip('"')

	geneType = otherInfo[2]
	geneType = geneType.strip('gene_type ')
	geneType = geneType.strip('"')

	output = chromosome + tab + geneId + tab + geneName + tab + featureType + tab + geneStart + tab + geneEnd + newline
	
	# create list of protein-coding genes
	if geneType == "protein_coding":
		outFile.write(output)

outFile.close()
inFile.close()