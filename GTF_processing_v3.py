# Date Created: 27 June 2018
# Date Last Modified: 28 June 2018
# Execution: 
# This program proceses GENCODE Comprehensive Gene Annotation GTF File for the start and end locations of protein-coding genes. 
# The output is a file containing locations of protein-coding genes, and maintains both GeneID (ENSG#) and name. 

#!/usr/bin/env python

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
        
	geneStart = int(columns[3])
	geneEnd = int(columns[4])
	
	# process other info for gene_id and gene_name given that gene_type is protein-coding
	otherInfo = columns[8]
	otherInfo = otherInfo.split(';')

	geneId = otherInfo[0]
	geneName = otherInfo[4]
	
	# parse file for protein-coding genes
	geneType = otherInfo[2]
	# if geneType == "gene_type "protein_coding"":
	
	outFile.write(geneId)
	outFile.write(tab)
	outFile.write(geneName)
	outFile.write(tab)
	outFile.write(geneStart)
	outFile.write(tab)
	outFile.write(geneEnd)
	outFile.write(newline)

# for each line: 
	# determine whether gene already has start and end specified in new file
	# check both GeneID and Gene Name
		# if no:
			# write new line in new file for gene
		# if yes:
			# check start position
			# if GTF_start < newFile_start, then newFile_start == GTF_start
			
			# check end position
			# if GTF_end > newFile_end, then newFile_end == GTF_end

outFile.close()
inFile.close()
