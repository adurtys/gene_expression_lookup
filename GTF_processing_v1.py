# Date Created: 27 June 2018
# Date Last Modified: 27 June 2018
# Execution: 
# This program proceses GENCODE Comprehensive Gene Annotation GTF File for the start and end locations of protein-coding genes. 
# The output is a file containing locations of protein-coding genes, and maintains both GeneID (ENSG#) and name. 

#!/usr/bin/env python

# read in the GTF file (tab-separated, header = true)
inFilename = "gencode.v19.annotation.gtf"
outFilename = "gene_annotations.txt"

inFile = open(inFilename, 'r')
outFile = open(outFilename, 'w')

inFile.readline() # skip header line

for line in inFile:
	line = line.rstrip('\r\n')
	columns = line.split() # split line on whitespace (tab), returns list of columns

	geneStart = columns[3] + "\t"
	geneEnd = columns[4] + "\t"
	otherInfo = columns[8] + "\n"

	outFile.write(geneStart)
	outFile.write(geneEnd)
	outFile.write(otherInfo)
	
	# process other info for gene_id and gene_name given that gene_type is protein-coding
	# otherInfo = otherInfo.split(';')



# parse file for protein-coding genes

# write new file with four columns (GeneID, Gene Name, Start, End)

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