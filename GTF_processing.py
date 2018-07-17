# Date Created: 27 June 2018
# Date Last Modified: 16 July 2018
# Execution: python GTF_processing.py gencodeFilename tstatFilename
# Description: This program proceses GENCODE Comprehensive Gene Annotation GTF File for the start and end locations of protein-coding genes. 
# 	The input is "gencode.v19.annotation.gtf", and the output is "gene_annotations.txt", a file with five columns (tab separated).
# 	The columns are chromosome name, GeneID, gene name, gene start location, and gene end location, respectively, for all protein-coding genes in the GENCODE file. 
# Run Time: 15 sec

#!/usr/bin/env python
import sys

# check to make sure file was run with correct number of arguments
if len(sys.argv) != 3:
	print "ERROR: Incorrect number of command-line arguments!"

# read in the GTEx file
tstatFilename = sys.argv[2]
tstatFile = open(tstatFilename, 'r')

# skip header line
tstatFile.readline()

# make list of all geneIds in tstat file
idsInGTEx = {}
for line in tstatFile:
	line = line.rstrip('\r\n')

	columns = line.split('\t')

	ensgId = columns[0]
	idsInGTEx[ensgId] = 0

tstatFile.close()

# read in the GTF file
inFilename = sys.argv[1]
inFile = open(inFilename, 'r')

# create dictionary containing geneId as key, and a list of relevant info as the value
genes = {}

# skip first five lines
for i in range(5):
	inFile.readline()

for line in inFile:
	line = line.rstrip('\r\n')
	
	# split line on tab, creating list of possible columns
	columns = line.split('\t')

	chromosome = columns[0]
	geneStart = int(columns[3])
	geneEnd = int(columns[4])

	# process otherInfo column by getting rid of key label and quotes
	otherInfo = columns[8]
	otherInfo = otherInfo.split(';')

	geneId = otherInfo[0]
	geneId = geneId.strip('gene_id ')
	geneId = geneId.strip('"')

	geneId = geneId.split('.')
	geneId = geneId[0]

	name = otherInfo[4]
	name = name.strip('gene_name ')
	name = name.strip('"')

	geneType = otherInfo[2]
	geneType = geneType.strip('gene_type ')
	geneType = geneType.strip('"')

	# create list to store relevant info
	geneInfo = []

	# check if the current geneId is in the gene dictionary
	if geneId in genes:
		# compare start locations
		if geneStart < genes[geneId][2]:
			# update start location
			genes[geneId][2] = geneStart
		# compare end locations
		if geneEnd > genes[geneId][3]:
			# update end location
			genes[geneId][3] = geneEnd

	# parse for genes that are either protein-coding or contained in the GTEx file
	elif (geneType == "protein_coding") or (geneId in idsInGTEx):
		# populate geneInfo with relevant info
		geneInfo.append(name)
		geneInfo.append(chromosome)
		geneInfo.append(geneStart)
		geneInfo.append(geneEnd)
		geneInfo.append(geneType)

		genes[geneId] = geneInfo

inFile.close()

numGenes = len(genes)
print "Writing", numGenes, "genes to the annotations file. These genes are either protein-coding or have tissue expression data in the GTEx file."

# create a new file for start and end positions of only protein-coding genes
outFilename = "gene_annotations.txt"
outFile = open(outFilename, 'w')

tab = "\t"
newline = "\n"

for gene in genes:
	# convert start and end positions back to strings
	genes[gene][2] = str(genes[gene][2])
	genes[gene][3] = str(genes[gene][3])

	output = gene + tab
	for i in range(len(genes[gene])):
		if i == len(genes[gene]) - 1:
			output += genes[gene][i] + newline
		else:
			output += genes[gene][i] + tab

	outFile.write(output)

outFile.close()