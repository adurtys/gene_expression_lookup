# Date Created: 11 July 2018
# Date Last Modified: 16 July 2018
# Execution: python nonOverlappingGenes.py tstatFilename geneAnnotationsFilename
# Description: Finds genes listed in the gene_annotations file (containing only the protein-coding genes from
# 	the GENCODE GTF file) but that do not have tissue expression data in the GTEx file, and vice versa. Prints out these 
# 	non-overlapping genes into separate files, namely "genesNotInGencode.txt" and "genesNotInGTEx.txt", respectively.
# Run Time: 10 sec

#!/usr/bin/env python
import sys

# check to make sure file was run with correct number of arguments
if len(sys.argv) != 3:
	print "ERROR: Incorrect number of command-line arguments!"

# read in GTEx t-statistics file
tstatFilename = sys.argv[1]
tstatFile = open(tstatFilename, 'r')

# skip header line
headerLine = tstatFile.readline()

# parse GTEx file for GeneIds
idsInTstatFile = []
for line in tstatFile:
	line = line.rstrip('\r\n')
	tissues = line.split('\t')

	# add ENSGID, found in first column, to the list
	idsInTstatFile.append(tissues[0])

tstatFile.close()

# read in gene_annotations file
annotationsFilename = sys.argv[2]
annotationsFile = open(annotationsFilename, 'r')

# parse file for GeneIds
idsInGencodeFile = []
for line in annotationsFile:
   	line = line.rstrip('\r\n')
	
	# split line on tab, creating list of possible columns
	columns = line.split('\t')
	geneId = columns[0]

	# add ENSGID to the list
	idsInGencodeFile.append(geneId)

annotationsFile.close()

newline = "\n"

# create output file for genes in GTEx file and not in GENCODE file
tstatGenesFilename = "genesNotInAnnotationsFile.txt"
tstatGenesFile = open(tstatGenesFilename, 'w')

genesNotInGencode = []
tstatOutput = ""
for item in idsInTstatFile:
	if (item not in idsInGencodeFile) and (item not in genesNotInGencode):
		genesNotInGencode.append(item)
		tstatOutput += item + newline

tstatGenesFile.write(tstatOutput)
tstatGenesFile.close()

# create output file for genes in GENCODE file and not in GTEx file
gencodeGenesFilename = "genesNotInGTEx.txt"
gencodeGenesFile = open(gencodeGenesFilename, 'w')

genesNotInGTEx = []
gencodeOutput = ""
for item in idsInGencodeFile:
	# only add new ids to outFile once
	if (item not in idsInTstatFile) and (item not in genesNotInGTEx):
		genesNotInGTEx.append(item)
		gencodeOutput += item + newline

gencodeGenesFile.write(gencodeOutput)
gencodeGenesFile.close()

print "There are", len(genesNotInGencode), "genes in the GTEx file but not in the gene annotations file."
print "There are", len(genesNotInGTEx), "genes in the gene annotations file but not in the GTEx file."