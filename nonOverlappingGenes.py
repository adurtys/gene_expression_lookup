# Date Created: 11 July 2018
# Date Last Modified: 26 July 2018
# Execution: python nonOverlappingGenes.py tstatFilename geneAnnotationsFilename
# argv1: filename for file containing tissue expression t-statistics
# argv2: filename for gene annotations file
# Description: finds genes listed in gene annotations file but which do not have tissue expression t-statistics in GTEx.
# 	Outputs these genes in "genesWithoutTstat.txt"
# Run Time: 10 sec

#!/usr/bin/env python
import sys

# check to make sure file was run with correct number of arguments
if len(sys.argv) != 3:
	print "ERROR (nonOverlappingGenes.py line 15): Incorrect number of command-line arguments!"

# read in t-statistics file
tstatFilename = sys.argv[1]
tstatFile = open(tstatFilename, 'r')

# skip header line
headerLine = tstatFile.readline()

# parse t-statistics file for gene IDs
idsInTstatFile = {}
for line in tstatFile:
	line = line.rstrip('\r\n')
	tissues = line.split('\t')
	geneId = tissues[0]

	# add gene ID to dictionary
	idsInTstatFile[geneId] = 0

tstatFile.close()

# read in gene_annotations file
annotationsFilename = sys.argv[2]
annotationsFile = open(annotationsFilename, 'r')

# parse gene annotations file for gene ids
idsInGeneAnnotationsFile = {}
for line in annotationsFile:
   	line = line.rstrip('\r\n')
	
	# split line on tab, creating list of possible columns
	columns = line.split('\t')
	geneId = columns[0]

	# add gene ID to dictionary
	idsInGeneAnnotationsFile[geneId] = 0

annotationsFile.close()

newline = "\n"

# create output file for genes in GENCODE file and not in GTEx file
outFilename = "genesWithoutTstat.txt"
outFile = open(outFilename, 'w')

genesWithoutTstat = {}
output = ""
for item in idsInGeneAnnotationsFile:
	# only add new ids to outFile once
	if (item not in idsInTstatFile) and (item not in genesWithoutTstat):
		genesWithoutTstat[item] = 0
		output += item + newline

outFile.write(output)
outFile.close()

print "There are", len(genesWithoutTstat), "genes in the gene annotations file but that do not have tissue expression t-statistics."