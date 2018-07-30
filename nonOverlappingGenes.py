# Date Created: 11 July 2018
# Date Last Modified: 26 July 2018
# Execution: python nonOverlappingGenes.py tstatFilename geneAnnotationsFilename
# argv1: filename for file containing tissue expression t-statistics
# argv2: filename for gene annotations file
# Description: finds genes listed in gene annotations file but which do not have tissue expression t-statistics in GTEx.
# 	Outputs these genes in "genesWithoutTstat.txt"
# Run Time: 1 sec

#!/usr/bin/env python
import sys

# check to make sure file was run with correct number of arguments
print "Arguments passed into nonOverlappingGenes.py:", sys.argv[1:]
if len(sys.argv) != 3:
	print "ERROR (nonOverlappingGenes.py line 15): Incorrect number of command-line arguments!"

tstatFilename = sys.argv[1]
annotationsFilename = sys.argv[2]

genesWithoutTstatsDict = {}

def main(tstatisticsFilename, geneAnnotationsFilename):
	# read in t-statistics file
	tstatFile = open(tstatisticsFilename, 'r')

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
	annotationsFile = open(geneAnnotationsFilename, 'r')

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

	output = ""
	for item in idsInGeneAnnotationsFile:
		# only add new ids to outFile once
		if (item not in idsInTstatFile) and (item not in genesWithoutTstatsDict):
			genesWithoutTstatsDict[item] = 0
			output += item + newline

	outFile.write(output)
	outFile.close()

	print "There are", len(genesWithoutTstatsDict), "genes in the gene annotations file but that do not have tissue expression t-statistics."

# create function to access genesWithoutTstat
def getGenesWithoutTstat(tstatisticsFilename, geneAnnotationsFilename):
	main(tstatisticsFilename, geneAnnotationsFilename)
	return genesWithoutTstatsDict.keys()

if __name__ == "__main__":
	main(tstatFilename, annotationsFilename)