# Date Created: 12 July 2018
# Date Last Modified: 13 July 2018
# Execution: python differenceInFiles.py geneIdListFilename
# geneIdListFilename is a string representing the name of the file containing ids for which to search in the
# non-overlapping genes file. 
# Description: Determines whether ids listed in the input file are in the non-overlapping genes file.

#!/usr/bin/env python
import sys

idsNotFoundFilename = sys.argv[0]
idsNotFoundFile = open(idsNotFoundFilename, 'r')

idsNotFound = []
for geneId in idsNotFoundFile:
	idsNotFound.append(geneId)

idsNotFoundFile.close()

nonOverlappingGenesFilename = "genesNotInGTEx.txt"
nonOverlappingGenesFile = open(nonOverlappingGenesFilename, 'r')

nonOverlappingGenes = []
for gene in nonOverlappingGenesFile:
	nonOverlappingGenes.append(gene)

nonOverlappingGenesFile.close()

errorIds = []
for ensgId in idsNotFound:
	if ensgId not in nonOverlappingGenes:
		print ensgId, "wasn't found in the non-overlapping genes file."
		errorIds.append(ensgId)

if len(errorIds) == 0:
	print "All ids are listed in the non-overlapping genes list."