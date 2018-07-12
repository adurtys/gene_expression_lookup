#!/usr/bin/env python

idsNotFoundFilename = "idsNotFound_fullTest.txt"
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
	print "All ids that weren't in the GTEx file are listed in the non-overlapping genes list."