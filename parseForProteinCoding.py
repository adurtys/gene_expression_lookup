# Date Created: 12 July 2018
# Date Last Modified: 13 July 2018
# Execution: python parseForProteinCoding.py geneIdList
# geneIdList is a string containing the name of the file containing the genes to look up
# Description: checks whether a list of ENSGIDs is protein-coding, according to GENCODE

#!/usr/bin/env python
import sys

idsNotFoundFilename = sys.argv[1]
idsNotFoundFile = open(idsNotFoundFilename, 'r')

idsNotFound = []
for gene in idsNotFoundFile:
	gene = gene.rstrip('\r\n')
	idsNotFound.append(gene)

idsNotFoundFile.close()

print "Searching for", len(idsNotFound), "ids."

geneAnnotationsFilename = "gene_annotations.txt"
geneAnnotationsFile = open(geneAnnotationsFilename, 'r')

proteinCodingGenes = []
for line in geneAnnotationsFile:
	line = line.rstrip('\r\n')

	columns = line.split('\t')

	ensgId = columns[0]
	ensgId = ensgId.split('.')
	ensgId = ensgId[0]

	proteinCodingGenes.append(ensgId)

geneAnnotationsFile.close()

print "There are", len(proteinCodingGenes), "protein coding genes in the gene annotations file."

# if the ids in idsNotFound are in the gene annotations file, then they are ids for protein-coding genes
notFoundButProtCoding = []
for idNotFound in idsNotFound:
	if idNotFound in proteinCodingGenes:
		notFoundButProtCoding.append(idNotFound)

if len(notFoundButProtCoding) == len(idsNotFound):
	print "All genes that weren't found were listed as protein-coding."
else:
	print len(notFoundButProtCoding), "genes in the inputted list of ENSGIds are protein-coding, according to GENCODE."

# # output these protein-coding genes as separate file
# newline = "\n"
# output = ""
# for protCodingGene in notFoundButProtCoding:
# 	output += protCodingGene + newline

# outFilename = "proteinCoding_missingGenes.txt"
# outFile = open(outFilename, 'w')
# outFile.write(output)
# outFile.close()