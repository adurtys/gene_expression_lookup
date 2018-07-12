#!/usr/bin/env python

idsNotFoundFilename = "idsNotFound_fullTest.txt"
idsNotFoundFile = open(idsNotFoundFilename, 'r')

idsNotFound = []
for gene in idsNotFoundFile:
	idsNotFound.append(gene)

idsNotFoundFile.close()

geneAnnotationsFilename = "gene_annotations.txt"
geneAnnotationsFile = open(geneAnnotationsFilename, 'r')

proteinCodingGenes = []
for line in geneAnnotationsFile:
	line = line.rstrip('\r\n')

	columns = line.split('\t')

	geneId = columns[1]
	geneId = geneId.split('.')
	ensgId = geneId[0]

	proteinCodingGenes.append(ensgId)

notFoundButProtCoding = []
for idNotFound in idsNotFound:
	if idNotFound in proteinCodingGenes:
		notFoundButProtCoding.append(idNotFound)

print "There are", len(notFoundButProtCoding), "genes that don't have tissue-expression data in GTEx file but are protein-coding according to GENCODE."
print notFoundButProtCoding