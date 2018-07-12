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

geneAnnotationsFile.close()

print "There are", len(proteinCodingGenes), "protein coding genes in the gene_annotations file."
print "There are", len(idsNotFound), "ids that weren't found in the GTEx file."

# if the ids in idsNotFound are in the gene annotations file, then they are ids for protein-coding genes

notFoundButProtCoding = []
for gene in proteinCodingGenes:
	if gene in idsNotFound:
		notFoundButProtCoding.append(gene)

print "There are", len(notFoundButProtCoding), "genes that don't have tissue-expression data in GTEx file but are protein-coding according to GENCODE."
print notFoundButProtCoding

