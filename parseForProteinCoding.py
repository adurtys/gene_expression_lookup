#!/usr/bin/env python

idsNotFoundFilename = "idsNotFound_fullTest.txt"
idsNotFoundFile = open(idsNotFoundFilename, 'r')

idsNotFound = []
for gene in idsNotFoundFile:
	gene = gene.rstrip('\r\n')
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

# newline = "\n"
# output = ""
# for protCodingGene in proteinCodingGenes:
# 	output += protCodingGene + newline

# outFilename = "protein_coding_genes.txt"
# outFile = open(outFilename, 'w')
# outFile.write(output)
# outFile.close()

print "There are", len(proteinCodingGenes), "protein coding genes in the gene annotations file."
print "There are", len(idsNotFound), "ids that weren't found in the GTEx file."

# if the ids in idsNotFound are in the gene annotations file, then they are ids for protein-coding genes

notFoundButProtCoding = []
for idNotFound in idsNotFound:
	if idNotFound in proteinCodingGenes:
		notFoundButProtCoding.append(idNotFound)

# for i in range(10):
# 	print idsNotFound[i]
# 	isProteinCoding = False
# 	if idsNotFound[i] in proteinCodingGenes:
# 		isProteinCoding = True
# 	print isProteinCoding

if len(notFoundButProtCoding) == len(idsNotFound):
	print "All genes that weren't found were listed as protein-coding."
else:
	print "There are", len(notFoundButProtCoding), "genes that don't have tissue-expression data in GTEx file but are protein-coding according to GENCODE."

print notFoundButProtCoding