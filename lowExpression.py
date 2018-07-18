#!/usr/bin/env python

genesToLookupFilename = "./genesInGTEx_v6p.txt"
expressionFilename = "/project/voight_datasets/GTEx_V6p/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct"

genesToLookupFile = open(genesToLookupFilename, 'r')
expressionFile = open(genesToLookupFilename, 'r')

genesToLookup = {}
for gene in genesToLookupFile:
	gene = gene.rstrip('\r\n')
	genesToLookup[gene] = 0
genesToLookupFile.close()

expressionFileVectors = [[] for gene in range(len(genesToLookup))]
geneIndex = 0
for expressionVector in expressionFile:
	expressionVector = expressionVector.rstrip('\r\n')
	tissues = expressionVector.split('\t')

	numColumns = len(tissues)

	geneId = tissues[0]
	geneId = geneId.split('.')
	geneId = geneId[0]

	if geneId in genesToLookup:
		for i in range(numColumns):
			expressionFileVectors[geneIndex].append(tissues[i])
		geneIndex += 1
expressionFile.close()

tab = "\t"
newline = "\n"

outFilename = "lowExpressionCheck.txt"
outFile = open(outFilename, 'w')

for i in range(len(expressionFileVectors)):
	allZeros = True
	expressionOutput = expressionFileVectors[i][0] + tab

	for j in range(len(expressionFileVectors[0])):
		if (j != 0) and (allZeros == True):
			if int(expressionFileVectors[i][j]) != 0:
				print expressionFileVectors[i][0], expressionFileVectors[i][j]
				allZeros = False

	print "exited inner for loop!"
	expressionOutput += str(allZeros) + newline
	outFile.write(expressionOutput)
outFile.close()