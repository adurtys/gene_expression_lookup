#!/usr/bin/env python

genesToLookupFilename = "./genesInGTEx_v6p.txt"
expressionFilename = "/project/voight_datasets/GTEx_V6p/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct"

genesToLookupFile = open(genesToLookupFilename, 'r')
expressionFile = open(expressionFilename, 'r')

genesToLookup = {}
for gene in genesToLookupFile:
	gene = gene.rstrip('\r\n')
	genesToLookup[gene] = 0
genesToLookupFile.close()

# skip first two lines
expressionFile.readline()
expressionFile.readline()

headerline = expressionFile.readline()
headerline = headerline.rstrip('\r\n')
headers = headerline.split('\t')

numColumns = len(headers)

expressionFileVectors = [[] for gene in range(len(genesToLookup))]

geneIndex = 0
for expressionVector in expressionFile:
	expressionVector = expressionVector.rstrip('\r\n')
	tissues = expressionVector.split('\t')

	# print "numTissueColumns:", len(tissues)

	geneId = tissues[0]
	geneId = geneId.split('.')
	geneId = geneId[0]

	geneName = tissues[1]

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
		# skip first two columns, containing geneId and geneName, respectively
		if (j > 1) and (allZeros == True):
			if int(expressionFileVectors[i][j]) != 0:
				print "expression isn't all zero for", expressionFileVectors[i][0], ":", expressionFileVectors[i][j]
				allZeros = False

	expressionOutput += str(allZeros) + newline
	outFile.write(expressionOutput)
outFile.close()