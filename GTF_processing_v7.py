# Date Created: 27 June 2018
# Date Last Modified: 28 June 2018
# Execution: python GTF_processing_v3
# This program proceses GENCODE Comprehensive Gene Annotation GTF File for the start and end locations of protein-coding genes. 
# The output is a file containing locations of protein-coding genes, and maintains both GeneID (ENSG#) and name. 

#!/usr/bin/env python

# read in the GTF file (TODO: for now, it is the shortened version)
inFilename = "firstLines.txt"
inFile = open(inFilename, 'r')

# create data structures to hold information pertaining to each gene
genes = []
geneNames = []
chromNums = []
geneTypes = []
geneStartLocations = []
geneEndLocations = []

# skip first five lines
for i in range(5):
        inFile.readline()

for line in inFile:
   	line = line.rstrip('\r\n')
	
	# split line on tab, return list of columns
	columns = line.split('\t')
    
	chromosome = columns[0]
	geneStart = int(columns[3])
	geneEnd = int(columns[4])
	
	# process otherInfo column by getting rid of key label and quotes
	otherInfo = columns[8]
	otherInfo = otherInfo.split(';')

	geneId = otherInfo[0]
	geneId = geneId.strip('gene_id ')
	geneId = geneId.strip('"')

	name = otherInfo[4]
	name = name.strip('gene_name ')
	name = name.strip('"')

	geneType = otherInfo[2]
	geneType = geneType.strip('gene_type ')
	geneType = geneType.strip('"')

	# parse for protein-coding genes
	if geneType == "protein_coding":
		# add info to their respective data structures
		genes.append(geneId)
		geneNames.append(name)
		chromNums.append(chromosome)
		geneTypes.append(geneType)
		geneStartLocations.append(geneStart)
		geneEndLocations.append(geneEnd)

inFile.close()

numGenes = len(genes)
print "There are", numGenes, "total protein-coding genes in this file."

# create new data structures to store information for unique protein-coding genes
uniqGenes = []
uniqGeneName = {}
uniqGeneChrom = {}
uniqGeneStart = {}
uniqGeneEnd = {}

for index in range(numGenes):
	# parse protein-coding gene list for unique genes
	if genes[index] not in uniqGenes:
		uniqGenes.append(genes[index])

	# edit dictionaries containing gene info
	key = genes[index]

	uniqGeneName[key] = geneNames[index]
	uniqGeneChrom[key] = chromNums[index]

	# add start location if not already in the dictionary
	if key not in uniqGeneStart:
		uniqGeneStart[key] = geneStartLocations[index]
	# only modify start location if new position is smaller than current one
	elif geneStartLocations[index] < uniqGeneStart[key]:
		uniqGeneStart[key] = geneStartLocations[index]

	# add end location if not already in the dictionary
	if key not in uniqGeneEnd:
		uniqGeneEnd[key] = geneEndLocations[index]
	# only modify end location if new position is larger than current one
	elif geneEndLocations[index] > uniqGeneEnd[key]:
		uniqGeneEnd[key] = geneEndLocations[key]

numUniqGenes = len(uniqGenes)
print "There are ", numUniqGenes, "unique protein-coding genes in this file."

for uniqId in uniqGenes:
	print uniqId, "\t", uniqGeneName[uniqId], "\t", uniqGeneChrom[uniqId], "\t", uniqGeneStart[uniqId], "\t", uniqGeneEnd[uniqId]


# create a new file for start and end positions of only protein-coding genes
# outFilename = "gene_annotations_v2.txt"
# outFile = open(outFilename, 'w')

# tab = "\t"
# newline = "\n"

# # convert start and end positions back to strings
# for geneId in geneStartDict:
# 	geneStartDict[geneId] = str(geneStartDict[geneId])
# for geneId in geneEndDict:
# 	geneEndDict[geneId] = str(geneEndDict[geneId])

# string = geneId + tab + geneNameDict[geneId] + tab + geneStartDict[geneId] ...

# print "writing file now"
# for geneId in uniqGenes:
# 	 outFile.write(geneId)
# 	 outFile.write(tab)
# 	 outFile.write(geneNameDict[geneId])
# 	 outFile.write(tab)
# 	 outFile.write(geneStartDict[geneStart])
# 	 outFile.write(tab)
# 	 outFile.write(geneEndDict[geneEnd])
# 	 outFile.write(newline)

# outFile.close()
