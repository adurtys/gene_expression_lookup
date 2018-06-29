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
chromNums = []
geneNames = []
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

	# add info to their respective data structures
	genes.append(geneId)
	chromNums.append(chromosome)
	geneNames.append(name)
	geneTypes.append(geneType)
	geneStartLocations.append(geneStart)
	geneEndLocations.append(geneEnd)

inFile.close()

print "There are ", len(genes), "total genes in this file."

# create new data structures to store information for unique protein-coding genes
uniqGenes = []
uniqGeneChrom = {}
uniqGeneName = {}
uniqGeneStart = {}
uniqGeneEnd = {}

# parse gene list for unique protein-coding genes
for id in range(len(genes)):
	if (geneTypes[id] == "protein_coding") and (genes[id] not in uniqGenes):
		uniqGenes.append(genes[id])

print "There are ", len(uniqGenes), "unique protein-coding genes in this file."

for uniqId in uniqGenes:
	print uniqId


		#	uniqGeneChrom[id] =  chromNums[id]
		#	uniqGeneName[id]


 	# edit dictionaries containing gene information using uniqGenes list as keys
	# if geneId not in chromosomeDict:
	# 	chromosomeDict[geneId] = chromosome

	# 	if geneId not in geneNameDict:
	# 		geneNameDict[geneId] = geneName
		
	# 	if geneId not in geneStartDict:
	# 		geneStartDict[geneId] = geneStart
	# 	# only modify geneStart if new starting position is smaller than current one
	# 	elif (geneId in geneStartDict) and (geneStart < geneStartDict[geneId]):
	# 		geneStartDict[geneId] = geneStart

	# 	if geneId not in geneEndDict:
	# 		geneEndDict[geneId] = geneEnd
	# 	# only modify geneEnd if new ending position is larger than the current one
	# 	elif (geneId in geneEndDict) and (geneEnd > geneEndDict[geneId]):
	# 		geneEndDict[geneId] = geneEnd

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
