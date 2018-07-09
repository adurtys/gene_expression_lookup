# Date Created: 27 June 2018
# Date Last Modified: 28 June 2018
# Execution: python GTF_processing_v3
# This program proceses GENCODE Comprehensive Gene Annotation GTF File for the start and end locations of protein-coding genes. 
# The output is a file containing locations of protein-coding genes, and maintains both GeneID (ENSG#) and name. 

#!/usr/bin/env python

# read in the GTF file
inFilename = "gencode.v19.annotation.gtf"
inFile = open(inFilename, 'r')

# create data structures to hold information pertaining to each gene
uniqGenes = []
chromosomeDict = {}
geneNameDict = {}
geneStartDict = {}
geneEndDict = {}

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

	geneName = otherInfo[4]
	geneName = geneName.strip('gene_name ')
	geneName = geneName.strip('"')

	geneType = otherInfo[2]
	geneType = geneType.strip('gene_type ')
	geneType = geneType.strip('"')

	# parse file and extract information for protein-coding genes
	if geneType == "protein_coding":
		# create list of unique protein-coding genes
		if geneId not in uniqGenes:
			print "a unique gene, ", geneId, ", has been found."
			uniqGenes.append(geneType)

		# edit dictionaries containing gene information using uniqGenes list as keys
		if geneId not in chromosomeDict:
			chromosomeDict[geneId] = chromosome

		if geneId not in geneNameDict:
			geneNameDict[geneId] = geneName
		
		if geneId not in geneStartDict:
			geneStartDict[geneId] = geneStart
		# only modify geneStart if new starting position is smaller than current one
		elif (geneId in geneStartDict) and (geneStart < geneStartDict[geneId]):
			geneStartDict[geneId] = geneStart

		if geneId not in geneEndDict:
			geneEndDict[geneId] = geneEnd
		# only modify geneEnd if new ending position is larger than the current one
		elif (geneId in geneEndDict) and (geneEnd > geneEndDict[geneId]):
			geneEndDict[geneId] = geneEnd

print "finished reading in file!"
inFile.close()

# create a new file for start and end positions of only protein-coding genes
outFilename = "gene_annotations_v2.txt"
outFile = open(outFilename, 'w')

tab = "\t"
newline = "\n"

# convert start and end positions back to strings
for geneId in geneStartDict:
	geneStartDict[geneId] = str(geneStartDict[geneId])
for geneId in geneEndDict:
	geneEndDict[geneId] = str(geneEndDict[geneId])

print "writing file now"
for geneId in uniqGenes:
	 outFile.write(geneId)
	 outFile.write(tab)
	 outFile.write(geneNameDict[geneId])
	 outFile.write(tab)
	 outFile.write(geneStartDict[geneStart])
	 outFile.write(tab)
	 outFile.write(geneEndDict[geneEnd])
	 outFile.write(newline)

outFile.close()
