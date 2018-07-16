# Date Created: 28 June 2018
# Date Last Modified: 13 July 2018
# Execution: python expression_lookup_final.py snpFilename numGenes distance threshold outFilename lostSnpsFilename missingGenesFilename
# 	numGenes is an int representing the number of closest genes on either side of the snp that should be analyzed 
# 	with respect to expression in the tissues listed in "GTEx.tstat.tsv" file.
# distance is a float representing the max distance (in kilobp) from the snp that the closest genes can be
# threshold is a float (decimal) representing the percentage of top rank-ordered t-statistics that should be considered
# 	as "highly expressed" for each tissue.
# outFilename is a string representing the name of the file that will contain the results of the gene expression lookup
# lostSnpsFilename is a string representing the name of 
# Description: This program finds the genes closest to the snp, then outputs a file of gene expression across tissues in the
# 	GTEx t-statistics file, where 1 = highly expressed, 0 = not highly expressed. This is done for each snp listed in the
# 	"snpFile.txt" file, which contains the snps in the format "chr#:location" to look up in the "gene_annotations.txt" file. 
# Run Time: 0.1 sec for each snp

#!/usr/bin/env python
import sys

# Input: a sorted list of locations, a position to be searched, the number of nearest genes to return, the distance
# within which to search for nearby genes (in kilobp), the first index (default = 0), and the last index (default = None) 
# Output: if the searched position was found, the function will return a list containing the locations of the searched number 
# and the next closest numbers, both greater and less than the number that was searched. If the searched position was not
# found, then the function will return a list of the nearest positions. The number of items in this list is the number of
# nearest genes specified by the user, as long as the genes are within the specified distance.
# Description: this function is a recursive binary search function that will return locations for the number of closest
# genes specified by the user, given that these genes are within 1mb from the snp. If the snp itself occurs at a location,
# then this location is included in the list of locations. numGenes is the number of locations less than and greater 
# than the value to be returned.
def binarySearch(sortedList, number, numGenes, distance, first = 0, last = None):
	locations = []

	# only search for genes within 1mb of the snp
	if distance < 0:
		print "ERROR:", distance, " is an invalid input value for distance. Please enter a valid distance from the snp."
	if distance > 1000:
		distance = 1000

	# convert distance (kilobp) to bp
	distance = distance * 1000

	# this is the first call of the function on the list
	if last is None:
		last = len(sortedList) - 1
	
	# stop recursive search when first > last
	if first > last:
		# create a list of closest numbers
		for i in range(numGenes):
			# check distance between next closest genes on either side
			# only add gene locations to the list if they are within the specified distance
			if (last - i) >= 0:
				distLeft = abs(sortedList[last - i] - number)
				if distLeft <= distance:
					locations.insert(0, sortedList[last - i])
			
			if (first + i) <= len(sortedList) - 1:
				distRight = abs(sortedList[first + i] - number)
				if distRight <= distance:
					locations.append(sortedList[first + i])

		print "Finished searching."
		print "locations:", locations
		return locations

	# mid is the index for searching the number
	mid = (first + last) // 2

	# number was found
	if number == sortedList[mid]:
		print "number was found."	
		# create list of number and closest numbers
		locations.append(sortedList[mid])
		for i in range(1, numGenes + 1):
			# check distance between next closest genes on either side
			# only add gene location if within the specified distance
			if (mid - i) >= 0:
				distLeft = abs(sortedList[mid - i] - number)
				if distLeft <= distance:
					locations.insert(0, sortedList[mid - i])
			
			if (mid + 1) <= len(sortedList) - 1:
				distRight = abs(sortedList[mid + i] - number)
				if distRight <= distance:
					locations.append(sortedList[mid + i])

		print "locations:", locations
		return locations
	
	# number is smaller than the number at mid
	if number < sortedList[mid]:
		return binarySearch(sortedList, number, numGenes, distance, first, mid - 1)
	else:
		# number is larger than the number at mid
		return binarySearch(sortedList, number, numGenes, distance, mid + 1, last)

# Input: a list
# Output: a list of items that are present in the list more than once
# Description: finds items that are present more than once in the inputted list 
def findDuplicates(anyList):
	seen = set()
	duplicates = []

	for item in anyList:
		if item in seen and item not in duplicates:
			duplicates.append(item)
		else:
			seen.add(item)

	return duplicates

# check to make sure file was run with correct number of arguments
if len(sys.argv) != 8:
	print "ERROR: Incorrect number of command-line arguments!"

# read in file containing snps
snpFilename = sys.argv[1]
snpFile = open(snpFilename, 'r')

# read in arguments from user
numGenes = int(sys.argv[2])
distance = float(sys.argv[3])
threshold = float(sys.argv[4])

# create new file that will contain the results of the expression lookup
outFilename = sys.argv[5]
outFile = open(outFilename, 'w')

# create new file that will contain snps that were lost due to lack of tissue expression data for their nearest gene in GTEx file
lostSnpsFilename = sys.argv[6]

# create new file that will contain genes associated with snps in the snpFile that didn't have tissue expression data
missingGenesFilename = sys.argv[7]

tab = "\t"
newline = "\n"

numSnps = 0
	
# create counter variable for snps that didn't have a nearby gene
numSnpsNoGenes = 0

# create lists for snps and genes without corresponding tissue expression data in GTEx file
snpsNoTissueExp = []
idsWithoutTissueExpData = []

for snp in snpFile:
	numSnps += 1
	snp = snp.rstrip('\r\n')

	# process the snp
	snp = snp.split(':')
	chromosome = snp[0]
	snpLocation = int(snp[1])

	# read in gene annotations file
	inFilename = "gene_annotations.txt"
	inFile = open(inFilename, 'r')

	# create data structures containing information pertaining to each gene
	nameDict = {}
	startDict = {}
	endDict = {}
	startLocations = []
	endLocations = []

	for line in inFile:
		line = line.rstrip('\r\n')
		
		# split line on tab, return list of columns
		columns = line.split('\t')
		
		geneChrom = columns[0]
		geneId = columns[1]
		geneName = columns[2]
		geneStart = int(columns[4])
		geneEnd = int(columns[5])

		# parse for chromosome specified by user
		if geneChrom == chromosome:
			# populate arrays to be used for searching
			startLocations.append(geneStart)
			endLocations.append(geneEnd)

			# populate the dictionaries
			nameDict[geneId] = geneName
			# start and end location dictionaries use GeneID as value
			startDict[geneStart] = geneId
			endDict[geneEnd] = geneId

	inFile.close()

	print "The chromosome to be searched is:", chromosome, "\nThere are", len(nameDict), "genes on this chromosome."

	# create list of nearby genes downstream from the snp
	sortedStartDictionary = sorted(startDict)
	print "Searching for genes downstream from the snp."
	nearbyStartLocations = binarySearch(sortedStartDictionary, snpLocation, numGenes, distance, first = 0, last = None)
	downstreamGenes = []
	for position in nearbyStartLocations:
		downstreamGenes.append(startDict[position])

	# create list of nearby genes upstream from the snp
	sortedEndDictionary = sorted(endDict)
	print "Searching for genes upstream from the snp."
	nearbyEndLocations = binarySearch(sortedEndDictionary, snpLocation, numGenes, distance, first = 0, last = None)
	upstreamGenes = []
	for posiition in nearbyEndLocations:
		upstreamGenes.append(endDict[posiition])

	print "Genes From Start Position:", downstreamGenes
	print "Genes From End Position:", upstreamGenes

	# determine the closest genes to the snp
	distanceFromSnpDict = {}
	for gene in downstreamGenes:
		for startLocation, geneId in startDict.items():
			if geneId == gene:
				distanceFromSnpDict[gene] = abs(startLocation - snpLocation)
		print gene, "is", distanceFromSnpDict[gene], "bp away from the snp."
	for gene in upstreamGenes:
		if gene not in distanceFromSnpDict:
			for endLocation, geneId in endDict.items():
				if geneId == gene:
					distanceFromSnpDict[gene] = abs(endLocation - snpLocation)
			print gene, "is", distanceFromSnpDict[gene], "bp away from the snp."
		else:
			# gene is already in the dictionary
			# modify distance in dictionary only if new distance is less than current distance
			for endLocation, geneId in endDict.items():
				if (geneId == gene) and (abs(endLocation - snpLocation) < distanceFromSnpDict[gene]):
					distanceFromSnpDict[gene] = abs(endLocation - snpLocation)
			print gene, "is", distanceFromSnpDict[gene], "bp away from the snp."

	closestDistances = sorted(distanceFromSnpDict.values())
	print "Closest Distances:", closestDistances
	# check if there are any duplicate distances in distanceFromSnpDict
	duplicates = findDuplicates(closestDistances)

	# create list of only numGenes closest genes to the snp
	genesForAnalysis = []

	if len(duplicates) == 0:
		print "The snp is not equidistant from genes."
		# no genes have identical distances from the snp
		for gene in distanceFromSnpDict:
			if len(closestDistances) > 1:
				# more than one gene is near the snp
				# only add the distances less than closestDistances[numGenes]
				if distanceFromSnpDict[gene] < closestDistances[numGenes]:
					genesForAnalysis.append(gene)
			elif len(closestDistances) == 1:
				# only one gene is near the snp
				# analyze this gene
				genesForAnalysis.append(gene)
			else:
				# len(closestDistances) == 0
				print "No gene is within 1mbp on either side of the snp."
	else:
		# genes have identical distances from snp
		print "The snp is equidistant from genes."

		for gene in distanceFromSnpDict:
			# one more entry in list of genes to analyze than if there weren't duplicated values
			if (distanceFromSnpDict[gene] in duplicates) and (distanceFromSnpDict[gene] < closestDistances[numGenes + 1]):
				genesForAnalysis.append(gene)

	# read in the normalized tissue expression file
	inFilename2 = "normalizedGTEx.tstat.txt"
	inFile2 = open(inFilename2, 'r')

	# store column labels in header row
	headerLine = inFile2.readline()
	headerLine = headerLine.rstrip('\r\n')
	headers = headerLine.split('\t')

	if len(genesForAnalysis) != 0:
		# there are genes to analyze
		numGenesForAnalysis = len(genesForAnalysis)
		print "Analyzing", numGenesForAnalysis, "genes:", genesForAnalysis # should be equal to numGenes

		ids = []

		# create a matrix to hold tissue expressions for the genes to be searched
		matrix = [[] for gene in range(numGenesForAnalysis)]

		totalGenes = 0
		geneIndex = 0
		for line in inFile2:
			line = line.rstrip('\r\n')
			tissues = line.split("\t")

			# first column contains GeneID, not tissue
			numTissues = len(tissues) - 1

			# parse for ENSGIDs in genesForAnalysis
			if tissues[0] in genesForAnalysis:
				ids.append(tissues[0])

				# store the ranks for this ID
				for i in range(numTissues):
					matrix[geneIndex].append(int(tissues[i + 1]))
				geneIndex += 1
			
			totalGenes += 1

		print "There are", len(ids), "ids for which to look up tissue expression:", ids

		# len(ids) should be the same as numGenesForAnalysis
		if len(ids) != numGenesForAnalysis:
			# add snp to list of snps whose nearest gene doesn't have tissue experssion data in the GTEx file
			snpName = snp[0] + ":" + snp[1]
			snpsNoTissueExp.append(snpName)

			for item in genesForAnalysis:
				# if the GeneId isn't in the GTEx file, add it to the list of ids without tissue expression data
				# if it isn't already in the list
				if (item not in ids) and (item not in idsWithoutTissueExpData):
					idsWithoutTissueExpData.append(item)

		numTissues = len(matrix[0])

		# determine whether expression meets threshold for high expression
		numTopRankingGenes = threshold * totalGenes
		critRank = totalGenes - numTopRankingGenes

		expressionMatrix = [[0 for tissue in range(numTissues)] for gene in range(numGenesForAnalysis)]

		# if the rank for expression of a gene in a particular tissue type is greater than critRank, then that gene
		# will be considered highly expressed in that tissue
		for i in range(numGenesForAnalysis):
			for j in range(numTissues):
				if matrix[i][j] >= critRank:
					expressionMatrix[i][j] = 1

		# write header line onto new file once (only for first snp)
		if numSnps == 1:
			headerLineOutput = headerLine + newline
			outFile.write(headerLineOutput)

		# create counter variable to determine how many tissue expression vectors have been written out
		numOutputs = 0

		for i in range(len(ids)):
			output = ids[i] + tab + nameDict[ids[i]] + tab
			for j in range(numTissues):
				if j < (numTissues - 1):
					output += str(expressionMatrix[i][j]) + tab
				else:
					output += str(expressionMatrix[i][j]) + newline

			outFile.write(output)
			numOutputs += 1

	else:
		# there are no genes to analyze
		numSnpsNoGenes += 1
		output = "No genes were found within 1mbp on either side of the snp that was searched." + newline
		outFile.write(output)

	print "Finished writing tissue expression lookup results onto output file."
	inFile2.close()

print "There were", numSnps, "snps to search in the snpFile.txt file."
print "There were", numSnpsNoGenes, "snps that did not have a nearby gene to analyze."
print "There were", len(snpsNoTissueExp), "snps that did not have tissue expression data for their nearest gene."
print "There were", len(idsWithoutTissueExpData), "GeneIDs that did not have tissue expression data."

# output lost snps to output file
lostSnpsFile = open(lostSnpsFilename, 'w')

lostSnpsOutput = ""
for lostSnp in snpsNoTissueExp:
	lostSnpsOutput += lostSnp + newline

lostSnpsFile.write(lostSnpsOutput)
lostSnpsFile.close()

# output genes without tissue expression data to output file
missingGenesFile = open(missingGenesFilename, 'w')

missingGenesOutput = ""
for missingGene in idsWithoutTissueExpData:
	missingGenesOutput += missingGene + newline

missingGenesFile.write(missingGenesOutput)
missingGenesFile.close()

outFile.close()
snpFile.close()