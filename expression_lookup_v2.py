# Date Created: 28 June 2018
# Date Last Modified: 6 July 2018
# Execution: python expression_lookup_v1.py snp numGenes distance threshold
# snp is a string of the format "chr#:location" for the snp to be looked up in "gene_annotations.txt" file.
# numGenes is an int representing the number of closest genes on either side of the snp that should be analyzed 
# with respect to expression in the tissues listed in "GTEx.tstat.tsv" file.
# distance is a float representing the max distance (in kilobp) from the snp that the closest genes can be
# threshold is a float (decimal) representing the percentage of top rank-ordered t-statistics that should be considered
# as "highly expressed" within a particular tissue.
# TODO: This program ...

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
		print number, "wasn't found in the list."

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
		print "number was found at index:", mid
		
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

# check to make sure file was run with correct number of arguments
if len(sys.argv) != 5:
	print "ERROR: Incorrect number of command-line arguments!"

# read in arguments from user
snp = sys.argv[1]
numGenes = int(sys.argv[2])
distance = float(sys.argv[3])
threshold = float(sys.argv[4])

# for testing!
print "snp:", snp, "numGenes:", numGenes, "distance:", distance, "threshold:", threshold

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

	# remove decimal at end of ENSGID
	ensgId = columns[1]
	ensgId = ensgId.split('.')
	geneId = ensgId[0]

	geneName = columns[2]
	# geneFeature = columns[3]
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

print "Downstream Genes:", downstreamGenes
print "Upstream Genes:", upstreamGenes

# read in the normalized tissue expression file
inFilename2 = "normalizedGTEx.tstat.txt"
inFile2 = open(inFilename2, 'r')

# store column labels in header row
headerLine = inFile2.readline()
headerLine = headerLine.rstrip('\r\n')
headers = headerLine.split('\t')


# calculate number of genes to analyze for tissue expression
uniqueDownstreamGenes = 0
for gene in downstreamGenes:
	if gene not in upstreamGenes:
		uniqueDownstreamGenes += 1

genesForAnalysis = len(upstreamGenes) + uniqueDownstreamGenes
print "genes for analysis:", genesForAnalysis

ids = []

# create a matrix to hold tissue expressions for the genes to be searched
matrix = [[] for gene in range(genesForAnalysis)]

totalGenes = 0
for line in inFile2:
	line = line.rstrip('\r\n')
	tissues = line.split("\t")

	# first column contains GeneID, not tissue
	numTissues = len(tissues) - 1

	# parse for ENSGIDs in upstreamGenes or downstreamGenes	
	if (tissues[0] in downstreamGenes) or (tissues[0] in upstreamGenes):
		ids.append(tissues[0])
		
		# store the ranks for this ID
		for i in range(genesForAnalysis):
			for j in range(numTissues):
				matrix[i].append(int(tissues[j + 1]))
			print j
	totalGenes += 1

# len(ids) should be the same as genesForAnalysis
print "there are", len(ids), "ids for which to look up tissue expression:", ids

numTissues = len(matrix[0])
print "numTissues:", numTissues

# determine whether expression meets threshold for high expression
numTopRankingGenes = threshold * totalGenes
critRank = totalGenes - numTopRankingGenes

expressionMatrix = [[0 for tissue in range(numTissues)] for gene in range(genesForAnalysis)]

# if the rank for expression of a gene in a particular tissue type is greater than critRank, then that gene
# will be considered highly expressed in that tissue
for i in range(genesForAnalysis):
	for j in range(numTissues):
		if matrix[i][j] >= critRank:
			expressionMatrix[i][j] = 1

# create new file containing expressionMatrix
outFilename = "geneExpressionLookupResults.txt"
outFile = open(outFilename, 'w')

tab = "\t"
newline = "\n"

# write header line onto new file
headerLineOutput = headerLine + newline
outFile.write(headerLineOutput)

for i in range(len(ids)):
	output = ids[i] + tab
	for j in range(numTissues):
		if j < (numTissues - 1):
			output += str(expressionMatrix[i][j]) + tab
		else:
			output += str(expressionMatrix[i][j]) + newline

	outFile.write(output)

outFile.close()
inFile2.close()