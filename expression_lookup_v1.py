# Date Created: 28 June 2018
# Date Last Modified: 2 July 2018
# Execution: python expression_lookup_v1.py snp numGenes distance
# snp is a string of the format "chr#:location" for the snp to be looked up in "gene_annotations.txt" file.
# numGenes is an int representing the number of closest genes on either side of the snp that should be analyzed 
# with respect to expression in the tissues listed in "GTEx.tstat.tsv" file.
# distance is a flot representing the max distance (in kilobp) from the snp that the closest genes can be
# TODO: This program ...

#!/usr/bin/env python

# read in arguments from user
snp = sys.argv[1]
numGenes = int(sys.argv[2])
distance = float(sys.argv[3])
if distance > 1000:
	distance = 1000

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

for line in file:
	line = line.rstrip('\r\n')
	
	# split line on tab, return list of columns
	columns = line.split('\t')
	
	geneChrom = columns[0]
	geneId = columns[1]
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

print "The chromosome to be searched is:", chromosome, ".\nThere are", len(nameDict), "genes on this chromosome."

# Input: a sorted list of locations, a location to be searched, the number of nearest genes to return, the distance
# within which to search for nearby genes, the first index (default = 0), and the last index (default = None) 
# Output: if the searched number was found, the function will return a list containing the locations of the searched number 
# and the next closest numbers, both greater and less than the number that was searched. If the searched number was not
# found, then the function will return a list of the nearest numbers. The number of items in this list is the number
# specified by the user, as long as the genes are within the specified distance.
# Description: this function is a recursive binary search function that will return locations for the number of closest
# genes specified by the user, given that these genes are within 1mb from the snp. If the snp itself occurs at a location,
# then this location is included in the list of locations. numGenes is the number of locations less than and greater 
# than the value to be returned.
def binarySearch(sortedList, number, numGenes, distance, first = 0, last = None):
	locations = []

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

	# number is larger than the number at mid
	return binarySearch(sortedList, number, numGenes, distance, mid + 1, last)

# create list of nearby genes downstream from the snp
sortedStartDictionary = sorted(startDict)
nearbyStartLocations = binarySearch(sortedStartDictionary, snpLocation, numGenes, distance, first = 0, last = None)
downstreamGenes = []
for position in nearbyStartLocations:
	downstreamGenes.append(startDict[position])

# create list of nearby genes upstream from the snp
sortedEndDictionary = storted(endDict)
nearbyEndLocations = binarySearch(sortedEndDictionary, snpLocation, numGenes, distance, first = 0, last = None)
upstreamGenes = []
for posiition in nearbyEndLocations:
	upstreamGenes.append(endDict[posiition])

# read in tissue expression file (tab-separated)




