# Date Created: 28 June 2018
# Date Last Modified: 17 July 2018
# Execution: python expression_lookup.py snpFilename geneAnnotationsFilename tstatFilename numGenes distance threshold
# 	outFilename nearestGenesFilename lostSnpsFilename missingGenesFilename
# 		numGenes is an int representing the number of closest genes on either side of the snp that should be analyzed 
# 			with respect to expression in the tissues listed in "GTEx.tstat.tsv" file.
# 		distance is a float representing the max distance (in kilobp) from the snp that the closest genes can be
# 		threshold is a float (decimal) representing the percentage of top rank-ordered t-statistics that should be considered
# 			as "highly expressed" for each tissue.
# 		outFilename is a string representing the name of the file that will contain the results of the gene expression lookup
# 		lostSnpsFilename is a string representing the name of 
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
if len(sys.argv) != 11:
	print "ERROR: Incorrect number of command-line arguments!"

# read in file containing snps
snpFilename = sys.argv[1]
snpFile = open(snpFilename, 'r')

# read in arguments from user
numGenes = int(sys.argv[4])
distance = float(sys.argv[5])
threshold = float(sys.argv[6])

# create new file that will contain the results of the expression lookup
outFilename = sys.argv[7]
outFile = open(outFilename, 'w')

# create new file that will contain the nearest genes associated with each snp
nearestGenesFilename = sys.argv[8]
nearestGenesFile = open(nearestGenesFilename, 'w')

# create new file that will contain snps that were lost due to lack of tissue expression data for their nearest gene in GTEx file
lostSnpsFilename = sys.argv[9]

# create new file that will contain genes associated with snps in the snpFile that didn't have tissue expression data
missingGenesFilename = sys.argv[10]

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
	inFilename = sys.argv[2]
	inFile = open(inFilename, 'r')

	# create dictionary containing relevant information pertaining to each gene
	genes = {}

	for line in inFile:
		line = line.rstrip('\r\n')
		
		# split line on tab, return list of columns
		columns = line.split('\t')

		geneId = columns[0]
		geneName = columns[1]
		geneChrom = columns[2]
		geneStart = int(columns[3])
		geneEnd = int(columns[4])

		geneInfo = []

		# parse for chromosome specified by user
		if geneChrom == chromosome:
			# populate geneInfo with relevant info
			geneInfo.append(geneName)
			geneInfo.append(geneStart)
			geneInfo.append(geneEnd)

			genes[geneId] = geneInfo

	inFile.close()

	print "The chromosome to be searched is:", chromosome, "\nThere are", len(genes), "genes on this chromosome."

	# create dictionaries of nearby genes with respect to gene start and end locations
	startLocations = {}
	endLocations = {}
	for gene in genes:
		startLocations[gene] = genes[gene][1]
		endLocations[gene] = genes[gene][2]

	# sort start and end location dictionaries
	sortedStartLocations = sorted(startLocations.values())
	sortedEndLocations = sorted(endLocations.values())

	# search for nearest genes with respect to start location
	print "Searching for genes by start locations."
	nearbyStartLocations = binarySearch(sortedStartLocations, snpLocation, numGenes, distance, first = 0, last = None)
	genesByStart = []

	for position in nearbyStartLocations:
		# obtain geneId from the position value
		genesByStart.append(startLocations.keys()[startLocations.values().index(position)])

	# search for nearest genes with respect to end location
	print "Searching for genes by end locations."
	nearbyEndLocations = binarySearch(sortedEndLocations, snpLocation, numGenes, distance, first = 0, last = None)
	genesByEnd = []

	for position in nearbyEndLocations:
		genesByEnd.append(endLocations.keys()[endLocations.values().index(position)])

	print "Genes From Start Position:", genesByStart
	print "Genes From End Position:", genesByEnd

	# determine the closest genes to the snp
	distanceFromSnpDict = {}
	for gene in genesByStart:
		for geneId, startLocation in startLocations.items():
			if geneId == gene:
				distanceFromSnpDict[gene] = abs(startLocation - snpLocation)
		print gene, "is", distanceFromSnpDict[gene], "bp away from the snp."
	for gene in genesByEnd:
		if gene not in distanceFromSnpDict:
			for geneId, endLocation in endLocations.items():
				if geneId == gene:
					distanceFromSnpDict[gene] = abs(endLocation - snpLocation)
			print gene, "is", distanceFromSnpDict[gene], "bp away from the snp."
		else:
			# gene is already in the dictionary
			# modify distance in dictionary only if new distance is less than current distance
			for geneId, endLocation in endLocations.items():
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
			else: # len(closestDistances) == 0
				print "No gene is within 1mbp on either side of the snp."
	else: # genes have identical distances from snp
		print "The snp is equidistant from genes."

		for i in range(len(closestDistances)):
			# append one more entry to genesForAnalysis (numGenes + 1) than if no duplicated values
			# (assumes only equidistant from two gene)
			if closestDistances[i] < closestDistances[numGenes + 1]:
				if closestDistances[i] in duplicates:
					# entry is a duplicated value --> add both entry and subsequent entry

					# find gene corresponding to the distance in distanceFromSnpDict
					index = distanceFromSnpDict.values().index(closestDistances[i])
					gene1 = distanceFromSnpDict.keys()[index]
					gene2 = distanceFromSnpDict.keys()[index + 1]

					genesForAnalysis.append(gene1)
					genesForAnalysis.append(gene2)
				else: 
					# entry is not a duplicated value --> add only the entry
					index = distanceFromSnpDict.values().index(closestDistances[i])
					gene = distanceFromSnpDict.keys()[index]

					genesForAnalysis.append(gene)

	# read in the normalized tissue expression file
	inFilename2 = sys.argv[3]
	inFile2 = open(inFilename2, 'r')

	# store column labels in header row
	headerLine = inFile2.readline()
	headerLine = headerLine.rstrip('\r\n')
	headers = headerLine.split('\t')

	# subtract 1 because first column contains GeneID
	numTissues = len(headers) - 1

	if len(genesForAnalysis) != 0:
		# there are genes to analyze
		numGenesForAnalysis = len(genesForAnalysis)
		print "Analyzing", numGenesForAnalysis, "genes:", genesForAnalysis

		idsForTissueExpression = []

		# create a matrix to hold tissue expressions for the genes to be searched
		matrix = [[] for gene in range(numGenesForAnalysis)]

		totalGenes = 0
		geneIndex = 0
		for line in inFile2:
			line = line.rstrip('\r\n')
			tissues = line.split('\t')

			# parse for ENSGIDs in genesForAnalysis
			if tissues[0] in genesForAnalysis:
				idsForTissueExpression.append(tissues[0])

				# store the ranks for this ID
				for i in range(numTissues):
					matrix[geneIndex].append(int(tissues[i + 1]))
				geneIndex += 1
			
			totalGenes += 1

		print "There are", len(idsForTissueExpression), "ids for which to look up tissue expression:", idsForTissueExpression

		if len(idsForTissueExpression) != numGenesForAnalysis:
			# add snp to list of snps whose nearest gene doesn't have tissue expression data in the GTEx file
			snpName = snp[0] + ":" + snp[1]
			snpsNoTissueExp.append(snpName)

			for item in genesForAnalysis:
				# if the GeneId isn't in the GTEx file, add it to the list of ids without tissue expression data
				# if it isn't already in the list
				if (item not in idsForTissueExpression) and (item not in idsWithoutTissueExpData):
					idsWithoutTissueExpData.append(item)

		# determine whether expression meets threshold for high expression
		numTopRankingGenes = threshold * totalGenes
		critRank = totalGenes - numTopRankingGenes

		expressionMatrix = [[0 for tissue in range(numTissues)] for gene in range(numGenesForAnalysis)]

		# if the rank for expression of a gene in a particular tissue type is greater than critRank, then that gene
		# will be considered highly expressed in that tissue
		for i in range(numGenesForAnalysis):
			# for every nearby gene that was found
			for j in range(numTissues):
				# edit experssionMatrix for every tissue
				if i < len(idsForTissueExpression):
					# gene has tissue expression data in the GTEx file
					if matrix[i][j] >= critRank:
						expressionMatrix[i][j] = 1
				else:
					# gene does not have tissue expression data in the GTEx file
					# print out zeros for every tissue
					expressionMatrix[i][j] = 0 # probably happens automatically but written explicitly

		# create expression vector - only output one line per snp
		# if more than one gene for analysis, combine the tissue expression vectors for both genes
		expressionVector = [0 for tissue in range(numTissues)]

		for i in range(numGenesForAnalysis):
			for j in range(numTissues):
				if matrix[i][j] == 1:
					# update any tissue's expression if the expression of any gene meets threshold
					expressionVector[j] = 1

		# write header line onto new file once (only for first snp)
		if numSnps == 1:
			headerLineOutput = headerLine + newline
			outFile.write(headerLineOutput)

		snpName = snp[0] + ":" + snp[1]

		# write out expression vector
		output = snpName + tab
		for i in range(numTissues):
			if i < (numTissues - 1):
				output += str(expressionVector[i]) + tab
			else: # add newline if last entry in the vector
				output += str(expressionVector[i]) + newline
		outFile.write(output)

		# write out file containing nearest genes corresponding to each snp
		nearestGenesOutput = snpName + tab
		for i in range(numGenesForAnalysis):
			if i < (genesForAnalysis - 1):
				nearestGenesOutput += genesForAnalysis[i] + tab + genes[genesForAnalysis[i]][0] + tab
			else: # add newline if last entry in the vector
				nearestGenesOutput += genesForAnalysis[i] + tab + genes[genesForAnalysis[i]][0] + newline

		nearestGenesFile.write(nearestGenesOutput)

	else:
		# there are no genes to analyze
		numSnpsNoGenes += 1

		snpName = snp[0] + ":" + snp[1]
		output = snpName + tab
		for i in range(numTissues):
			if i < (numTissues - 1):
				output += "0" + tab
			else:
				output += "0" + newline

		outFile.write(output)

		nearestGenesOutput = "No genes were found within 1mbp on either side of the snp that was searched." + newline
		nearestGenesFile.write(nearestGenesOutput)

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
nearestGenesFile.close()
snpFile.close()