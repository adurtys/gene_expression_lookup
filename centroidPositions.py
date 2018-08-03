# Date Created: 24 July 2018
# Date Last Modified: 3 August 2018
# Execution: python centroidPositions.py groupedSnpsFilename geneAnnotationsFilename
# Description: TODO

#!/usr/bin/env python
import sys

# read in command-line arguments
groupedSnpsFilename = sys.argv[1]
geneAnnotationsFilename = sys.argv[2]

# create dictionary containing relevant information pertaining to each gene in the gene annotations file (key = geneId; value = list of info [geneChrom, geneName, geneStart, geneEnd])
geneAnnotations = {}

# populate this dictionary by parsing gene annotations file
geneAnnotationsFile = open(geneAnnotationsFilename, 'r')

for line in geneAnnotationsFile:
	line = line.rstrip('\r\n')
	
	# split line on tab, return list of columns
	columns = line.split('\t')
	geneId = columns[0]
	geneName = columns[1]
	geneChrom = columns[2]
	geneStart = int(columns[3])
	geneEnd = int(columns[4])

	geneInfo = []

	# populate geneInfo with relevant info
	geneInfo.append(geneChrom)
	geneInfo.append(geneName)
	geneInfo.append(geneStart)
	geneInfo.append(geneEnd)

	geneAnnotations[geneId] = geneInfo

geneAnnotationsFile.close()
print "Finished reading in gene annotations file."

# create dictionary that will store groups as keys and the snps in each group as values
snpGroupsDict = {}

# parse groupedSnpsFile
groupedSnpsFile = open(groupedSnpsFilename, 'r')
for line in groupedSnpsFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	groupNumber = int(columns[0])
	snp = columns[1]

	# create list of snps in each group
	snpsInEachGroup = []

	if groupNumber in snpGroupsDict: # group is already stored in the dictionary
		# add the snp to the list of snps contained in that group
		snpGroupsDict[groupNumber].append(snp)
	else: # group isn't already in the dictionary
		# add snp to the list
		snpsInEachGroup.append(snp)
		snpGroupsDict[groupNumber] = snpsInEachGroup
groupedSnpsFile.close()
print "There are", len(snpGroupsDict), "groups on this chromosome."

# create dictionary that will store the group number as the key and the centroid snp as its value
centroidSnps = {}

for group in snpGroupsDict:
	groupSnpLocations = []

	for groupedSnp in snpGroupsDict[group]:
		# process the snp
		groupedSnp = groupedSnp.split(':')
		chromosome = groupedSnp[0]
		snpLocation = int(groupedSnp[1])

		groupSnpLocations.append(snpLocation)

	minSnpLocation = min(groupSnpLocations)
	maxSnpLocation = max(groupSnpLocations)
	centroidSnpLocation = (minSnpLocation + maxSnpLocation) / 2.0

	centroidSnps[group] = centroidSnpLocation

tab = "\t"
newline = "\t"

output = ""

outFilename = "centroidSnps_" + groupedSnpsFilename
outFile = open(outFilename, 'w')

for group in centroidSnps:
	output += group + tab + str(centroidSnps[group]) + newline

outFile.write(output)


	




