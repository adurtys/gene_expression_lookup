# Date Created: 24 July 2018
# Date Last Modified: 3 August 2018
# Execution: python centroidPositions.py grouped_snps_filename nearestGenesFilename
# Description: TODO

#!/usr/bin/env python
import sys, collections

# read in command-line arguments
groupedSnpsFilename = sys.argv[1]
nearestGenesFilename = sys.argv[2]

# create dictionary containing the nearest gene to each snp (key = chromosome, value = nearest gene)
nearestGeneDict = {}

# parse nearestGenesFilename
nearestGenesFile = open(nearestGenesFilename, 'r')
for line in nearestGenesFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	snp = columns[0]
	nearestGene = columns[1] # TODO: FIX THIS for if snp is equidistant!

	nearestGeneDict[snp] = nearestGene
nearestGenesFile.close()

# create dictionary that will store groups as keys and the snps in each group as values
snpGroupsDict = {}

# parse groupedSnpsFile
groupedSnpsFile = open(groupedSnpsFilename, 'r')
for line in groupedSnpsFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	groupNumber = int(columns[0])
	snp = columns[1]

	# process the snp
	snp = snp.split(':')
	chromosome = snp[0]
	snpLocation = int(snp[1])
	snpName = chromosome + ":" + str(snpLocation)

	# create list of snps in each group
	snpsInEachGroup = []

	if groupNumber in snpGroupsDict: # group is already stored in the dictionary
		# add the snp to the list of snps contained in that group
		snpGroupsDict[groupNumber].append(snpName)
	else: # group isn't already in the dictionary
		# add snp to the list
		snpsInEachGroup.append(snpName)
		snpGroupsDict[groupNumber] = snpsInEachGroup

groupedSnpsFile.close()

for group in snpGroupsDict:
	# create dictionary containing the nearest genes for each snp in the group
	nearestGeneToSnps = {}

	for snpInGroup in snpGroupsDict[group]:
		if snpInGroup in nearestGeneDict:
			nearestGeneToSnps[snpInGroup] = nearestGeneDict[snpInGroup]
		else:
			print "ERROR (centroidPositions.py line 67): snp doesn't have a nearest gene in the nearest genes file."
			# TODO: handle this case!

	# the values in the nearestGeneToSnpsDict are the nearest genes to each snp in the group
	nearestGenes = nearestGeneToSnps.values()

	uniqNearestGenes = collections.Counter(nearestGenes).keys()
	print "There are", len(uniqNearestGenes), "for snps in group", group

	# if len(uniqNearestGenes) != 1:
		# snps in the group don't all have the same gene
		# TODO: fix this














	




