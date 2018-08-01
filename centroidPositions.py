# Date Created: 24 July 2018
# Date Last Modified: 31 July 2018
# Execution: python centroidPositions.py grouped_snps_filename
# Description: TODO

import sys, expression_lookup3

groupedSnpsFilename = sys.argv[1]
groupedSnpsFile = open(groupedSnpsFilename, 'r')

# create dictionary that will store groups as keys and the snps in each group as values
snpGroupsDict = {}

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

# for each group number
	# for each snp
		# figure out nearest gene
	# if every snp in group has same gene, then fine
	# if every snp has different nearest genes
		