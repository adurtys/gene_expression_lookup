# Date Created: 24 July 2018
# Date Last Modified: 24 July 2018
# Execution: python centroidPositions.py grouped_snps_filename
# Description: TODO

import sys
groupedSnpsFilename = sys.argv[1]
groupedSnpsFile = open(groupedSnpsFilename, 'r')

snpGroups = {}

for line in groupedSnpsFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	groupNumber = columns[0]
	snp = columns[1]

	snp = snp.split(':')
	chromosome = snp[0]
	snpLocation = int(snp[1])

	if groupNumber in snpGroups:
		snpGroups[groupNumber].append(snpLocation)
	else:
		snpGroups[groupNumber] = [snpLocation]

groupedSnpsFile.close()

# for each group number
	# for each snp
		# figure out nearest gene
	# if every snp in group has same gene, then fine
	# if every snp has different nearest genes
		