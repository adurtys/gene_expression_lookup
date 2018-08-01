# Date Created: 25 July 2018
# Date Last Modified: 31 July 2018
# Execution: python nearestGeneSearch.py positionsOnGenome snpLocation numGenes distanceFromSnp
# argv1: positionsOnGenome is a list containing positions on the genome (locations where genes start or end)
# argv2: snpLocation is the location of the snp for which to find the nearest gene
# argv3: numGenes is the number of nearest genes from the snp for which to search
# argv4: distanceFromSnp is the distance from the snp for which to search for genes (in kbp)
# Description: 

#!/usr/bin/env python
import sys, expression_lookup3

if len(sys.argv) != 5:
	print "ERROR (nearestGeneSearch.py line 12): Incorrect number of command-line arguments!"

# read in command line arguments
positionsOnGenome = sys.argv[1]
snpLocation = int(sys.argv[2])
numGenes = int(sys.argv[3])
distanceFromSnp = int(sys.argv[4])

# create list that sorts the positions specified in positionsOnGenome
sortedPositionsOnGenome = sorted(positionsOnGenome)

# convert distanceFromSnp (kbp) to bp
distance = distance * 1000

locationsOfNearestGenes = []
