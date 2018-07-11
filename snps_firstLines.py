# Date Created: 11 July 2018
# Date Last Modified: 11 July 2018
# Execution: python snps_firstLines.py numLinesToPrint
# Description: Writes the first numLinesToPrint lines of "index_snps_used.txt" to "snpFile.txt"

#!/usr/bin/env python

import sys

numLinesToPrint = int(sys.argv[1])

# read in entire snps file
inFilename = "index_snps_used.txt"
inFile = open(inFilename, 'r')

outFilename = "snpFile.txt"
outFile = open(outFilename, 'w')

lines = inFile.readlines()

for i in range(numLinesToPrint):
	outFile.write(lines[i])

outFile.close()
inFile.close()