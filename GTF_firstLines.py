# Date Created: 27 June 2018
# Date Last Modified: 28 June 2018
# TODO: Execution: python GTF_firstLines.py
# TODO: This program ...

#!/usr/bin/env python

# read in the GTF file
inFilename = "gencode.v19.annotation.gtf"
outFilename = "firstLines.txt"

inFile = open(inFilename, 'r')
outfile = open(outFilename, 'w')

lines = inFile.readlines()

for i in range(50):
	outFile.write(lines[i])

inFile.close()
outFile.close()