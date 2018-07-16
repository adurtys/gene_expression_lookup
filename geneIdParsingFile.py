#!/usr/bin/env python

inFilename = "geneExpressionLookupResults_index_snps_used.txt"
inFile = open(inFilename, 'r')

for line in inFile:
	line = line.rstrip('\r\n')

	columns = line.split('')