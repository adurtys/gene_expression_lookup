# Date Created: 20 August 2018
# Date Last Modified: 20 August 2018
# Execution: python expression_ML_table.py expressionVectorTableFilename
# argv1: filename for file containing centroid snp gene expression lookup results
# Description: TODO

#!/usr/bin/env python
import sys

expressionVectorTableFilename = sys.argv[1]
expressionVectorTableFile = open(expressionVectorTableFilename, 'r')

outFilename = "expression_ML_table.txt"
outFile = open(outFilename, 'w')

# make dictionary storing expression vectors --> key = snp, value = expression vector
expressionVectorDictionary = {}

# store column lables in header row
headerLine = expressionVectorTableFile.readline()
headerLine = headerLine.rstrip('\r\n')
headers = headerLine.split('\t')

# exclude first column when counting the number of tissues present in the file
numTissues = len(headers) - 1

# populate expression vector dictionary
for line in expressionVectorTableFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	snp = columns[0]
	
	vector = []
	for i in range(numTissues):
		vector.append(columns[i + 1])

	expressionVectorDictionary[snp] = vector

# sort dictionary to maintain original order
sortedSnps = sorted(expressionVectorDictionary)

tab = "\t"
newline = "\n"

# create output file
output = ""
for i in range(len(sortedSnps)):
	key = sortedSnps[i]
	output += key + tab

	for j in range(numTissues):
		if j < (numTissues - 1):
			output += expressionVectorDictionary[key][j] + tab
		else: # create new line at the end of each vector
			output += expressionVectorDictionary[key][j] + newline

outFile.write(output)

expressionVectorTableFile.close()
outFile.close()