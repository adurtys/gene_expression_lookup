# Date Created: 20 August 2018
# Date Last Modified: 20 August 2018
# Execution: python expression_ML_table.py groupedSnpTypesFilename expressionVectorTableFilename
# argv1: filename for file containing snp type for each group of snps
# argv2: filename for file containing grouped snp gene expression lookup results (determined by the centroid snp's expression for each group)
# Description: TODO
# Run Time: TODO

#!/usr/bin/env python
import sys

groupedSnpTypesFilename = sys.argv[1]
groupedSnpTypesFile = open(groupedSnpTypesFilename, 'r')

# create dictionary containing snp type for each group of snps --> key = snp group, value = vector containing snp group, snp type, and snp source
snpTypeDict = {}

# ignore header line in groupedSnpTypesFilename
groupedSnpTypesFile.readline()

lineNum = 0
for line in groupedSnpTypesFile:
	lineNum += 1
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	snp = columns[0]
	snpType = columns[1] # "index" or "control"

	# determine snp type (lipid or T2D-like, testing or training) --> TODO: CHECK THAT THIS WAS DONE CORRECTLY!
	if lineNum in range(2, 1682):
		snpSource = "lipid_testing"
	elif lineNum in range(1682, 5026):
		snpSource = "lipid_training"
	elif lineNum in range(5026, 6066):
		snpSource = "T2Dlike_testing"
	elif lineNum in range(6066, 8098):
		snpSource = "T2Dlike_training"

	info = []
	info.append(snp)
	info.append(snpType)
	info.append(snpSource)

	snpTypeDict[snp] = info

groupedSnpTypesFile.close()

expressionVectorTableFilename = sys.argv[2]
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

# create counter variables for different types of snps
numIndexSnps = 0
numControlSnps = 0
numSnpsNoType = 0
numLipidTestingSnps = 0
numLipidTrainingSnps = 0
numT2DLikeTestingSnps = 0
numT2DLikeTrainingSnps = 0

# populate expression vector dictionary
for line in expressionVectorTableFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	snp = columns[0]

	# determine snp type from dictionary
	if snp not in snpTypeDict:
		print "ERROR: cannot determine snp type for", snp, "because this snp isn't in any of the lipid or T2D-like testing or training files." # TODO: Deal with this
		snpType = "none"
		numSnpsNoType += 1
	else:
		snpType = snpTypeDict[snp][1]
		if snpType == "index":
			numIndexSnps += 1
		else: # snpType == "control"
			numControlSnps += 1

	if snpTypeDict[snp][2] == "lipid_testing":
		numLipidTestingSnps += 1
	elif snpTypeDict[snp][2] == "lipid_training":
		numLipidTrainingSnps += 1
	elif snpTypeDict[snp][2] == "T2Dlike_testing":
		numT2DLikeTestingSnps += 1
	else: # snpSource == "T2Dlike_training"
		numT2DLikeTrainingSnps += 1
	
	vector = []
	for i in range(numTissues):
		vector.append(columns[i + 1])

	# add snp and type (index for all?? (TODO)) to beginning of vector
	vector.insert(0, snp)
	vector.insert(1, snpType)

	expressionVectorDictionary[snp] = vector

totalNumSnps = len(expressionVectorDictionary)
if (numIndexSnps + numControlSnps != totalNumSnps):
	print "ERROR: control + index != total"
if (numLipidTrainingSnps + numLipidTrainingSnps + numT2DLikeTrainingSnps + numT2DLikeTestingSnps != totalNumSnps):
	print "ERROR: lipidTesting + lipidTraining + T2Dtesting + T2Dtraining != total"

tab = "\t"
newline = "\n"

# create new header line
newHeaderLine = "snp" + tab + "type" + tab

# create output file
output = ""

for i in range(numSnps):
	key = expressionVectorDictionary.keys()[i]
	output += key + tab

	for j in range(numTissues):
		if j < (numTissues - 1):
			output += expressionVectorDictionary[key][j] + tab
		else: # create new line at the end of each vector
			output += expressionVectorDictionary[key][j] + newline

outFile.write(output)

expressionVectorTableFile.close()
outFile.close()