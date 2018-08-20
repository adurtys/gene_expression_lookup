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

	snpSource = ""
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

# store column lables in header row
headerLine = expressionVectorTableFile.readline()
headerLine = headerLine.rstrip('\r\n')
headers = headerLine.split('\t')

# exclude first column when counting the number of tissues present in the file
numTissues = len(headers) - 1

# create data structures for the different snp types --> key = snp, value = expression vector
numIndexSnps = 0
numControlSnps = 0
lipidTestingSnps = {}
lipidTrainingSnps = {}
T2DLikeTestingSnps = {}
T2DLikeTrainingSnps = {}

# populate expression vector dictionary
for line in expressionVectorTableFile:
	line = line.rstrip('\r\n')
	columns = line.split('\t')

	snp = columns[0]

	# determine snp type from dictionary --> discard snps not in snpTypeDict
	if snp in snpTypeDict:
		vector = []
		for i in range(numTissues):
			vector.append(columns[i + 1])

		snpType = snpTypeDict[snp][1]
		if snpType == "index":
			numIndexSnps += 1
		else: # snpType == "control"
			numControlSnps += 1

		# add snp and type to beginning of vector
		vector.insert(0, snp)
		vector.insert(1, snpType)

		if snpTypeDict[snp][2] == "lipid_testing":
			lipidTestingSnps[snp] = vector
		elif snpTypeDict[snp][2] == "lipid_training":
			lipidTrainingSnps[snp] = vector
		elif snpTypeDict[snp][2] == "T2Dlike_testing":
			T2DLikeTestingSnps[snp] = vector
		else: # snpSource == "T2Dlike_training"
			T2DLikeTrainingSnps[snp] = vector

# create output files
lipidTestingOutFilename = "lipid_testing_expression_ML_table.txt"
lipidTrainingOutFilename = "lipid_training_expression_ML_table.txt"
T2DLikeTestingOutFilename = "T2D_like_testing_expression_ML_table.txt"
T2DLikeTrainingOutFilename = "T2D_like_training_expression_ML_table.txt"

# open the four output files
lipidTestingOutFile = open(lipidTestingOutFilename, 'w')
lipidTrainingOutFile = open(lipidTrainingOutFilename, 'w')
T2DLikeTestingOutFile = open(T2DLikeTestingOutFilename, 'w')
T2DLikeTrainingOutFile = open(T2DLikeTrainingOutFilename, 'w')

tab = "\t"
newline = "\n"

# create new header line
newHeaderLine = "snp" + tab + "type" + tab
for i in range(numTissues):
	if i < (numTissues - 1):
		newHeaderLine += headers[i + 1] + tab
	else: # create new line at end
		newHeaderLine += headers[i + 1] + newline

# create lipid testing output file
lipidTestingOutput = newHeaderLine
for snp in lipidTestingSnps:
	snp = lipidTestingSnps[snp][0]
	snpType = lipidTestingSnps[snp][1]

	lipidTestingOutput += snp + tab + snpType + tab

	for i in range(numTissues):
		if i < (numTissues - 1):
			lipidTestingOutput += lipidTestingSnps[i + 2] + tab
		else: # create new line at the end of each vector
			lipidTestingOutput += lipidTestingSnps[i + 2] + newline

lipidTestingOutFile.write(lipidTestingOutput)
lipidTestingOutFile.close()

# create lipid training output file
lipidTrainingOutput = newHeaderLine
for snp in lipidTrainingSnps:
	snp = lipidTrainingSnps[snp][0]
	snpType = lipidTrainingSnps[snp][1]

	lipidTrainingOutput += snp + tab + snpType + tab

	for i in range(numTissues):
		if i < (numTissues - 1):
			lipidTrainingOutput += lipidTrainingSnps[i + 2] + tab
		else: # create new line at the end of each vector
			lipidTrainingOutput += lipidTrainingSnps[i + 2] + newline

lipidTrainingOutFile.write(lipidTrainingOutput)
lipidTrainingOutFile.close()

# create T2D-like testing output file
T2DLikeTestingOutput = newHeaderLine
for snp in T2DLikeTestingSnps:
	snp = T2DLikeTestingSnps[snp][0]
	snpType = T2DLikeTestingSnps[snp][1]

	T2DLikeTestingOutput += snp + tab + snpType + tab

	for i in range(numTissues):
		if i < (numTissues - 1):
			T2DLikeTestingOutput += T2DLikeTestingSnps[i + 2] + tab
		else: # create new line at the end of each vector
			T2DLikeTestingOutput += T2DLikeTestingSnps[i + 2] + newline

T2DLikeTestingOutFile.write(T2DLikeTestingOutput)
T2DLikeTestingOutFile.close()

# create T2D-like training output file
T2DLikeTrainingOutput = newHeaderLine
for snp in T2DLikeTrainingSnps:
	snp = T2DLikeTrainingSnps[snp][0]
	snpType = T2DLikeTrainingSnps[snp][1]

	T2DLikeTrainingOutput += snp + tab + snpType + tab

	for i in range(numTissues):
		if i < (numTissues - 1):
			T2DLikeTrainingOutput += T2DLikeTrainingSnps[i + 2] + tab
		else: # create new line at the end of each vector
			T2DLikeTrainingOutput += T2DLikeTrainingSnps[i + 2] + newline

T2DLikeTrainingOutFile.write(T2DLikeTrainingOutput)
T2DLikeTrainingOutFile.close()

expressionVectorTableFile.close()