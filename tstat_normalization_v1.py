# Date Created: 3 July 2018
# Date Last Modified: 3 July 2018
# Execution: python tstat_normalization_v1.py
# TODO: This program rank-normalizes the t statistics in the "GTEx.tstat.tsv" file

#!/usr/bin/env python

# Input: a list
# Output: a list of items that are present in the list more than once
# Description: finds items that are present more than once in the inputted list 
def findDuplicates(anyList):
	seen = set()
	duplicates = []

	for item in anyList:
		if item in seen and item not in duplicates:
			duplicates.append(item)
		else:
			seen.add(item)

	return duplicates

# Input: a list, and 
# Output: a list of the 
def numDuplicates(anyList, item):
	seen = set()
	duplicates = []
	numDupDict = {}
	for i in range(len(anyList)):
		if anyList[i] not in seen:
			seen.add(i)
		# item is a duplicate
		else:
			if anyList[i] not in duplicates:
				duplicates.append(anyList[i])
				numDupDict[anyList[i]] = 1
			else:
				numDupDict[anyList[i]] += 1

	if item in numDupDict:
		return numDupDict[item]
	return 0


def numOccurrences(anyList, item):
	numOccurrences = {}

	for i in range(len(anyList)):
		if anyList[i] not in numOccurrences:
			numOccurrences[anyList[i]] = 1
		else: 
			# item is a duplicate
			numOccurrences[anyList[i]] += 1

	if item not in numOccurrences:
		return 0
	return numOccurrences[item]

def numOccurrences(anyList):
	numOccurrences = {}

	for i in range(len(anyList)):
		if anyList[i] not in numOccurrences:
			numOccurrences[anyList[i]] = 1
		else:
			# item is a duplicate
			numOccurrences[anyList[i]] += 1

	return numOccurrences



def rankList(anyList):
	sortedList = sorted(anyList)
	duplicates = findDuplicates(anyList)
	
	# create a dictionary with numbers in the unsorted list as keys and ranks as values
	ranksDict = {}
	for item in anyList:
		# if no duplicated values, then rank is the index of each item in the sorted list
		if len(duplicates) == 0:
			ranksDict[item] = sortedList.index(item)
		# duplicated values are present in the list
		else:
			for number in duplicates:
				# if the number is duplicated, its rank is the average of its indexes in the sorted list
				if item == number:
					# determine the number of time each non-unique value occurs in the list
					numDup = numDuplicates(anyList, number)

					# determine the first index of the number
					firstIndex = sortedList.index(item)

					# calculate the average of the indexes
					sum = float(firstIndex)
					for i in range(numDup):
						sum += firstIndex + i + 1
					averageRank = sum / (numDup + 1)

					ranksDict[item] = averageRank

				else:
					# if the number is not duplicated, its rank is its index in the sorted list
					ranksDict[item] = sortedList.index(item)

	# create a list containing the ranks for the numbers in the unsorted list
	ranks = []
	for num in anyList:
		rank = ranksDict[num]
		ranks.append(rank)

	return ranks

# read in the tstat file
inFilename = "GTEx.tstat.tsv"
inFile = open(inFilename, 'r')

# create new file containing rank-normalized tstats
outFilename = "normalizedGTEx.tstat.txt"
outFile = open(outFilename, 'w')

tab = "\t"
newline = "\n"

# store column labels in header row
headerLine = inFile.readline()
headerLine = headerLine.rstrip('\r\n')
headers = headerLine.split('\t')

numColumns = len(headers)

# initialize matrix with numColumns
matrix = [[0 for x in range(numColumns)]]

# initialize a temporary column to store values to be rank-ordered
tempColumn = []

numRowsSoFar = 0
for line in inFile:
	line = line.rstrip('\r\n')
	tissues = line.split('\t')

	# keep first column as the ENSGID
	id = tissues[0]
	matrix[numRowsSoFar][0] = tissues[0]
	
	# for every subsequent column, rank order the column, then append it to the matrix


	numRowsSoFar += 1







prevOutput = ""
output = ""
# start at the second column (skip first column, containing ENSGID)
# FOR TESTING, numColumns = 3!!!!!! (also didn't rank yet)
for i in range(1, 3):
	print "i:", i
	tempColumn = []

	for line in output:
		prevOutput = output.rstrip('\r\n')

	for line in inFile:
		line = line.rstrip('\r\n')
		tissues = line.split('\t')

		tempColumn.append(tissues[i])

	print "tempColumn[1]:", tempColumn[1]
	for tstat in tempColumn:
		tstat = str(tstat)
		output = prevOutput + tstat + tab + newline
	
	# rankedColumn = rankList(tempColumn)
	# for rank in rankedColumn:
	# 	rank = str(rank)
	# 	output += rank + newline
	# 	prevOutput = output

outFile.write(output)


outFile.close()
inFile.close()