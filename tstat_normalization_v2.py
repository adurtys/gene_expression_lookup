# Date Created: 3 July 2018
# Date Last Modified: 4 July 2018
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

# Input: a list, an item to be searched within the list (optional)
# Output: If item is not specified, returns dictionary containing each item and its respective number of occurrences in the list.
# If item is specified, returns an integer for the number of occurrences of that item in the list. 
# Description: Determines the number of occurences for items in the list. If a particular item is specified, returns the number
# of occurrences of that item in the list. Otherwise, returns a dictionary containing the number of occurences of every item
# in the list. 
def numOccurrences(anyList, item = None):
	numOccurrences = {}

	for i in range(len(anyList)):
		if anyList[i] not in numOccurrences:
			numOccurrences[anyList[i]] = 1
		else: 
			# item is a duplicate
			numOccurrences[anyList[i]] += 1

	if item == None:
		# no item was specified to be parsed in the list
		return numOccurrences
	else:
		# parse list for the number of occurences of the specified item
		if item not in numOccurrences:
			return 0
	
	return numOccurrences[item]

# Input: a list
# Output: a rank-ordered version of the original list
# Description: Determines the ranks of every item in a list, then returns a rank-ordered list, replacing each item in the
# original list with its rank.
def rankList(anyList):
	# sort the inputted list
	sortedList = sorted(anyList)

	# create dictionary to store each item in the list and its corresponding rank
	ranksDict = {}

	# determine whether there are any duplicated items
	duplicates = findDuplicates(anyList)
	if len(duplicates) == 0:
		# list has no duplicated items
		# ranks are simply the indexes of the items in the sorted list
		for item in anyList:
			ranksDict[item] = sortedList.index(item)
	else:
		# list has duplicated items
		for item in anyList:
			# the ranks of unique items are their indexes in the sorted list
			if numOccurrences(anyList, item) == 1:
				ranksDict[item] = sortedList.index(item)
			else:
				# the rank of non-unique items is the average of their indexes in the sorted list
				numItem = numOccurrences(anyList, item)

				firstIndex = float(sortedList.index(item))

				# determine the average of the indexes in the sorted list
				sumRanksSoFar = firstIndex
				for i in range(1, numItem):
					sumRanksSoFar += firstIndex + i

				averageRank = sumRanksSoFar / numItem

				ranksDict[item] = averageRank

	# create the list replacing the items in the original list with their ranks
	rankedList = []
	for item in anyList:
		rank = ranksDict[item]
		rankedList.append(rank)

	return rankedList


# read in the tstat file
# strip file of new lines
# split based on tabs
# first column contains ensgId
# create matrix to store all info
# for each line in the file, read in each column, starting from the second, and fill out matrix
# rank-order each column



# read in the tstat file
inFilename = "GTEx.tstat.tsv"
inFile = open(inFilename, 'r')

# # create new file containing rank-normalized tstats
# outFilename = "normalizedGTEx.tstat.txt"
# outFile = open(outFilename, 'w')

# tab = "\t"
# newline = "\n"

# store column labels in header row
headerLine = inFile.readline()
headerLine = headerLine.rstrip('\r\n')
headers = headerLine.split('\t')

numCol = len(headers)

# initialize matrix with numColumns
matrix = [[0 for x in range(numCol)]]

# initialize another matrix of the same size to store the ranks
rankedMatrix = [[0 for x in range(numCol)]]

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