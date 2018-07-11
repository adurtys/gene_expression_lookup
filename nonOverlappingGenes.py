# Date Created: 11 July 2018
# Date Last Modified: 11 July 2018
# Execution: python nonOverlappingGenes.py
# Description:

#!/usr/bin/env python

# read in GTEx t-statistics file
tstatFilename = "GTEx.tstat.tsv"
tstatFile = open(tstatFilename, 'r')

# skip header line
headerLine = tstatFile.readline()

# parse GTEx file for GeneIds
idsInTstatFile = []
for line in tstatFile:
	line = line.rstrip('\r\n')
	tissues = line.split('\t')

	# add ENSGID, found in first column, to the list
	idsInTstatFile.append(tissues[0])

tstatFile.close()

# read in GTF file
gencodeFilename = "gencode.v19.annotation.gtf"
gencodeFile = open(gencodeFilename, 'r')

# skip first five lines
for i in range(5):
        gencodeFile.readline()

# parse GENCODE GTF file for GeneIds
idsInGencodeFile = []
for line in gencodeFile:
   	line = line.rstrip('\r\n')
	
	# split line on tab, creating list of possible columns
	columns = line.split('\t')
	
	# process otherInfo column by getting rid of key label and quotes
	otherInfo = columns[8]
	otherInfo = otherInfo.split(';')

	geneId = otherInfo[0]
	geneId = geneId.strip('gene_id ')
	geneId = geneId.strip('"')

	# add ENSGID to the list
	idsInGencodeFile.append(geneId)

gencodeFile.close()

newline = "\n"

# create output file for genes in GTEx file and not in GENCODE file
tstatGenesFilename = "genesNotInGencode.txt"
tstatGenesFile = open(tstatGenesFilename, 'w')

genesNotInGencode = []
tstatOutput = ""
for item in idsInTstatFile:
	if item not in idsInGencodeFile:
		genesNotInGencode.append(item)
		tstatOutput += item + newline

tstatGenesFile.write(tstatOutput)
tstatGenesFile.close()

# create output file for genes in GENCODE file and not in GTEx file
gencodeGenesFilename = "genesNotInGTEx.txt"
gencodeGenesFile = open(gencodeFilename, 'w')

genesNotInGTEx = []
gencodeOutput = ""
for item in idsInGencodeFile;
	if item not in idsInTstatFile:
		genesNotInGTEx.append(item)
		gencodeOutput += item + newline

gencodeGenesFile.write(gencodeOutput)
gencodeGenesFile.close()

print "There are", len(genesNotInGencode), "genes in the GTEx file but not in the GENCODE file."
print "There are", len(genesNotInGTEx), "genes in the GENCODE file but not in the GTEx file."