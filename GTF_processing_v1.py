# Date Created: 27 June 2018
# Date Last Modified: 27 June 2018
# Execution: 
# This program proceses GENCODE Comprehensive Gene Annotation GTF File for the start and end locations of protein-coding genes. 
# The output is a file containing locations of protein-coding genes, and maintains both GeneID (ENSG#) and name. 

#!/usr/bin/env python

# read in the GTF file
	# file is tab-separated, header = True

# parse file for protein-coding genes

# write new file with four columns (GeneID, Gene Name, Start, End)

# for each line: 
	# determine whether gene already has start and end specified in new file
	# check both GeneID and Gene Name
		# if no:
			# write new line in new file for gene
		# if yes:
			# check start position
			# if GTF_start < newFile_start, then newFile_start == GTF_start
			
			# check end position
			# if GTF_end > newFile_end, then newFile_end == GTF_end