#! /usr/bin/env python

#This script makes new sets of alignments according to the output of alstats.py
#with some sequences and individuals excluded, and has the option of only including
#one or two sequences per individual.  When it does not use all of the sequences in the list,
#it chooses them randomly, so that the script can be run multiple times to create
#multiple alignments with a randomly-chosen subset of individuals.
#newal.py version 1.0 15 Jan. 2015 Abby Moore
'''
Example command:
newal.py ~/Documents/finalseqs1/ aa .fa ~/Documents/finalseqs1output/seqs_taxa_best.txt 2 ~/Documents/finalseqs2/ 2
newal.py ~/Documents/finalseqs1/ aa .fa ~/Documents/finalseqs1output_3LLC/seqs_pops_best.txt 2 ~/Documents/finalseqs2_3LLC/pops3LLC/ 3LLC
newal.py InFolder InFilePre InFilePost AlListIn SeqsperInd OutFolder OutFilePre
'''

import sys #We want to be able to send error messages to the screen.
from collections import defaultdict #We want to be able to make dictionaries with multiple levels.
import random #We want to be able to generate (pseudo-)random numbers.

Usage = '''
newal.py version 1.0 reads the files created by alstats.py and makes new, pruned
alignment files, with a random subset of individuals if all individuals in the alstats.py
output file are not to be included.
Input should be in the following format:
newal.py [folder containing alignment files] [prefix for old alignment files] [suffix
for old alignment files] [output file from alstats.py that should be used to make the new
alignments] [number of sequences per set (taxon, population, or individual, depending
on the alstats.py output file selected) to be used in the alignment--all if all are to be used]
[folder into which new alignment files should be put] [prefix for new alignments--none if no prefix] 
'''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 8:
	sys.exit("ERROR!  newal.py requires 7 additional arguments and you supplied %d.\n%s" % (len(sys.argv)-1, Usage))
InFolder = sys.argv[1]
InFilePre = sys.argv[2]
InFilePost = sys.argv[3]
AlListIn = sys.argv[4]
SeqsperInd = sys.argv[5]
OutFolder = sys.argv[6]
OutFilePre = sys.argv[7]

#making sure that the folder names end with a slash
if InFolder[-1] != "/":
	InFolder += "/"
print ("Sequence files will be read from the folder %s.\n" % (InFolder))
sys.stderr.write("Sequence files will be read from the folder %s.\n" % (InFolder))
if OutFolder[-1] != "/":
	OutFolder += "/"

#Setting up the dictionaries and things we will fill:
OTUDict = { } #The dictionary that says which column goes with which OTU (individual/
#population/taxon
OTUList = [ ] #The list of OTUs, so they can always be in the same order.
SeqNameDict = defaultdict(dict) #The dictionary of the sequence names that we want to use for each locus
SeqDict = { } #The actual dictionary of the sequences we want to use, with keys being locus names
#and values being the list of sequences for each locus.

InFile = open(AlListIn, 'rU')
LineNum = 0
for Line in InFile:
	Line = Line.strip("\n").strip("\r").split("\t")
	#If this is the first line, then read the column headings
	if LineNum == 0:
		OTUNum = 0
		for OTU in Line:
			if OTUNum != 0:
				OTUDict[OTUNum] = OTU
				OTUList.append(OTU)
			OTUNum += 1
	#Otherwise, read the sequence names
	else:
		Locus = Line[0]
		OTUNum = 0
		for OTUEntry in Line:
			if OTUNum != 0:
				OTUSeqs = OTUEntry.split(",")
				if SeqsperInd != 'all':
					#Making SeqsperInd an int again, because somehow it was not.
					SeqsperInd = int(SeqsperInd)
					if SeqsperInd < len(OTUSeqs):
						OTUSeqsn = random.sample(OTUSeqs, SeqsperInd)
						SeqNameDict[Locus][OTUDict[OTUNum]] = OTUSeqsn
					else:
						SeqNameDict[Locus][OTUDict[OTUNum]] = OTUSeqs
				else:
					SeqNameDict[Locus][OTUDict[OTUNum]] = OTUSeqs
			OTUNum += 1
	LineNum += 1
InFile.close()

print("Information about %d loci in %d OTUs (individuals, populations, or taxa) was read from the infile %s.\n" \
% (LineNum-1, OTUNum-1, AlListIn))
sys.stderr.write("Information about %d loci in %d OTUs (individuals, populations, or taxa) was read from the infile %s.\n" \
% (LineNum-1, OTUNum-1, AlListIn))

#random.sample(population, k)

#Now to read the alignment files for the individual loci and to find the sequences we want.
for Locus in SeqNameDict.keys():
	#First, reading the files:
	InFileName = InFolder + InFilePre + Locus + InFilePost
	InFile = open(InFileName, 'rU')
	SeqDictTemp = { }
	#And making a dictionary of all of the sequences in that file:
	for Line in InFile:
		Line = Line.strip('\n').strip('\r')
		if Line[0] == ">":
			SeqName = Line[1:]
		else:
			SeqDictTemp[SeqName] = Line
	InFile.close()
	#Now to make a list of the sequences we want to keep:
	SeqListTemp = [ ]
	for OTU in OTUList:
		SeqList = SeqNameDict[Locus][OTU]
		if SeqList != ['']:
			for SeqName in SeqNameDict[Locus][OTU]:
				SeqListTemp.append(">"+SeqName+"\n"+SeqDictTemp[SeqName]+"\n")
	#And save that list as an entry in a dictionary (key: locus name, value: set of sequences)
	SeqDict[Locus] = SeqListTemp

print("Sequence alignments for %d loci have been read.\n" % (len(SeqDict)))
sys.stderr.write("Sequence alignments for %d loci have been read.\n" % (len(SeqDict)))

#Now to write the lists of sequences to their files:
for Locus in SeqDict.keys():
	OutFileName = OutFolder + OutFilePre + Locus + ".fa"
	OutFile = open(OutFileName, 'w')
	for Seq in SeqDict[Locus]:
		OutFile.write(Seq)
	OutFile.close()

print("The pruned sequence alignments for these loci were written to files with names of the form %s.\n" % \
(OutFileName))
sys.stderr.write("The pruned sequence alignments for these loci were written to files with names of the form %s.\n" % \
(OutFileName))
