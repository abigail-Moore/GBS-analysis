#! /usr/bin/env python

#seq_group_concat.py version 1.0 31 May 2014 Abby Moore

#This script looks through the SeqGroups files produced by blast_parse_output.py 
#to see how many groups a given template sequence belongs to.  Perhaps some of these
#template sequences are too similar to other template sequences and the results
#would be better if they were removed.

import sys #We want to be able to get information from the command line
import re #We will need to be able to use regular expressions to parse the file names
from collections import defaultdict #We want to be able to make dictionaries with multiple levels

Usage = '''
seq_group_concat.py version 1.0
This script parses the SeqGroups files from blast_parse_output.py to look for
sequences that are present an unusually high number of times.
seq_group_concat.py [outfile name for list of groups] [fasta file with template sequences]
[name of new fasta file with sequences with many matches excluded]
[list of input files (SeqGroups files)]
'''

#format of SeqGroups files:
'''
1012_0 [0]: one sequence in pair
365_3 [1]: the other sequence in pair
3.0 [2]: the number of times this pair occurs (with each pair listed twice)
'''

if len(sys.argv) < 5:
	print(Usage)
else:
	GroupsOutFileName = sys.argv[1]
	TemplateSeqFile = sys.argv[2]
	SeqsOutFileName = sys.argv[3]
	InFileList = sys.argv[4:]

TDict = defaultdict(dict) #the keys will be the template sequences and the values will be the number of times
#they are present total and in each population
ManyMatches = { }
PopList = [ ] #List of the populations
SeqDict = { }
FileRe1 = r".*/(\w+)_SeqGroups\.txt"
FileRe2 = r"(\w+)_SeqGroups\.txt"
FileSub = r"\1"

#Reading the files with the lists of matching sequences
for InFileName in InFileList:
	#This takes care of the case when the whole file path is given
	PopName = re.sub(FileRe1, FileSub, InFileName)
	#This takes care of the case when only the file name is given
	if InFileName == PopName:
		PopName = re.sub(FileRe2, FileSub, InFileName)
	PopList.append(PopName)
	InFile = open(InFileName, 'rU')
	LineNum = 1
	for Line in InFile:
		if LineNum > 1:
			Line = Line.strip('\r').strip('\n').split('\t')
			#We only need to look at the first sequence, because each pair is listed twice.
			TSeqName = Line[0]
			NumMatches = int(float(Line[2]))
			try:
				TDict[TSeqName]['total'] += NumMatches
			except KeyError:
				TDict[TSeqName]['total'] = NumMatches
			try:
				TDict[TSeqName][PopName] += NumMatches
			except KeyError:
				TDict[TSeqName][PopName] = NumMatches
			LineNum += 1
		else:
			LineNum += 1
	InFile.close()
	print("%d lines were read from the file %s.\n" % (LineNum, InFileName))
	sys.stderr.write("%d lines were read from the file %s.\n" % (LineNum, InFileName))

#sorting the list of template sequences
TSeqList = TDict.keys()
TSeqList.sort()
#sorting the list of populations
PopList.sort()
#writing the list that will be outputted
ToWrite = [ ]
Line = "Seq Name\tTotal"
for PopName in PopList:
	Line += "\t"+PopName
Line += "\n"
ToWrite.append(Line)
for TSeqName in TSeqList:
	Line = TSeqName+"\t"+str(TDict[TSeqName]['total'])
	if TDict[TSeqName]['total'] > 99:
		ManyMatches[TSeqName] = TDict[TSeqName]['total']
	for PopName in PopList:
		try:
			Line += "\t"+str(TDict[TSeqName][PopName])
		except KeyError:
			Line += "\t0"
	Line += "\n"
	ToWrite.append(Line)
#writing the list to the output file
OutFile = open(GroupsOutFileName, 'w')
for Line in ToWrite:
	OutFile.write(Line)
OutFile.close()

print("Information about %d sequences was written to the file %s.\n" % (len(TSeqList), GroupsOutFileName))
sys.stderr.write("Information about %d sequences was written to the file %s.\n" % (len(TSeqList), GroupsOutFileName))

print("%d sequences had 100 or more matches.\nThey will now be removed from the database.\n" % (len(ManyMatches)))
sys.stderr.write("%d sequences had 100 or more matches.\nThey will now be removed from the database.\n" % (len(ManyMatches)))

#Reading the old template sequences from their fasta file.
InFile = open(TemplateSeqFile, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r')
	#Each sequence is two lines, first the name, with a > symbol first,
	if Line[0] == ">":
		SeqName = Line[1:]
	#then the sequence itself.
	else:
		Seq = Line
		SeqDict[SeqName] = Seq
InFile.close()

print("%d template sequences were read from the file %s.\n" % (len(SeqDict), TemplateSeqFile))
sys.stderr.write("%d template sequences were read from the file %s.\n" % (len(SeqDict), TemplateSeqFile))

#writing a new file of template sequences without the ones that have more than 100 matches.
OutList = [ ]
for SeqName in SeqDict:
	if SeqName in ManyMatches.keys():
		"skip"
	else:
		Line = ">"+SeqName+"\n"+SeqDict[SeqName]+"\n"
		OutList.append(Line)
OutFile = open(SeqsOutFileName, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()

print("The template sequences with fewer than 100 matches were written to the new\n \
file %s.  A new BLAST database needs to be made from this file, \nand the files of \
uniqued sequences need to be BLASTed against this once again.\n" % (SeqsOutFileName))
sys.stderr.write("The template sequences with fewer than 100 matches were written to the new\n \
file %s.  A new BLAST database needs to be made from this file, \n and the files of \
uniqued sequences need to be BLASTed against this once again.\n" % (SeqsOutFileName))
