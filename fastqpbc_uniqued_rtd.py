#! /usr/bin/env python

''' "Uniqued" file format (with sequence and Phred score trimmed to 5 bp):

GATCC	9	EEHFG	Mlk38_20,Mlk38_19	7,2	.	.

GATCC	44	DDFGF	Mlk38_20	44	.	.

There does not appear to be anything special at the end of the file.'''

#This program reads the files of with the sequences sorted by individual by fastqcomb_pbc.py.
#It makes two dictionaries of the unique sequences with information about how often and in which individuals they occur.
#These dictionaries are written to files of "uniqued" sequences, in the correct format for the rtd pipeline.
#First, it makes a dictionary that only includes the first 30 bp of the sequences for analysis by the rtd pipeline.
#Second, it makes a dictionary  that includes the entire sequences, for the second set of anlayses.

import re #We will need the regular expressions module
import sys #We want to be able to output error messages to the screen.
import gzip #We want to be able to open zipped files.

from collections import defaultdict #We need defaultdict to be able to create dictionaries with multiple levels.

#dezip is from preprocess_radtag_lane.py from the rtd program
def dezip(values_in):
    '''opposite of zip(), i.e. 
    >>> dezip([("a",1),("b",2),("c",3)])
    (["a","b","c"],[1,2,3])'''

    if isinstance(values_in,tuple):
        values_in = [values_in]
	
    lol = []
    for i in range(len(values_in[0])):
        lol.append([])
    for l in values_in:
        for it,li in zip(l,lol):
	    li.append(it)
    return tuple(lol)

#First, get the files names from the command line or print some instructions if not enough names are found..

Usage = '''
fastqpbc_uniqued_rtd.py - version 1.0, 9 March 2014
Reads separate files of the sequences present in each individual (output from fastqcomb_pbc.py)
and combines them into an input file for the rtd pipeline.
Usage:
	fastqpbc_uniqued_rtd.py (length of sequences in uniqued file) 
	(output file name--without uniqued.gz) (sequence file names)
The file names can be either the full path, e.g., /home/moorea/... 
or the file name only if you are in that directory.

Sequence file names can either be a list of all sequence file names,
or something of the form "outfile_*.fastq.gz".
And the sequence files can either be be zipped (with a ".gz" ending) or not.
'''

if len(sys.argv) <3:
	print Usage
else:
	SeqLen = int(sys.argv[1])
	OutFileRoot = sys.argv[2]
	InFileList = sys.argv[3:]

print("fastqpbc_uniqued_rtd.py")

#These are the default values for testing.
#SeqLen = 30	
#InFileList = ["/home/abby/Documents/Moore1/230_1_all.fastq.gz","/home/abby/Documents/Moore1/230_3_all.fastq.gz","/home/abby/Documents/Moore1/230_5_all.fastq.gz","/home/abby/Documents/Moore1/230_6_all.fastq.gz","/home/abby/Documents/Moore1/230_7_all.fastq.gz"]
#OutFileRoot = "/home/abby/Documents/sandbox/outfile"

OutFileNameShort = OutFileRoot+str(SeqLen)+".uniqued.gz"
OutFileNameAll = OutFileRoot+"_all.uniqued.gz"

NumInFiles = 0
NumSeqsTot = 0
NumLimit = 1 #This can be increased if you want to only include sequences that are present more than a certain number of times.


#the dictionary with multiple levels to hold the uniqued sequences--first 30 bp
UniquedSeqsShort = defaultdict(dict) 
UniquedSeqsAll = defaultdict(dict) #dictionary for entire sequences
sys.stderr.write("The first %d bases of each sequence will be written to the file %s.\n" % (SeqLen, OutFileNameShort))
print("The first %d bases of each sequence will be written to the file %s.\n" % (SeqLen, OutFileNameShort))
sys.stderr.write("The entire sequences will be written to the file %s.\n" % (OutFileNameAll))
print("The entire sequences will be written to the file %s.\n" % (OutFileNameAll))

#%%%%%%%%Reading the files of sequences and putting them into the dictionary.%%%%%%%%%%%%%%
for InFileName in InFileList:
	if InFileName.endswith('.gz'):
		InFile = gzip.open(InFileName)
	else:
		InFile = open(InFileName,'rU')
	LineNum = 0
	NumSeqs = 0
	for Line in InFile:
		Line = Line.strip('\n').strip('\r')
		SeqLine = (LineNum + 3) % 4
		if LineNum == 0: #The first line of the file should have the individual name only.
			Ind = Line
			LineNum += 1
		elif SeqLine == 0 or SeqLine == 2: #The position line and the line with the barcode length and individual number are not copied
			LineNum += 1
		elif SeqLine == 1: #The first bases of the sequence line are saved
			SeqShort = Line[:SeqLen]
			SeqAll = Line
			LineNum += 1
			NumSeqs += 1
		elif SeqLine == 3: #The first bases of the quality score line are saved and potentially transformed
			#(although this is commented out).
			QualShort = Line[:SeqLen]
			QualAll = Line
#The transformation of the quality score is from preprocess_radtag_lane.py from the rtd program and allows one
#to make a composite quality score from the individual quality scores.  But this adds quite a bit of information
#to the dictionary and I do not want it to get too big, so I am leaving this out for now.
#If this is put back in, the output file will have to be changed to make the composite quality score and
#transform it back into a normal quality score
			#Qual = [ord(c)-33 for c in Line] #***33 is for Phred33 scores.  If you are using Phred64 scores, this must be changed.***
			LineNum += 1
			#first, the file of the full sequences:
			try: #First, I assume the sequence is already in the list, and just add that this individual to the sequence's entry.
				UniquedSeqsAll[SeqAll]['count'][Ind] += 1
			except KeyError: #If the sequence is not found, then I add an entry for that sequence
				UniquedSeqsAll[SeqAll]['count'] = defaultdict(int)
				UniquedSeqsAll[SeqAll]['count'][Ind] += 1
				UniquedSeqsAll[SeqAll]['quality'] = QualAll
			#then the file of the shortened sequences for the first analysis:
			try: #First, I assume the sequence is already in the list, and just add that this individual to the sequence's entry.
				UniquedSeqsShort[SeqShort]['count'][Ind] += 1
				#UniquedSeqs[Seq]['sum_quality'] += Qual #And potentially the quality score.
			except KeyError: #If the sequence is not found, then I add an entry for that sequence
				UniquedSeqsShort[SeqShort]['count'] = defaultdict(int)
				UniquedSeqsShort[SeqShort]['count'][Ind] += 1
				UniquedSeqsShort[SeqShort]['quality'] = QualShort
	sys.stderr.write("%s has been read.\n" % InFileName)
	print("%s has been read.\n" % InFileName)
	InFile.close()
	NumInFiles += 1
	NumSeqsTot += NumSeqs

sys.stderr.write("%d sequences were read from %d files.\n" % (NumSeqsTot, NumInFiles))
print("%d sequences were read from %d files.\n" % (NumSeqsTot, NumInFiles))

if len(InFileList) != NumInFiles:
	sys.stderr.write("Error: %d files were in the list of infiles, but only %d files were read.\n" % (len(InFileList), NumInFiles))
	print("Error: %d files were in the list of infiles, but only %d files were read.\n" % (len(InFileList), NumInFiles))

#%%%%%%%%%%%%Reformatting the dictionaries and writing them to files.%%%%%%%%%%%%%%%%%%%%%%%%
LineNum = 0
CommonSeqs = 0
FinalSeqs = [ ] #This will be the list of reformatted, uniqued sequences
NumUSeqsDict = {} #This dictionary stores information about the number of times sequences are present in a given number of individuals.
for SeqShort in UniquedSeqsShort.keys():
	#making the line to be written to the list
	NewSeq = SeqShort
	NewQual = UniquedSeqsShort[SeqShort]['quality']
	#The following three lines are modified from preprocess_radtags_lanes.py from the rtd pipeline.
	ind,indcount = dezip(sorted([(k,v) for k,v in UniquedSeqsShort[SeqShort]['count'].items()],reverse=True,key = lambda x:x[1]))
	SeqTotal = sum(indcount)
	NumInds = len(ind)
	#writing the line to the list
	#(assuming more than the lower limit of individuals have that sequence)
	if NumInds >= NumLimit:
		NewLine = '%s\t%s\t%s\t%s\t%s\t.\t.\n' % \
		(NewSeq,SeqTotal,NewQual, ','.join(ind), ','.join([str(i) for i in indcount]))
		FinalSeqs.append(NewLine)
		CommonSeqs += 1
	#Collecting statistics on how many individuals have each sequence
	try:
		NumUSeqsDict[NumInds] += 1
	except KeyError:
		NumUSeqsDict[NumInds] = 1

#These two lines can be added back in if only the sequences present in more than a given number of individuals will be included.		
#sys.stderr.write("%d sequences were found in at least %d individuals.\n" % (CommonSeqs, NumLimit))
#print("%d sequences were found in at least %d individuals.\n" % (CommonSeqs, NumLimit))

#Printing a list of how many individuals have each sequence

sys.stderr.write("number of individuals\tnumber of sequences present in that many individuals\n")
print("number of individuals\tnumber of sequences present in that many individuals\n")
for key in NumUSeqsDict:
	print("%d\t%d\n" % (key, NumUSeqsDict[key]))
	sys.stderr.write("%d\t%d\n" % (key, NumUSeqsDict[key]))

#The reformatted, uniqued sequences are written to a file.
OutFile = gzip.open(OutFileNameShort, 'w')
for Line in FinalSeqs:
	OutFile.write(Line)
OutFile.close()

	#Writing a file that the rtd pipeline may or may not need
StrNums = "%d\t%d" % (CommonSeqs, len(FinalSeqs))
CacheOutFileName = OutFileNameShort+".rc.cache"
OutFile = open(CacheOutFileName, 'w')
OutFile.write(StrNums)
OutFile.close()

sys.stderr.write("%d unique sequences of %d bases were found and written to %s.\n" % (CommonSeqs, SeqLen, OutFileNameShort))
sys.stderr.write("Statistics written to %s.\n" % (CacheOutFileName))
print("%d unique sequences of %d bases were found and written to %s.\n" % (CommonSeqs, SeqLen, OutFileNameShort))
print("Statistics written to %s.\n" % (CacheOutFileName))

#Now the entire sequences:
LineNum = 0
CommonSeqs = 0
FinalSeqs = [ ] #This will be the list of reformatted, uniqued sequences
NumUSeqsDict = {} #This dictionary stores information about the number of times sequences are present in a given number of individuals.
for SeqAll in UniquedSeqsAll.keys():
	#making the line to be written to the list
	NewSeq = SeqAll
	NewQual = UniquedSeqsAll[SeqAll]['quality']
	#The following three lines are modified from preprocess_radtags_lanes.py from the rtd pipeline.
	ind,indcount = dezip(sorted([(k,v) for k,v in UniquedSeqsAll[SeqAll]['count'].items()],reverse=True,key = lambda x:x[1]))
	SeqTotal = sum(indcount)
	NumInds = len(ind)
	#writing the line to the list
	#(assuming more than the lower limit of individuals have that sequence)
	if NumInds >= NumLimit:
		NewLine = '%s\t%s\t%s\t%s\t%s\t.\t.\n' % \
		(NewSeq,SeqTotal,NewQual, ','.join(ind), ','.join([str(i) for i in indcount]))
		FinalSeqs.append(NewLine)
		CommonSeqs += 1
	#Collecting statistics on how many individuals have each sequence
	try:
		NumUSeqsDict[NumInds] += 1
	except KeyError:
		NumUSeqsDict[NumInds] = 1

#These two lines can be added back in if only the sequences present in more than a given number of individuals will be included.		
#sys.stderr.write("%d sequences were found in at least %d individuals.\n" % (CommonSeqs, NumLimit))
#print("%d sequences were found in at least %d individuals.\n" % (CommonSeqs, NumLimit))

#Printing a list of how many individuals have each sequence
sys.stderr.write("number of individuals\tnumber of sequences present in that many individuals\n")
print("number of individuals\tnumber of sequences present in that many individuals\n")
for key in NumUSeqsDict:
	print("%d\t%d\n" % (key, NumUSeqsDict[key]))
	sys.stderr.write("%d\t%d\n" % (key, NumUSeqsDict[key]))

#The reformatted, uniqued sequences are written to a file.
OutFile = gzip.open(OutFileNameAll, 'w')
for Line in FinalSeqs:
	OutFile.write(Line)
OutFile.close()

#Writing a file that the rtd pipeline may or may not need
StrNums = "%d\t%d" % (CommonSeqs, len(FinalSeqs))
CacheOutFileName = OutFileNameAll+".rc.cache"
OutFile = open(CacheOutFileName, 'w')
OutFile.write(StrNums)
OutFile.close()

sys.stderr.write("%d unique full sequences were found and written to %s.\n" % (CommonSeqs, OutFileNameAll))
sys.stderr.write("Statistics written to %s.\n" % (CacheOutFileName))
sys.stderr.write("\n\n")
print("%d unique sequences were found and written to %s.\n" % (CommonSeqs, OutFileNameAll))
print("Statistics written to %s.\n" % (CacheOutFileName))
print("\n\n")
