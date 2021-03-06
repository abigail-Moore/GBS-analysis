#! /usr/bin/env python

#blast_parse_output.py version 2 22 Aug. 2014 Abby Moore
#This script is for analyzing the output from BLASTing the full set
#of unique sequences from a couple of individuals against the database of cleaned
#template sequences.  The results of this analysis help to remove chimeric sequences
#and others that have too many hits.

'''
The columns in the BLAST tabular output (-outfmt '6 std qlen slen'):
[0] qseqid
[1] sseqid
[2] pident
[3] length (meaning length of aligned portion)
[4] mismatch (number of mismatches)
[5] gapopen (number of gap openings)
[6] qstart
[7] qend
[8] sstart
[9] send
[10] evalue
[11] bitscore
[12] qlen
[13] slen
'''

Usage = '''
blast_parse_output.py version 2.0 calculates statistics from the output of running
a BLAST search of all sequences against the cleaned template sequences
Input is of the form:
blast_parse_output.py [output prefix--for statistics/lists, not for sequences]
[sequence output folder] [BLAST output filename] [sequence file name]
'''

from collections import defaultdict #We want to be able to make dictionaries with multiple levels.
import sys #We want to be able to output error messages to the screen.
from Bio.Seq import Seq #We need to be able to reverse-complement sequences.
from Bio.Alphabet import IUPAC

if len(sys.argv) < 5:
	print Usage
else:
	OutFilePre = sys.argv[1]
	SeqFolder = sys.argv[2]
	BLFileName = sys.argv[3]
	SeqFileName = sys.argv[4]

print("blast_parse_output.py")

#the dictionaries we will fill out:
SeqDict = { } #dictionary for the sequences
SeqInfo = defaultdict(list) #This is the information sorted by query sequence.
#This dictionary is for determining whether multiple template sequences that match the same
#query sequence always form discrete (non-overlapping) groups, or whether the groups overlap.
#But the groups do appear to overlap, so it is not that interesting, and should perhaps be removed
#from future versions.
SeqGroups = defaultdict(dict)
#The keys of the following dictionary are the number of times a query sequence is matched,
#while the values are the number of times a query sequence with that number of matches is found.
NumMatchesDict = { }
FRDict = { } #dictionary that tells whether each sequence needs to be reverse-complemented
#before being added to the alignment

GroupsDict = { } #a dictionary where the keys are groups and the values are lists of 
#sequences in that group (only the ones that belong to that group only)
AmbigDict = { } #dictionary where the keys are sequences that belong to multiple groups and the values are lists of
#groups in which that sequence is present

#checking to make sure SeqFolder ends with a /
if SeqFolder[-1] != "/":
	SeqFolder = SeqFolder+"/"

#reading the sequence file
InFile = open(SeqFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r')
	#Each sequence is two lines, first the name, with a > symbol first,
	if Line[0] == ">":
		SeqName = Line[1:]
	#then the sequence itself.
	else:
		Sequ = Line
		SeqDict[SeqName] = Sequ
InFile.close()
print ("%d sequences were read from the file %s.\n" % (len(SeqDict), SeqFileName))
sys.stderr.write("%d sequences were read from the file %s.\n" % (len(SeqDict), SeqFileName))

#reading the BLAST output and saving the information into a dictionary where the keys
#are the query sequence names and the values are lists of matching template sequences
InFile = open(BLFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	SeqName = Line[0]
	Match = Line[1]
	SimPercent = float(Line[2])
	AlLength = float(Line[3])
	if int(Line[12]) <= int(Line[13]):
		SeqLength = float(Line[12])
	else:
		SeqLength = float(Line[13])
	AlPercent = AlLength/SeqLength
	if (SimPercent >= 90) and (AlPercent >= 0.90):
		try:
			SeqInfo[SeqName].append(Match)
		except KeyError:
			SeqInfo[SeqName] = [Match]
		if int(Line[8]) < int(Line[9]):
			FRDict[SeqName] = "F"
		else:
			FRDict[SeqName] = "R"
InFile.close()

print("%d of these sequences have matches in the file %s.\n" % (len(SeqInfo), BLFileName))
sys.stderr.write("%d of these sequences have matches in the file %s.\n" % (len(SeqInfo), BLFileName))


#Now the lists need to be made into sets, so that each template sequence is present only once.
for SeqName in SeqInfo:
	if len(SeqInfo[SeqName]) > 1:
		ListTemp = list(set(SeqInfo[SeqName]))
		SeqInfo[SeqName] = ListTemp
		#The template sequences need to be added to the dictionary of
		#template sequences with matches.
		for FirstSeq in ListTemp:
			for SecondSeq in ListTemp:
				if FirstSeq != SecondSeq:
					try:
						SeqGroups[FirstSeq][SecondSeq] += 0.5
					except KeyError:
						SeqGroups[FirstSeq][SecondSeq] = 0.5
					try:
						SeqGroups[SecondSeq][FirstSeq] += 0.5
					except KeyError:
						SeqGroups[SecondSeq][FirstSeq] = 0.5
	#Statistics on the number of matches also need to be collected.
	NumMatches = len(SeqInfo[SeqName])
	try:
		NumMatchesDict[NumMatches] += 1
	except KeyError:
		NumMatchesDict[NumMatches] = 1
		
for SeqName in SeqInfo:
	#Classifying the sequences with one match by alignment (template sequence) name
	if len(SeqInfo[SeqName]) == 1:
		try:
			GroupsDict[SeqInfo[SeqName][0]].append(SeqName)
		except KeyError:
			GroupsDict[SeqInfo[SeqName][0]] = [SeqName]
	#adding the sequences with multiple matches to a dictionary to be dealt with later
	elif len(SeqInfo[SeqName]) > 1:
		AmbigDict[SeqName] = SeqInfo[SeqName]
print("Of these, %d have multiple matches and need to be examined separately.\n" % (len(AmbigDict)))
sys.stderr.write("Of these, %d have multiple matches and need to be examined separately.\n" % (len(AmbigDict)))

#Now to write everything to files:
ToWrite = [ ]
Line = "Seq_Name\tNum_90_90\tSeqs_90_90\n"
ToWrite.append(Line)
SortedSeqs = SeqInfo.keys()
SortedSeqs.sort()
for SeqName in SortedSeqs:
	Line = SeqName+"\t"+str(len(SeqInfo[SeqName]))+"\t"+', '.join(SeqInfo[SeqName])+"\n"
	ToWrite.append(Line)
OutFileName = OutFilePre+"_SeqStats.txt"
OutFile = open(OutFileName, 'w')
for Line in ToWrite:
	OutFile.write(Line)
OutFile.close()

print("The list of sequences and their matching templates was written to the file %s.\n" \
	% (OutFileName))
sys.stderr.write("The list of sequences and their matching templates was written to the file %s.\n" \
	% (OutFileName))

ToWrite = [ ]
Line = "Num_Matches\tNum_Seqs\n"
ToWrite.append(Line)
OutFileName = OutFilePre+"_NumMatchesDict.txt"
for Num in NumMatchesDict:
	Line = str(Num)+"\t"+str(NumMatchesDict[Num])+"\n"
	ToWrite.append(Line)
OutFile = open(OutFileName, 'w')
for Line in ToWrite:
	OutFile.write(Line)
OutFile.close()

print("Information about how many sequences had a given number of matches was written to \
the file %s.\n" % (OutFileName))
sys.stderr.write("Information about how many sequences had a given number of matches was written to \
the file %s.\n" % (OutFileName))

ToWrite = [ ]
OutFileName = OutFilePre+"_SeqGroups.txt"
Line = "Seq\tMatch\tNum_Matches\tAll sequences should be present twice.\n"
ToWrite.append(Line)
SortedFirstSeqs = SeqGroups.keys()
SortedFirstSeqs.sort()
for FirstSeq in SortedFirstSeqs:
	SortedSecondSeqs = SeqGroups[FirstSeq].keys()
	SortedSecondSeqs.sort()
	for SecondSeq in SortedSecondSeqs:
		Line = FirstSeq+"\t"+SecondSeq+"\t"+str(SeqGroups[FirstSeq][SecondSeq])+"\n"
		ToWrite.append(Line)
OutFile = open(OutFileName, 'w')
for Line in ToWrite:
	OutFile.write(Line)
OutFile.close()

print("The list of pairs of templates that both matched the same query sequence was written\
 to the file %s.\n" % (OutFileName))
sys.stderr.write("The list of pairs of templates that both matched the same query sequence was written\
 to the file %s.\n" % (OutFileName))

#Now to go through the list of sequences and print them to alignments.
#This is a bit excessive, but this way multiple BLAST files can be read, one after
#another.
#It seems me to that there is no reason for sequences with only one match to be
#written to files like this, since the whole process will be repeated with
#the final set of sequences anyway, so I am commenting this out.
'''
for GroupName in GroupsDict:
	OutFileName = SeqFolder+GroupName+".fa"
	OutFile = open(OutFileName, 'a')
	for SeqName in GroupsDict[GroupName]:
		if FRDict[SeqName] == "R":
			SeqF = SeqDict[SeqName]
			SeqF = Seq(SeqF, IUPAC.unambiguous_dna)
			SeqRC = str(SeqF.reverse_complement())
			SeqDict[SeqName] = SeqRC
		Line = ">"+SeqName+"\n"+SeqDict[SeqName]+"\n"
		OutFile.write(Line)
	OutFile.close()

print("The sequences that have only one match have been written to their groups' files,\
 which have names of the form %s.\n" % (OutFileName))
sys.stderr.write("The sequences that have only one match have been written to their groups' files,\
 which have names of the form %s.\n" % (OutFileName))
'''

#Now to make a list of the sequences with multiple matches:
ToWrite = [ ]
OutFileName = OutFilePre+"_AmbigSeqs.txt"
Line = "Seq\tMatching Seqs\n"
ToWrite.append(Line)
for SeqName in AmbigDict:
	Line = SeqName+"\t"+", ".join(AmbigDict[SeqName])+"\n"
	ToWrite.append(Line)
OutFile = open(OutFileName, 'w')	
for Line in ToWrite:
	OutFile.write(Line)
OutFile.close()

print("The list of ambiguous sequences and their matches was written to the file %s.\n" \
	% (OutFileName))
sys.stderr.write("The list of ambiguous sequences and their matches was written to the file %s.\n" \
	% (OutFileName))
