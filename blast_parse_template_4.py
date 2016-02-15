#! /usr/bin/env python

#blast_parse_template.py version 3.0 30 April 2014, Abby Moore
#This script is supposed to reads the output of the BLAST analysis that tries to align
#the template sequences from the rtd alignment against each other (to get rid of
#reverse-complement sequences and others that are too similar).  It is not called
#directly, but instead is called as part of the script written by bam30_uniq_split_5.py
#Version 2 differs from version 1 in that it considers sequences to match when they are more than 80%
#identical over more than 75 bp, instead of requiring them to match over their entire lengths.
#It also looks at a second file in which the sequences were blasted against primers
#and removes all sequences that match the primers (since they match over enough of their
#lengths that the remaining sequences would be very short in most cases).
#Now, this version allows sequences that are more than 90% identical over 40 bp.  Hopefully
#this should deal with some of the problem of multiple matching of shorter query sequences.
#There should still potentially be some sort of criteria to select which template sequence is
#chosen to go into the final database.
#I am trying to have a criterion to reverse some of the sequences, since many of the
#matches are just because a given sequence is present in both orientations in the database.

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

import sys #We want to be able to output error messages to the screen.
from collections import defaultdict
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

Usage = '''
blast_parse_template.py version 3.0, 30 April 2014
This script parses the output from the BLAST alignment of the template
sequences to determine which template sequences are duplicates.
Input is of the form:
blast_parse_template.py [output filename for unique sequences] 
[output folder for sets of matching sequences] [BLAST output filename--templates] 
   [BLAST output filename--primers] [sequence filename--fasta]
'''

if len(sys.argv) < 4:
	print Usage
else:
	OutFileUniq = sys.argv[1]
	OutFilePre = sys.argv[2]
	BLFileName = sys.argv[3]
	BLPFileName = sys.argv[4]
	SeqFileName = sys.argv[5]

print("blast_parse_template_4.py")

#The dictionaries we will fill out as we read the input files:
TSeqDict = { } #The dictionary of template sequences
PrimerList = [ ] #The list of sequences that still include primer fragments
BLDict = defaultdict(dict) #The dictionary of all of the information from the BLAST file

#The dictionaries and lists we will make when we analyze the data:
UniqSeqs = [ ] #The list of unique sequences
NonUniq = { } #Dictionary of non-unique sequences and the groups in which they fall
GroupDict = { } #Dictionary of the groups and the sequences they contain
ScriptList = [ ] #The shell script that will align the sequences within groups.
NonUniqALList = [ ]
RevComp = 0 #The number of sequences that are duplicates because they are just the
#reverses of other sequences we had.

#First, make a dictionary of all the sequences from the file of sequences.
InFile = open(SeqFileName, 'rU')
LineNum = 0
for Line in InFile:
	Line = Line.strip('\n').strip('\r')
	#Each sequence is two lines, first the name, with a > symbol first,
	if Line[0] == ">":
		SeqName = Line[1:]
	#then the sequence itself.
	else:
		Seq = Line
		TSeqDict[SeqName] = Seq
	LineNum += 1
InFile.close()

if LineNum/2 != len(TSeqDict):
	print ("Error, %d sequences were read, but only %d were put in the dictionary.  Perhaps there are duplicates.\n" % (LineNum/2, len(TSeqDict)))
	sys.stderr.write("Error, %d sequences were read, but only %d were put in the dictionary.  Perhaps there are duplicates.\n" % (LineNum/2, len(TSeqDict)))

print ("%d sequences were read from the input file %s.\n" % (LineNum/2, SeqFileName))
sys.stderr.write("%d sequences were read from the input file %s.\n" % (LineNum/2, SeqFileName))

#Now we need to see which sequences match the primers, so those can be removed.
InFile = open(BLPFileName, 'rU')
LineNum = 0
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	SSeq = Line[1]
	#if the portion that matches the primer is more than 20 bp
	if (int(Line[3]) >= 20):
		PrimerList.append(SSeq)
InFile.close()
PrimerList = list(set(PrimerList))

print ("%d of these sequences must be ignored because they still include primer sequences, \
other sequences have microsatellites and are also not picked up by BLAST.\n" % \
(len(PrimerList)))
sys.stderr.write("%d of these sequences must be ignored because they still include primer sequences, \
other sequences have microsatellites and are also not picked up by BLAST.\n" % \
(len(PrimerList)))

SeqList = [ ]
#Now read the output from BLASTing the template sequences against one another:
OldSeq = ""
InFile = open(BLFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	Qseq = Line[0]
	Sseq = Line[1]
	SeqList.append(Qseq)
	if (Qseq in PrimerList) == False:
		#If this is a new sequence, start the count for matching sequences.
		if Qseq != OldSeq:
			BLDict[Qseq]['NumMatches'] = 0
			OldSeq = Qseq
		#If the subject and the query sequence are not the same but are sufficiently
		#similar to one another:
		if Qseq != Sseq and (Sseq in PrimerList) == False and ((float(Line[2]) >= 90.0) and (int(Line[3]) > 30)):
			BLDict[Qseq]['NumMatches'] += 1
			BLDict[Qseq][Sseq] = defaultdict(dict)
			BLDict[Qseq][Sseq]['pident'] = float(Line[2])
			BLDict[Qseq][Sseq]['pmatching'] = float(Line[3])/float(Line[12])
			if (int(Line[8]) < int(Line[9])):
				BLDict[Qseq][Sseq]['ForR'] = 'F'
			else:
				BLDict[Qseq][Sseq]['ForR'] = 'R'
InFile.close()

print ("Information about %d sequences was read from infile %s.\n" % (len(BLDict), BLFileName))
sys.stderr.write("Information about %d sequences was read from infile %s.\n" % (len(BLDict), BLFileName))

#Now to go through the dictionary of sequences to determine what to do with them.
NumUniq = 0
for Qseq in BLDict:
	if BLDict[Qseq]['NumMatches'] == 0:
		UniqSeqs.append(Qseq)
		NumUniq += 1
	elif BLDict[Qseq]['NumMatches'] == 1:
		for Sseq in BLDict[Qseq]:
			if ((Sseq != 'NumMatches') and (Sseq in UniqSeqs) == False):
				#print ("%s  %s  %f  %f  %s" % (Qseq, Sseq, BLDict[Qseq][Sseq]['pident'], BLDict[Qseq][Sseq]['pmatching'], BLDict[Qseq][Sseq]['ForR'])) 
				if (BLDict[Qseq][Sseq]['pident'] > 0.9) and \
				(BLDict[Qseq][Sseq]['pmatching'] > 0.98) and \
				(BLDict[Qseq][Sseq]['ForR'] == 'R'):
					#Add the sequence to the list of unique sequences
					UniqSeqs.append(Qseq)
					#And remove the match from BLDict
					RevComp += 1
				else:
					try:
						NonUniq[Qseq] = NonUniq[Sseq]
						GroupDict[NonUniq[Sseq]].append(Qseq)
					except KeyError:
						GroupDict[Qseq] = [Qseq]
						NonUniq[Qseq] = Qseq
	elif BLDict[Qseq]['NumMatches'] > 1:
		GroupListTemp = [ ]
		#add the sequences to the lists of sequences with matches
		for Sseq in BLDict[Qseq]:
			if Sseq != 'NumMatches':
				try:
					GroupListTemp.append(NonUniq[Sseq])
				except KeyError:
					"Do nothing"
		GroupListTemp = list(set(GroupListTemp))
		if len(GroupListTemp) == 0:
			#Make a new group named after that sequence
			GroupDict[Qseq] = [Qseq]
			NonUniq[Qseq] = Qseq
		#If the sequence only belongs to one group,
		elif len(GroupListTemp) == 1:
			#give it that group
			NonUniq[Qseq] = GroupListTemp[0]
			#and add it to that group's list of sequences
			GroupDict[GroupListTemp[0]].append(Qseq)
		elif len(GroupListTemp) > 1:
			#The new name for the group is the first group name on the list
			GroupName = GroupListTemp[0]
			#Go through the other group names,
			for GNameTemp in GroupListTemp:
				if GNameTemp != GroupName:
					#add the members of those groups to the first group
					GroupDict[GroupName] += GroupDict[GNameTemp]
					#and delete the separate entries for those groups
					del GroupDict[GNameTemp]
					#Also give all of those sequences the new group
					for SeqName in GroupDict[GroupName]:
						NonUniq[SeqName] = GroupName
			#and finally add the sequence to that group
			NonUniq[Qseq] = GroupName
			GroupDict[GroupName].append(Qseq)

if (NumUniq + RevComp != len(UniqSeqs)):
	sys.stderr.write("ERROR!!!  %d unique sequences and %d sequences with only a \n\
	reverse-complement were found, but %d were written to the list of sequences.\n" % \
	(NumUniq, RevComp, len(UniqSeqs)))
				
print("%d sequences did not match other sequences closely, while %d only had one match, which was their reverse-complement.\n" \
	 % (NumUniq, RevComp))
print("%d sequences in %d groups are more complicated and need to be examined further.\n" \
	% (len(NonUniq), len(GroupDict)))
sys.stderr.write("%d sequences did not match other sequences closely, while %d only had one match, which was their reverse-complement.\n" \
	 % (NumUniq, RevComp))
sys.stderr.write("%d sequences in %d groups are more complicated and need to be examined further.\n" \
	% (len(NonUniq), len(GroupDict)))

OutFile = open(OutFileUniq, 'w')
for SeqName in UniqSeqs:
	Line = ">"+SeqName+"\n"+TSeqDict[SeqName]+"\n"
	OutFile.write(Line)
OutFile.close()

print("%d unique sequences were written to the outfile %s.\n" % (len(UniqSeqs), OutFileUniq))
sys.stderr.write("%d unique sequences were written to the outfile %s.\n" % (len(UniqSeqs), OutFileUniq))

Line = "#! /bin/bash\n\n"
ScriptList.append(Line)

for GroupName in GroupDict:
	OutList = [ ]
	for SeqName in GroupDict[GroupName]:
		Line = ">"+SeqName+"\n"+TSeqDict[SeqName]+"\n"
		OutList.append(Line)
		Line = GroupName+"\t"+SeqName+"\n"
		NonUniqALList.append(Line)
	OutFileName = OutFilePre+GroupName+".fa"
	OutFile = open(OutFileName, 'w')
	for Line in OutList:
		OutFile.write(Line)
	OutFile.close()
	Line = "muscle3.8.31_i86linux64 -in "+OutFileName+" -out "+OutFilePre+"a"+GroupName+".fa\n"
	ScriptList.append(Line)

print ("%d groups of matching sequences were written, with file names similar to %s.\n" % (len(GroupDict), OutFileName))
sys.stderr.write("%d groups of matching sequences were written, with file names similar to %s.\n" % (len(GroupDict), OutFileName))

OutFileName = OutFilePre+"Seq_List.txt"
OutFile = open(OutFileName, 'w')
for Line in NonUniqALList:
	OutFile.write(Line)
OutFile.close()

print("A list of the groups and the sequences they contain has been written to the file %s.\n" % (OutFileName))
print("This list can be edited so that it contains only the sequences of interest\n\
	and can then be parsed by the next script.\n")
sys.stderr.write("A list of the groups and the sequences they contain has been written to the file %s.\n" % (OutFileName))
sys.stderr.write("This list can be edited so that it contains only the sequences of interest\n\
	and can then be parsed by the next script.\n")

OutFileName = OutFilePre+"Alignment_Script.sh"
OutFile = open(OutFileName, 'w')
for Line in ScriptList:
	OutFile.write(Line)
OutFile.close()

print("The script to align these files in Muscle is %s.\n" % (OutFileName))
print("You will need to give yourself permission to execute this file (chmod u+x) before you can execute it.\n")
sys.stderr.write("The script to align these files in Muscle is %s.\n" % (OutFileName))
sys.stderr.write("You will need to give yourself permission to execute this file (chmod u+x) before you can execute it.\n")

