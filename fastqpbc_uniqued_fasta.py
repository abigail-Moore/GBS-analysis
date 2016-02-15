#! /usr/bin/env python

''' "Uniqued" file format (with sequence and Phred score trimmed to 5 bp):

GATCC	9	EEHFG	Mlk38_20,Mlk38_19	7,2	.	.

GATCC	44	DDFGF	Mlk38_20	44	.	.

There does not appear to be anything special at the end of the file.'''

#This program reads the files of with the sequences sorted by individual by fastqcomb_pbc.py.
#It then writes a fasta file with the unique sequences.  These sequences are labeled as follows:
#sequence number (sequential number from 1 to the number of sequences)_individual name_number of sequences in that individual
#_2nd individual's name_number of sequences in 2nd individual_...

import re #We will need the regular expressions module
import sys #We want to be able to output error messages to the screen.
import gzip #We want to be able to open zipped files.

from collections import defaultdict #We need defaultdict to be able to create dictionaries with multiple levels.

#First, get the files names from the command line or print some instructions if not enough names are found..

Usage = '''
fastqpbc_uniqued_fasta.py - version 1.0, 13 March 2014
(from fastqpbc_unqiued_3.py)
Reads separate files of the sequences present in each individual (output from fastqcomb_pbc.py)
and combines them into a fasta file that can be aligned using blat or blast to the template sequences.
Usage:
	fastqpbc_uniqued_fasta.py (output file name--without suffix) 
	(sequence file names)
The file names can be either the full path, e.g., /home/moorea/... 
or the file name only if you are in that directory.

Sequence file names can either be a list of all sequence file names,
or something of the form "outfile_*.fastq.gz".
And the sequence files can either be be zipped (with a ".gz" ending) or not.
'''

if len(sys.argv) <3:
	print Usage
else:
	OutFileRoot = sys.argv[1]
	InFileList = sys.argv[2:]

print ("fastqpbc_uniqued_fasta.py")
	

#Files we will write and dictionaries that will hold the information
OutFileNameFa = OutFileRoot+"_uniq.fa"

UniquedSeqs = defaultdict(dict) #dictionary for the uniqued fasta sequences
NumInFiles = 0
NumSeqsTot = 0

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
		elif SeqLine == 1: #The sequence itself is read.
			Seq = Line
			try: #First, I assume the sequence is already in the list, and just add that this individual to the sequence's entry.
				UniquedSeqs[Seq][Ind] += 1
			except KeyError: #If the individual not found, then I add an entry for that individual
				UniquedSeqs[Seq][Ind] = 1
			LineNum += 1
			NumSeqs += 1
		else: #All other lines of the sequence (position, quality, and +) are ignored.
			LineNum += 1
	sys.stderr.write("%d sequences have been read from the file %s.\n" % (NumSeqs, InFileName))
	print ("%d sequences have been read from the file %s.\n" % (NumSeqs, InFileName))
	InFile.close()
	NumInFiles += 1
	NumSeqsTot += NumSeqs

sys.stderr.write("%d sequences were read from %d files.\n" % (NumSeqsTot, NumInFiles))
print("%d sequences were read from %d files.\n" % (NumSeqsTot, NumInFiles))

if len(InFileList) != NumInFiles:
	sys.stderr.write("Error: %d files were in the list of infiles, but only %d files were read.\n" % (len(InFileList), NumInFiles))
	print("Error: %d files were in the list of infiles, but only %d files were read.\n" % (len(InFileList), NumInFiles))

#%%%%%%%%%%%%Reformatting the dictionary and writing it to a file.%%%%%%%%%%%%%%%%%%%%%%%%

#The uniqued sequences are written to a fasta file for alignment with blat:
FinalSeqs = [ ]
SeqNum = 0
for Seq in UniquedSeqs.keys():
	SeqNum += 1
	NameList = [ ]
	for Ind in UniquedSeqs[Seq].keys():
		Name = Ind+"_"+str(UniquedSeqs[Seq][Ind])
		NameList.append(Name)
	Line = ">"+str(SeqNum)+"_"+"_".join(NameList)+"\n"+Seq+"\n"
	FinalSeqs.append(Line)

OutFile = open(OutFileNameFa, 'w')
for Line in FinalSeqs:
	OutFile.write(Line)
OutFile.close()

if SeqNum != len(FinalSeqs):
	print ("Error!  %d sequences were found, but only %d were written to the file of sequences!\n" % (SeqNum, len(FinalSeqs)))
	sys.stderr.write("Error!  %d sequences were found, but only %d were written to the file of sequences!\n" % (SeqNum, len(FinalSeqs)))

print ("%d unique sequences were written to the file %s.\n" % (SeqNum, OutFileNameFa))
sys.stderr.write("%d unique sequences were written to the file %s.\n" % (SeqNum, OutFileNameFa))
